#!/usr/bin/env python3
"""
Nightly pipeline:
1) Download a public copy of the source Google Sheet as CSV.
2) Run PubMed searches for each placeholder cell.
3) Emit:
   - A computed results CSV
   - A JSON + GitHub Pages HTML view with search links
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from typing import List, Optional
from urllib.parse import urlencode
import urllib.request
import xml.etree.ElementTree as ET

import datetime
import html as htmllib


DEFAULT_SOURCE_SHEET_ID = "1VidNfYpvIzB7SA3SePyhHUil0j-XP1TtiakLfWl24pg"
DEFAULT_SOURCE_GID = "363407775"
PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
PUBMED_WEB_URL = "https://pubmed.ncbi.nlm.nih.gov/"
HEATMAP_LOW_COLOR = (250, 252, 255)  # near-white
HEATMAP_HIGH_COLOR = (15, 56, 149)  # strong indigo


PLACEHOLDER_VALUES = {"", "Loading...", "loading"}
NON_RESULT_VALUES = {"#N/A", "N/A", "-", "NA"}


@dataclass(frozen=True)
class Config:
    source_sheet_id: str
    source_gid: str
    data_dir: Path
    public_dir: Path
    max_rows: int
    request_delay: float
    api_key: Optional[str]
    summarize: bool
    summary_model: str
    summary_max_abstracts: int
    summary_top_cells: int
    summary_api_key: Optional[str]
    summary_api_base: str


def _coerce_path(value: str) -> Path:
    return Path(value).expanduser().resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-sheet-id", default=DEFAULT_SOURCE_SHEET_ID)
    parser.add_argument("--source-gid", default=DEFAULT_SOURCE_GID)
    parser.add_argument("--output-dir", default="data")
    parser.add_argument("--public-dir", default="public")
    parser.add_argument("--max-rows", type=int, default=0, help="Optional row limit for quick runs in tests")
    parser.add_argument("--request-delay", type=float, default=float(os.getenv("PUBMED_REQUEST_DELAY", "0.34")))
    parser.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    parser.add_argument("--force-refresh", action="store_true", help="Ignore cached values and recompute all query cells")
    parser.add_argument(
        "--summarize",
        action="store_true",
        help="Download abstracts and generate summaries for searchable results",
    )
    parser.add_argument(
        "--summary-model",
        default=os.getenv("OPENAI_MODEL", "gpt-4o-mini"),
        help="LLM model for summaries (used when --summarize and OPENAI_API_KEY are present)",
    )
    parser.add_argument(
        "--summary-max-abstracts",
        type=int,
        default=int(os.getenv("PUBMED_SUMMARY_MAX_ABSTRACTS", "8")),
        help="Maximum abstracts fetched per summarized query",
    )
    parser.add_argument(
        "--summary-top-cells",
        type=int,
        default=int(os.getenv("PUBMED_SUMMARY_TOP_CELLS", "40")),
        help="Maximum number of cell summaries generated per run",
    )
    return parser.parse_args()


def sheet_csv_url(sheet_id: str, gid: str) -> str:
    return (
        "https://docs.google.com/spreadsheets/d/"
        f"{sheet_id}/export?format=csv&gid={gid}"
    )


def fetch_csv_rows(url: str) -> List[List[str]]:
    """Download the public sheet export and parse it as rows."""
    with urllib.request.urlopen(url, timeout=60) as resp:
        text = resp.read().decode("utf-8-sig", errors="replace")
    reader = csv.reader(text.splitlines())
    return [row for row in reader]


def normalize_rows(rows: List[List[str]]) -> List[List[str]]:
    max_cols = max((len(r) for r in rows), default=0)
    normalized = []
    for row in rows:
        if len(row) < max_cols:
            row = row + ["" for _ in range(max_cols - len(row))]
        normalized.append(row)
    return normalized


def read_local_csv(path: Path) -> List[List[str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as f:
        return [row for row in csv.reader(f)]


def read_local_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    return payload if isinstance(payload, dict) else {}


def apply_cached_results(rows: List[List[str]], cached_rows: List[List[str]]) -> List[List[str]]:
    """
    Pre-fill row values from a previously generated results file where source cells are placeholders.
    """
    if not rows or not cached_rows or len(rows) < 3 or len(cached_rows) < 3:
        return rows

    source_headers = rows[0]
    cached_headers = cached_rows[0]
    if not source_headers or not cached_headers:
        return rows

    cached_header_index: dict[str, int] = {}
    for idx, header in enumerate(cached_headers):
        header_value = (header or "").strip()
        if header_value:
            cached_header_index[header_value.lower()] = idx

    def row_key(row: List[str]) -> str:
        return (row[0] or "").strip().lower() or "__compound-axis__"

    cached_by_key: dict[str, List[str]] = {}
    for cached_row in cached_rows[2:]:
        cached_by_key[row_key(cached_row)] = cached_row

    for row in rows[2:]:
        cached_row = cached_by_key.get(row_key(row))
        if not cached_row:
            continue

        max_col = min(len(source_headers), len(cached_row), len(cached_headers))
        if len(row) < max_col:
            row.extend([""] * (max_col - len(row)))

        for col in range(2, max_col):
            current_value = (row[col] or "").strip()
            if current_value not in PLACEHOLDER_VALUES:
                continue
            header = (source_headers[col] or "").strip()
            if not header:
                continue
            cached_col = cached_header_index.get(header.lower())
            if cached_col is None or cached_col >= len(cached_row):
                continue
            cached_value = (cached_row[cached_col] or "").strip()
            if cached_value and cached_value not in PLACEHOLDER_VALUES:
                row[col] = cached_value

    return rows


def build_query(term: str, compound: str) -> str:
    """Keep terms stable and mirror sheet-like text queries in PubMed."""
    safe_term = f'"{term}"[All Fields]'
    safe_comp = f'"{compound}"[All Fields]'
    return f"{safe_term} AND {safe_comp}"


def build_term_only_query(term: str) -> str:
    """Query for a term by itself."""
    return f'"{term}"[All Fields]'


def build_compound_only_query(compound: str) -> str:
    """Query for a compound by itself."""
    return f'"{compound}"[All Fields]'


def _strip_and_normalize_term(value: str) -> str:
    value = (value or "").strip()
    if not value:
        return ""
    value = value.replace("\u2013", "-").replace("\u2014", "-")
    value = value.replace("\u2019", "'").replace("\u201c", '"').replace("\u201d", '"')
    value = re.sub(r"\s+", " ", value).strip()
    return value


def _query_term_variants(value: str) -> List[str]:
    """Generate a small set of broader query term variants."""
    base = _strip_and_normalize_term(value)
    if not base:
        return []

    variants: List[str] = []
    seen = set()

    def add(candidate: str) -> None:
        candidate = _strip_and_normalize_term(candidate)
        if not candidate or candidate in seen:
            return
        seen.add(candidate)
        variants.append(candidate)

    add(base)

    # Remove common wrappers and punctuation clutter.
    add(re.sub(r'["\'`]+', "", base))
    add(re.sub(r"[()\[\]{}]", "", base))
    add(re.sub(r"[;,]|//", " ", base))
    add(base.replace("/", " "))
    add(base.replace("-", " "))

    # Try removing parenthetical qualifiers, which often block strict term hits.
    no_parens = re.sub(r"\s*\([^)]*\)\s*", " ", base)
    add(no_parens)

    words = base.split()
    if len(words) >= 3:
        add(" ".join(words[:2]))
        add(" ".join(words[-2:]))

    # Remove a few generic trailing qualifiers that can overconstrain phrases.
    for suffix in ("syndrome", "condition", "disease"):
        lower_words = [w.lower() for w in words]
        if lower_words and lower_words[-1] == suffix:
            add(" ".join(words[:-1]))

    return variants[:4]


def _build_fielded_query(text: str, quoted: bool = True) -> str:
    field = text.replace('"', "").strip()
    if not field:
        return ""
    if quoted:
        return f'"{field}"[All Fields]'
    if re.search(r"[^0-9A-Za-z\\s-]", field):
        return ""
    return f"{field}[All Fields]"


def build_term_query_variants(term: str) -> List[str]:
    """Build alternate query candidates for a single term."""
    candidates = [_build_fielded_query(variant, quoted=True) for variant in _query_term_variants(term)]
    candidates.extend(_build_fielded_query(variant, quoted=False) for variant in _query_term_variants(term))
    return list(dict.fromkeys([c for c in candidates if c]))


def build_query_variants(term: str, compound: str) -> List[str]:
    """Build a small ordered set of fallback pairwise query variants."""
    term_candidates = _query_term_variants(term)
    compound_candidates = _query_term_variants(compound)
    candidates: List[str] = []

    # Keep strictest behavior first, then broaden progressively.
    combos = (
        (True, True),
        (False, True),
        (True, False),
        (False, False),
    )
    for term_variant in term_candidates[:3]:
        term_quoted = _build_fielded_query(term_variant, quoted=True)
        term_unquoted = _build_fielded_query(term_variant, quoted=False)
        for compound_variant in compound_candidates[:2]:
            compound_quoted = _build_fielded_query(compound_variant, quoted=True)
            compound_unquoted = _build_fielded_query(compound_variant, quoted=False)
            if not compound_quoted and not compound_unquoted:
                continue
            for quoted_term, quoted_comp in combos:
                left = term_quoted if quoted_term else term_unquoted
                right = compound_quoted if quoted_comp else compound_unquoted
                if not left or not right:
                    continue
                candidates.append(f"{left} AND {right}")

    return list(dict.fromkeys(candidates))


def _api_request_json(url: str, params: dict, timeout: int = 30) -> dict | str:
    if api_key := params.pop("_api_key", None):
        params["api_key"] = api_key
    if "api_key" in params and not params["api_key"]:
        params.pop("api_key", None)
    full_url = f"{url}?{urlencode(params)}"
    with urllib.request.urlopen(full_url, timeout=timeout) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    return raw if params.get("retmode") == "xml" or params.get("retmode") == "text" else json.loads(raw)


def pubmed_search_count(query: str, api_key: Optional[str]) -> str:
    payload = _esearch_payload(query, api_key, retmax=0)
    return payload["esearchresult"]["count"]


def _esearch_payload(query: str, api_key: Optional[str], retmax: int = 0) -> dict:
    params = {
        "db": "pubmed",
        "retmode": "json",
        "retmax": retmax,
        "term": query,
    }
    return _api_request_json(PUBMED_ESEARCH_URL, params)  # type: ignore[return-value]


def pubmed_search_count_and_ids(
    query: str,
    api_key: Optional[str],
    retmax: int,
) -> tuple[str, List[str]]:
    payload = _esearch_payload(query, api_key, retmax=retmax)
    esearch = payload.get("esearchresult", {})
    count = str(esearch.get("count", "0"))
    ids = list(esearch.get("idlist", []))
    return count, ids


def fetch_pubmed_abstracts(pmids: Iterable[str], api_key: Optional[str]) -> List[dict]:
    ids = [pid.strip() for pid in pmids if str(pid).strip()]
    if not ids:
        return []
    chunks = [",".join(ids)]
    abstracts: List[dict] = []
    for chunk in chunks:
        params = {
            "db": "pubmed",
            "retmode": "xml",
            "rettype": "abstract",
            "id": chunk,
            "_api_key": api_key,
        }
        payload = _api_request_json(PUBMED_EFETCH_URL, params)  # type: ignore[assignment]
        abstracts.extend(_extract_abstracts_from_xml(payload))
    return abstracts


def _extract_abstracts_from_xml(xml_text: str) -> List[dict]:
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return []

    papers: List[dict] = []
    for article in root.findall(".//PubmedArticle"):
        pmid = ""
        pmid_node = article.find(".//PMID")
        if pmid_node is not None and pmid_node.text:
            pmid = pmid_node.text.strip()

        title_node = article.find(".//ArticleTitle")
        title = "".join(title_node.itertext()).strip() if title_node is not None else ""

        abstract_parts = []
        for node in article.findall(".//Abstract/AbstractText"):
            heading = (node.attrib.get("Label") or "").strip()
            text = "".join(node.itertext()).strip()
            if not text:
                continue
            if heading:
                abstract_parts.append(f"{heading}: {text}")
            else:
                abstract_parts.append(text)
        abstract = "\n".join(abstract_parts).strip()
        if not abstract:
            continue
        papers.append({"pmid": pmid, "title": title, "abstract": abstract})
    return papers


def _normalize_summary_payload(payload: dict) -> str:
    choices = payload.get("choices") or []
    if not choices:
        return ""
    message = choices[0].get("message") or {}
    content = message.get("content") or ""
    if isinstance(content, str):
        return content.strip()
    return ""


def generate_query_summary(
    query: str,
    abstracts: List[dict],
    api_key: Optional[str],
    base_url: str,
    model: str,
) -> str:
    if not api_key:
        return "Summary unavailable (OPENAI_API_KEY not set)."
    if not abstracts:
        return "No abstracts available for summarization."
    if base_url.endswith("/"):
        base_url = base_url[:-1]

    snippets: List[str] = []
    for item in abstracts:
        title = item.get("title") or "Unknown title"
        abstract = item.get("abstract") or ""
        if abstract:
            snippets.append(f"Title: {title}\nAbstract: {abstract}")
        if len(snippets) >= 4:
            break

    prompt = (
        "Summarize the evidence in the abstracts below for the PubMed query. "
        "Keep it short, evidence-focused, and cautious about inference.\n\n"
        f"Query: {query}\n\n"
    )
    prompt += "\n\n---\n\n".join(snippets)

    body = {
        "model": model,
        "temperature": 0.2,
        "max_tokens": 240,
        "messages": [
            {
                "role": "system",
                "content": (
                    "You are a scientific assistant. Summarize biomedical findings briefly and clearly. "
                    "Mention the type of evidence represented by the abstracts and avoid overclaiming."
                ),
            },
            {"role": "user", "content": prompt},
        ],
    }
    url = f"{base_url}/chat/completions"
    req = urllib.request.Request(
        url,
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as resp:
        response_raw = resp.read().decode("utf-8", errors="replace")
    response = json.loads(response_raw)
    return _normalize_summary_payload(response)


def load_summary_cache(payload: dict) -> dict[str, str]:
    rows = payload.get("rows")
    if not isinstance(rows, list):
        return {}

    cache: dict[str, str] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        for cell in row.get("cells", []) if isinstance(row.get("cells"), list) else []:
            if not isinstance(cell, dict):
                continue
            query = (cell.get("query") or "").strip()
            summary = (cell.get("summary") or "").strip()
            if query and summary:
                cache[query] = summary
    return cache


def fetch_with_retry(query: str, api_key: Optional[str], delay: float, retries: int = 3) -> str:
    attempt = 0
    while True:
        try:
            count = pubmed_search_count(query, api_key)
            return count
        except Exception as err:  # pragma: no cover - operational retry logic
            attempt += 1
            if attempt > retries:
                return f"ERROR: {err.__class__.__name__}"
        time.sleep(min(2**attempt * delay, 5.0))


def fetch_with_zero_fallback(
    query_variants: List[str],
    api_key: Optional[str],
    delay: float,
) -> tuple[str, str]:
    """Try each query until a non-zero count is found; return count + winning query."""
    if not query_variants:
        return "0", ""
    last_count = ""
    for candidate in query_variants:
        count = fetch_with_retry(candidate, api_key, delay)
        if not str(count).startswith("ERROR:"):
            last_count = count
        if count != "0" and not str(count).startswith("ERROR:"):
            return count, candidate
    return last_count, query_variants[0]


def is_search_cell(value: str) -> bool:
    v = (value or "").strip()
    if v in NON_RESULT_VALUES:
        return False
    if v in PLACEHOLDER_VALUES:
        return True
    return False


def make_pubmed_url(query: str) -> str:
    return f"{PUBMED_WEB_URL}?{urlencode({'term': query})}"


def parse_count(value: str) -> Optional[int]:
    try:
        v = str(value).replace(",", "").strip()
        return int(v)
    except ValueError:
        return None


def compact_count(value: str) -> str:
    count = parse_count(value)
    if count is None or count < 0:
        return value
    if count >= 10000:
        return f"{count // 1000}k"
    return str(count)


def compute_heat_bounds(rows: List[List[str]], headers: List[str], start_col: int = 3) -> tuple[int, int]:
    values: List[int] = []
    for row in rows:
        for col in range(start_col, len(headers)):
            if col < len(row):
                count = parse_count(row[col])
                if count is not None and count >= 0:
                    values.append(count)
    if not values:
        return 0, 0
    return min(values), max(values)


def summarize_pairwise_cells(rows: List[List[str]], headers: List[str], max_items: int = 15) -> tuple[dict[tuple[int, int], str], str]:
    candidates = []

    for row_idx, row in enumerate(rows[2:], start=2):
        term = (row[0] or "").strip()
        is_compound_only_row = row_idx == 2 and not term
        if not term and not is_compound_only_row:
            continue
        if is_compound_only_row:
            continue

        row_term = term or "compound-only"
        row_context = (row[1] or "").strip()
        for col in range(3, len(headers)):
            compound = (headers[col] or "").strip()
            if not compound:
                continue
            raw = row[col] if col < len(row) else ""
            count = parse_count(raw)
            if count is None:
                continue
            query = build_query(row_term, compound)
            candidates.append(
                {
                    "row_idx": row_idx,
                    "col": col,
                    "count": count,
                    "query": query,
                    "term": row_term,
                    "context": row_context,
                    "compound": compound,
                }
            )

    if not candidates:
        return {}, ""

    ranked = sorted((c for c in candidates if c["count"] > 0), key=lambda item: item["count"], reverse=True)
    if not ranked:
        return {}, ""

    top_count = min(max_items, len(ranked))
    cutoff = ranked[top_count - 1]["count"]
    highlights = [c for c in ranked if c["count"] >= cutoff]

    summary_rows = []
    for item in highlights:
        term_label = item["term"]
        if item["context"]:
            term_label = f"{term_label} ({item['context']})"
        summary_rows.append(
            {
                "term": term_label,
                "compound": item["compound"],
                "count": item["count"],
                "query": item["query"],
                "row_idx": item["row_idx"],
                "col": item["col"],
            }
        )

    summary_html = []
    for item in summary_rows:
        summary_html.append(
            "<li>"
            f"{htmllib.escape(item['term'])} × {htmllib.escape(item['compound'])}: "
            f"{item['count']:,} hits (query: {htmllib.escape(item['query'])})"
            "</li>"
        )

    summary_map: dict[tuple[int, int], str] = {}
    for item in summary_rows:
        summary_map[(item["row_idx"], item["col"])] = (
            f"{item['term']} × {item['compound']}: {item['count']:,} hits.\n"
            f"PubMed query: {item['query']}"
        )

    return summary_map, "".join(summary_html)


def heatstyle(count: int, min_count: int, max_count: int) -> tuple[str, str]:
    if max_count <= 0:
        return "", "#0f172a"
    if count <= 0:
        return "", "#0f172a"
    if max_count == min_count:
        ratio = 1.0
    else:
        min_log = math.log1p(min_count)
        max_log = math.log1p(max_count)
        count_log = math.log1p(count)
        if max_log == min_log:
            ratio = 0.0
        else:
            ratio = (count_log - min_log) / (max_log - min_log)
    ratio = min(max(ratio, 0.0), 1.0)
    r = round(HEATMAP_LOW_COLOR[0] + (HEATMAP_HIGH_COLOR[0] - HEATMAP_LOW_COLOR[0]) * ratio)
    g = round(HEATMAP_LOW_COLOR[1] + (HEATMAP_HIGH_COLOR[1] - HEATMAP_LOW_COLOR[1]) * ratio)
    b = round(HEATMAP_LOW_COLOR[2] + (HEATMAP_HIGH_COLOR[2] - HEATMAP_LOW_COLOR[2]) * ratio)
    # Switch to white text when the background is dark enough.
    text_color = "#ffffff" if ratio >= 0.55 else "#111827"
    return f"background-color: rgb({r},{g},{b}); color: {text_color};", text_color


def run_search_grid(
    rows: List[List[str]],
    config: Config,
    force_refresh: bool = False,
    summary_cache: Optional[dict[str, str]] = None,
) -> tuple[List[List[str]], List[List[dict]], str]:
    """
    Return:
      - updated rows with counts
      - list of row-wise dicts for JSON export
      - generated timestamp label
    """
    timestamp = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat(timespec="seconds")
    headers = rows[0] if rows else []
    data_rows = rows[2:]  # keeps row format from observed sheet layout and preserves row 2 axis entries
    result_rows = [row[:] for row in rows]
    json_rows = []
    if len(headers) < 3:
        headers.extend([""] * (3 - len(headers)))
    headers[2] = "term-only"

    if not data_rows:
        return result_rows, [], timestamp

    summary_cache = dict(summary_cache or {})
    summary_budget = [max(0, config.summary_top_cells)]

    def make_summary(query: str, count_value: str, source: str) -> str:
        if summary_budget[0] <= 0 or not config.summarize:
            return ""
        if not config.summary_api_key:
            summary_budget[0] -= 1
            return "Summary unavailable (OPENAI_API_KEY not set)."

        count = parse_count(count_value)
        if count is None or count <= 0:
            return ""

        if query in summary_cache:
            return summary_cache[query]

        try:
            _, pmids = pubmed_search_count_and_ids(
                query,
                config.api_key,
                min(config.summary_max_abstracts, max(1, count)),
            )
            if not pmids:
                return f"No abstracts available for query: {query}"
            abstracts = fetch_pubmed_abstracts(pmids, config.api_key)
            if not abstracts:
                return f"No abstract text available for query: {query}"
            summary = generate_query_summary(
                query,
                abstracts,
                config.summary_api_key,
                config.summary_api_base,
                config.summary_model,
            )
            if not summary:
                return f"Summary unavailable for query: {query}"
            summary_cache[query] = summary
            summary_budget[0] -= 1
            return summary
        except Exception as err:  # pragma: no cover - resilience around network/API errors
            return f"Summary unavailable ({err.__class__.__name__})"

    data_rows_to_process = data_rows[: config.max_rows] if config.max_rows > 0 else data_rows
    for data_index, source_row in enumerate(data_rows_to_process, start=2):
        term = (source_row[0] or "").strip()
        is_compound_axis_row = data_index == 2 and not term
        if not term and not is_compound_axis_row:
            continue

        term_label = term or "compound-only"
        json_cells = []
        row_out = result_rows[data_index]

        if len(row_out) <= 2:
            row_out.extend([""] * (3 - len(row_out)))

        if not is_compound_axis_row:
            term_query_candidates = build_term_query_variants(term)
            base_cell_value = (row_out[2] if len(row_out) > 2 else "").strip()
            is_base_to_compute = force_refresh or is_search_cell(base_cell_value)
            base_query = build_term_only_query(term)
            if is_base_to_compute and term_query_candidates:
                base_count, base_query = fetch_with_zero_fallback(
                    term_query_candidates,
                    config.api_key,
                    config.request_delay,
                )
                row_out[2] = base_count
                base_source = "queried" if base_query == term_query_candidates[0] else "fallback"
            else:
                base_count = base_cell_value
                base_source = "kept"

            base_summary = make_summary(base_query, str(base_count), base_source)
            json_cells.append(
                {
                    "column": headers[2],
                    "value": base_count,
                    "query": base_query,
                    "link": make_pubmed_url(base_query),
                    "source": base_source,
                    "summary": base_summary,
                }
            )
            if is_base_to_compute:
                time.sleep(config.request_delay)
        else:
            row_out[2] = row_out[2] if len(row_out) > 2 else ""

        for col_idx in range(3, len(headers)):
            if col_idx >= len(row_out):
                row_out.extend([""] * (col_idx + 1 - len(row_out)))
            cell_value = (row_out[col_idx] or "").strip()

            compound = (headers[col_idx] or "").strip()
            if not compound:
                continue

            is_pair_to_compute = force_refresh or is_search_cell(cell_value)

            if is_compound_axis_row:
                compound_queries = build_term_query_variants(compound)
                query = compound_queries[0] if compound_queries else build_compound_only_query(compound)
                if is_pair_to_compute:
                    count, query = fetch_with_zero_fallback(
                        compound_queries if compound_queries else [query],
                        config.api_key,
                        config.request_delay,
                    )
                    row_out[col_idx] = count
                    time.sleep(config.request_delay)
                else:
                    count = cell_value
            elif not is_pair_to_compute:
                # Keep pre-existing values (including known N/A or 0 totals).
                query = build_query(term, compound)
                count = cell_value
                source = "kept"
            else:
                query_variants = build_query_variants(term, compound)
                count, query = fetch_with_zero_fallback(
                    query_variants,
                    config.api_key,
                    config.request_delay,
                )
                source = "queried" if query == query_variants[0] else "fallback"
                row_out[col_idx] = count
                time.sleep(config.request_delay)

            if is_compound_axis_row:
                source = "compound-only" if query == (compound_queries[0] if compound_queries else build_compound_only_query(compound)) else "fallback"

            summary = make_summary(query, str(count), source)
            json_cells.append(
                {
                    "column": headers[col_idx],
                    "value": count,
                    "query": query,
                    "link": make_pubmed_url(query),
                    "source": source,
                    "summary": summary,
                }
            )
            if not is_compound_axis_row:
                time.sleep(config.request_delay)

        json_rows.append({"term": term_label, "base": source_row[1] if len(source_row) > 1 else "", "cells": json_cells})

    return result_rows, json_rows, timestamp


def write_csv(path: Path, rows: List[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def render_html(public_dir: Path, headers: List[str], rows: List[List[str]], json_rows: List[dict], timestamp: str) -> None:
    public_dir.mkdir(parents=True, exist_ok=True)
    data_rows = rows[2:]
    heat_min, heat_max = compute_heat_bounds(data_rows, headers, start_col=2)
    has_heat_values = heat_max > 0
    summary_map, summary_html = summarize_pairwise_cells(rows, headers)
    json_summary_map: dict[tuple[int, int], str] = {}
    for row_idx, row_payload in enumerate(json_rows, start=2):
        for cell in row_payload.get("cells", []):
            if not isinstance(cell, dict):
                continue
            col_name = (cell.get("column") or "").strip()
            if not col_name:
                continue
            try:
                col = headers.index(col_name)
            except ValueError:
                continue
            summary = (cell.get("summary") or "").strip()
            if summary:
                json_summary_map[(row_idx, col)] = summary

    def render_cell(
        value: str,
        link: Optional[str],
        query: Optional[str] = None,
        is_heat_cell: bool = False,
        row_idx: Optional[int] = None,
        col: Optional[int] = None,
        summary: Optional[str] = None,
    ) -> str:
        display_value = compact_count(value)
        if value in NON_RESULT_VALUES or not value:
            return "<td></td>"
        if value.startswith("ERROR:"):
            return f"<td class=\"error\">{value}</td>"
        cell_style = ""
        link_style = ""
        cell_class = ""
        summary_text = summary
        is_highlight = False
        if is_heat_cell and row_idx is not None and col is not None:
            summary_text = json_summary_map.get((row_idx, col), summary_text)
            is_highlight = summary_text is not None
            if summary_text is None and query:
                count = parse_count(value)
                if count is None:
                    summary_text = f"Search: {query}"
                elif count == 0:
                    summary_text = f"No hits for query: {query}"
                else:
                    summary_text = f"{count:,} hits for query: {query}"
        if is_highlight and summary_text:
            cell_class = " class=\"high-score\""
        if is_heat_cell and has_heat_values:
            count = parse_count(value)
            if count is not None:
                style, text_color = heatstyle(count, heat_min, heat_max)
                if style:
                    cell_style = f" style=\"{style}\""
                    link_style = f" style=\"color:{text_color}\""
        if not link:
            if summary_text:
                return f"<td{cell_class} title=\"{htmllib.escape(summary_text)}\"{cell_style}>{display_value}</td>"
            return f"<td{cell_style}{cell_class}>{display_value}</td>"
        if summary_text:
            return f"<td{cell_class} title=\"{htmllib.escape(summary_text)}\"{cell_style}><a href=\"{link}\" target=\"_blank\" rel=\"noopener noreferrer\"{link_style}>{display_value}</a></td>"
        return f"<td{cell_style}><a href=\"{link}\" target=\"_blank\" rel=\"noopener noreferrer\"{link_style}>{display_value}</a></td>"

    row_terms = list(enumerate(rows[2:], start=2))
    # headers[0] is blank; headers[2] holds the first query modifier.
    term_headers = [h for h in headers[3:] if h]
    base_header = "term-only"
    html_headers = ["term", base_header]
    html_headers.extend(term_headers)
    column_count = max(len(html_headers), 1)
    uniform_width = 100.0 / column_count
    header_cells = [
        "<th class=\"count-col axis-header axis-term\">{}</th>".format(htmllib.escape(html_headers[0])),
        "<th class=\"count-col\">{}</th>".format(htmllib.escape(html_headers[1] if len(html_headers) > 1 else "")),
    ]
    for header in html_headers[2:]:
        header_cells.append(
            "<th class=\"rotate\"><span class=\"angle\">{}</span></th>".format(htmllib.escape(header))
        )
    header_row = f"<tr>{''.join(header_cells)}</tr>"

    body_rows = []
    for data_index, row in row_terms:
        term = row[0] if len(row) > 0 else ""
        is_compound_only_row = data_index == 2 and not term
        if not term and not is_compound_only_row:
            continue

        row_term = term or "compound-only"
        cells = []
        cells.append(f"<td class=\"axis-term\">{row_term}</td>")

        base_col = row[2] if len(row) > 2 else ""
        if is_compound_only_row:
            cells.append("<td></td>")
        elif base_col and base_col not in PLACEHOLDER_VALUES:
            query = build_term_only_query(term)
            cells.append(
                render_cell(
                    base_col,
                    make_pubmed_url(query),
                    query=query,
                    summary=json_summary_map.get((data_index, 2)),
                    is_heat_cell=True,
                    row_idx=data_index,
                    col=2,
                )
            )
        elif base_col:
            query = build_term_only_query(term)
            cells.append(
                render_cell(
                    base_col,
                    None,
                    query=query,
                    summary=json_summary_map.get((data_index, 2)),
                    is_heat_cell=True,
                    row_idx=data_index,
                    col=2,
                )
            )
        else:
            cells.append("<td></td>")
        for col in range(3, len(headers)):
            header = headers[col]
            if not header:
                continue
            value = row[col] if col < len(row) else ""
            query = build_compound_only_query(header) if is_compound_only_row else build_query(row_term, header)
            cells.append(
                render_cell(
                    value,
                    make_pubmed_url(query),
                    query=query,
                    summary=json_summary_map.get((data_index, col)),
                    is_heat_cell=True,
                    row_idx=data_index,
                    col=col,
                )
            )
        body_rows.append(f"<tr>{''.join(cells)}</tr>")

    html = f"""<!doctype html>
<html lang=\"en\">
  <head>
    <meta charset=\"utf-8\" />
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
    <title>Terpene PubMed Search Results</title>
    <style>
      :root {{
        color-scheme: light;
      }}
      body {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif; margin: 1rem; }}
      .container {{ max-width: 100%; overflow: hidden; }}
      .heat-wrap {{ overflow-x: auto; -webkit-overflow-scrolling: touch; padding-bottom: 0.25rem; }}
      .heat-table {{ border-collapse: collapse; width: 100%; min-width: 900px; table-layout: fixed; font-size: 14px; }}
      th, td {{ border: 1px solid #ddd; padding: 4px 6px; text-align: left; white-space: nowrap; width: {uniform_width:.4f}%; max-width: {uniform_width:.4f}%; }}
      th {{ background: #f5f5f5; position: sticky; top: 0; }}
      th.rotate {{ height: 130px; text-align: left; white-space: nowrap; padding: 0; overflow: visible; }}
      th.rotate .angle {{
        display: inline-block;
        transform: rotate(-45deg);
        transform-origin: left top;
        position: relative;
        left: 0.8rem;
        top: 2.4rem;
      }}
      .count-col {{ text-align: center; }}
      .axis-term {{ position: sticky; left: 0; z-index: 2; background: #fff; white-space: normal; }}
      .axis-term a {{ display: inline; }}
      .heat-legend {{ color: #444; font-size: 12px; margin-bottom: .5rem; }}
      .error {{ color: #b91c1c; }}
      .high-score {{ outline: 2px solid rgba(59, 130, 246, 0.55); }}
      .high-score-list {{ margin: 0 0 1rem; padding-left: 1.2rem; max-height: 160px; overflow: auto; padding-bottom: 0.15rem; }}
      .footer {{ margin-top: 1rem; color: #666; font-size: 12px; }}
      .summary-heading {{ font-weight: 600; margin: 0 0 .35rem; }}
      @media (max-width: 960px) {{
        body {{ margin: 0.75rem; }}
        .heat-table {{ min-width: 760px; font-size: 12px; }}
        th, td {{ padding: 4px 5px; }}
        .heat-legend {{ font-size: 11px; }}
        th.rotate {{ height: 100px; }}
        th.rotate .angle {{ left: .5rem; top: 2rem; }}
      }}
      @media (max-width: 640px) {{
        body {{ margin: 0.5rem; }}
        h1 {{ font-size: 1.35rem; margin: 0.25rem 0; }}
        .heat-legend {{ font-size: 10px; }}
        .high-score-list {{ max-height: 110px; }}
        p, .high-score-list, .footer, .heat-legend {{ font-size: 11px; }}
        .heat-table {{ min-width: 680px; }}
        th.rotate {{ height: auto; white-space: nowrap; }}
        th.rotate .angle {{
          transform: none;
          position: static;
          left: 0;
          top: 0;
          display: inline-block;
          padding: 3px 0;
        }}
      }}
      a {{ color: inherit; text-decoration: none; }}
      a:hover {{ text-decoration: underline; }}
    </style>
  </head>
  <body>
    <main class=\"container\">
      <h1>Terpene PubMed Search Results</h1>
      <p>Generated: {timestamp} UTC</p>
      <p class=\"heat-legend\">Heat map: darker blue means higher PubMed hit counts (term-only + pairwise terpene columns).</p>
      <p class="summary-heading">Top higher-scoring pairwise matches:</p>
      <ul class=\"high-score-list\">
        {summary_html}
      </ul>
      <p><a href=\"results.csv\">Download CSV</a> · <a href=\"results.json\">Download JSON</a></p>
    <div class=\"heat-wrap\">
    <table class=\"heat-table\">
        <thead>
          {header_row}
        </thead>
        <tbody>
          {''.join(body_rows)}
        </tbody>
      </table>
    </div>
      <p class=\"footer\">Source: <a href=\"https://docs.google.com/spreadsheets/d/1VidNfYpvIzB7SA3SePyhHUil0j-XP1TtiakLfWl24pg\" target=\"_blank\" rel=\"noopener noreferrer\">Google Sheet</a></p>
    </main>
  </body>
</html>
"""

    (public_dir / "index.html").write_text(html, encoding="utf-8")
    (public_dir / "README.md").write_text(
        "# Terpene PubMed Search\n\nThis repository is updated nightly by GitHub Actions.\n",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    config = Config(
        source_sheet_id=args.source_sheet_id,
        source_gid=args.source_gid,
        data_dir=_coerce_path(args.output_dir),
        public_dir=_coerce_path(args.public_dir),
        max_rows=args.max_rows,
        request_delay=args.request_delay,
        api_key=args.api_key,
        summarize=args.summarize,
        summary_model=args.summary_model,
        summary_max_abstracts=args.summary_max_abstracts,
        summary_top_cells=args.summary_top_cells,
        summary_api_key=os.getenv("OPENAI_API_KEY"),
        summary_api_base=os.getenv("OPENAI_API_BASE", "https://api.openai.com/v1"),
    )

    rows = normalize_rows(fetch_csv_rows(sheet_csv_url(config.source_sheet_id, config.source_gid)))
    rows = apply_cached_results(rows, read_local_csv(config.data_dir / "results.csv"))
    previous_results = read_local_json(config.public_dir / "results.json")
    previous_summary_cache = load_summary_cache(previous_results)

    updated_rows, json_rows, timestamp = run_search_grid(
        rows,
        config,
        force_refresh=args.force_refresh,
        summary_cache=previous_summary_cache,
    )
    write_csv(config.data_dir / "results.csv", updated_rows)
    write_json(
        config.public_dir / "results.json",
        {
            "generated_at": timestamp,
            "source_sheet_id": config.source_sheet_id,
            "source_gid": config.source_gid,
            "rows": json_rows,
        },
    )
    render_html(config.public_dir, updated_rows[0], updated_rows, json_rows, timestamp)

    # Keep a copy for quick verification in PRs.
    write_csv(config.public_dir / "results.csv", updated_rows)


if __name__ == "__main__":
    main()
