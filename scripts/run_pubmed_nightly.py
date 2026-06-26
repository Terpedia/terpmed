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
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional
from urllib.parse import urlencode

import datetime
import html as htmllib

import urllib.request


DEFAULT_SOURCE_SHEET_ID = "1VidNfYpvIzB7SA3SePyhHUil0j-XP1TtiakLfWl24pg"
DEFAULT_SOURCE_GID = "363407775"
PUBMED_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
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


def pubmed_search_count(query: str, api_key: Optional[str]) -> str:
    params = {
        "db": "pubmed",
        "retmode": "json",
        "retmax": 0,
        "term": query,
    }
    if api_key:
        params["api_key"] = api_key

    # Keep implementation dependency-free for GitHub Actions image compatibility.
    url = f"{PUBMED_BASE_URL}?{urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    try:
        payload = json.loads(raw)
        return payload["esearchresult"]["count"]
    except Exception as err:  # pragma: no cover - defensive for malformed responses
        raise RuntimeError(f"Unexpected response for query {query!r}: {raw[:200]}") from err


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


def run_search_grid(rows: List[List[str]], config: Config, force_refresh: bool = False) -> tuple[List[List[str]], List[List[dict]], str]:
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
            base_query = build_term_only_query(term)
            base_count = fetch_with_retry(base_query, config.api_key, config.request_delay)
            row_out[2] = base_count
            base_source = "queried"
            json_cells.append(
                {
                    "column": headers[2],
                    "value": base_count,
                    "query": base_query,
                    "link": make_pubmed_url(base_query),
                    "source": base_source,
                }
            )
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
                query = build_compound_only_query(compound)
                if is_pair_to_compute:
                    count = fetch_with_retry(query, config.api_key, config.request_delay)
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
                query = build_query(term, compound)
                count = fetch_with_retry(query, config.api_key, config.request_delay)
                source = "queried"
                row_out[col_idx] = count
                time.sleep(config.request_delay)

            if is_compound_axis_row:
                source = "compound-only"
            json_cells.append(
                {
                    "column": headers[col_idx],
                    "value": count,
                    "query": query,
                    "link": make_pubmed_url(query),
                    "source": source,
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

    def render_cell(
        value: str,
        link: Optional[str],
        is_heat_cell: bool = False,
        row_idx: Optional[int] = None,
        col: Optional[int] = None,
    ) -> str:
        display_value = compact_count(value)
        if value in NON_RESULT_VALUES or not value:
            return "<td></td>"
        if value.startswith("ERROR:"):
            return f"<td class=\"error\">{value}</td>"
        cell_style = ""
        link_style = ""
        cell_class = ""
        summary = None
        if is_heat_cell and row_idx is not None and col is not None:
            summary = summary_map.get((row_idx, col))
        if summary:
            cell_class = " class=\"high-score\""
        if is_heat_cell and has_heat_values:
            count = parse_count(value)
            if count is not None:
                style, text_color = heatstyle(count, heat_min, heat_max)
                if style:
                    cell_style = f" style=\"{style}\""
                    link_style = f" style=\"color:{text_color}\""
        if not link:
            if summary:
                return f"<td{cell_class} title=\"{htmllib.escape(summary)}\"{cell_style}>{display_value}</td>"
            return f"<td{cell_style}{cell_class}>{display_value}</td>"
        if summary:
            return f"<td{cell_class} title=\"{htmllib.escape(summary)}\"{cell_style}><a href=\"{link}\" target=\"_blank\" rel=\"noopener noreferrer\"{link_style}>{display_value}</a></td>"
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
        "<th class=\"count-col\">{}</th>".format(htmllib.escape(html_headers[0])),
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
        cells.append(f"<td>{row_term}</td>")

        base_col = row[2] if len(row) > 2 else ""
        if is_compound_only_row:
            cells.append("<td></td>")
        elif base_col and base_col not in PLACEHOLDER_VALUES:
            query = build_term_only_query(term)
            cells.append(render_cell(base_col, make_pubmed_url(query), is_heat_cell=True, row_idx=data_index, col=2))
        elif base_col:
            cells.append(render_cell(base_col, None, is_heat_cell=True, row_idx=data_index, col=2))
        else:
            cells.append("<td></td>")
        for col in range(3, len(headers)):
            header = headers[col]
            if not header:
                continue
            value = row[col] if col < len(row) else ""
            query = build_compound_only_query(header) if is_compound_only_row else build_query(row_term, header)
            cells.append(render_cell(value, make_pubmed_url(query), is_heat_cell=True, row_idx=data_index, col=col))
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
      .heat-legend {{ color: #444; font-size: 12px; margin-bottom: .5rem; }}
      .error {{ color: #b91c1c; }}
      .high-score {{ outline: 2px solid rgba(59, 130, 246, 0.55); }}
      .high-score-list {{ margin: 0 0 1rem; padding-left: 1.2rem; }}
      .footer {{ margin-top: 1rem; color: #666; font-size: 12px; }}
      @media (max-width: 960px) {{
        body {{ margin: 0.75rem; }}
        .heat-table {{ min-width: 760px; font-size: 12px; }}
        th, td {{ padding: 4px 5px; }}
        .heat-legend {{ font-size: 11px; }}
      }}
      @media (max-width: 640px) {{
        body {{ margin: 0.5rem; }}
        h1 {{ font-size: 1.35rem; margin: 0.25rem 0; }}
        .heat-table {{ min-width: 680px; }}
        th.rotate {{ height: auto; white-space: normal; }}
        th.rotate .angle {{
          transform: none;
          position: static;
          display: block;
          left: 0;
          top: 0;
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
      <p>Top higher-scoring pairwise matches:</p>
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
    )

    rows = normalize_rows(fetch_csv_rows(sheet_csv_url(config.source_sheet_id, config.source_gid)))

    updated_rows, json_rows, timestamp = run_search_grid(rows, config, force_refresh=args.force_refresh)
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
