# Terpene PubMed Search Automator

This repo runs the PubMed search cells server-side on a nightly
schedule and publishes results to GitHub Pages.

## What it does

- Pulls the source sheet (`.csv` export) and keeps a working copy in `data/results.csv`.
- Computes search counts for every row by itself in two places:
  - `term-only` column (column B/index 2): `"{term}"[All Fields]`
  - `compound-only` row (`row 2`): each terpene queried as `"{terpene}"[All Fields]`
- For all pairwise term+terpene cells, it computes log-scaled heatmap intensities directly in the HTML table.
- By default, existing numeric cells are kept and only placeholder cells are refreshed.
  Use `--force-refresh` to recompute every searchable cell in the matrix in one run.
- You can optionally summarize results by adding `--summarize`. This:
  - pulls top PMIDs per query,
  - fetches their abstracts from PubMed,
  - and calls an OpenAI-compatible LLM to generate a short evidence summary.
  The summary appears in cell hover tooltips and is included in `public/results.json`.
- Publishes:
  - `public/index.html` (human-readable results table with links to PubMed searches),
  - `public/results.csv` (full exported matrix),
  - `public/results.json` (machine-readable matrix payload).

## Source sheet

Configured by defaults in `scripts/run_pubmed_nightly.py`:

- Sheet ID: `1VidNfYpvIzB7SA3SePyhHUil0j-XP1TtiakLfWl24pg`
- GID: `363407775`

## Running locally

```bash
python scripts/run_pubmed_nightly.py
```

Optional flags:

- `--source-sheet-id`, `--source-gid`
- `--output-dir` (default: `data`)
- `--public-dir` (default: `public`)
- `--request-delay` (seconds between PubMed calls)
- `--api-key` / `NCBI_API_KEY` (optional, for higher rate limits)
- `--summarize` (download abstracts and generate summaries with OpenAI-style API)
- `--summary-model` / `OPENAI_MODEL` (default `gpt-4o-mini`)
- `--summary-max-abstracts` / `PUBMED_SUMMARY_MAX_ABSTRACTS` (default `8`)
- `--summary-top-cells` / `PUBMED_SUMMARY_TOP_CELLS` (default `40`)
- `OPENAI_API_KEY` and optional `OPENAI_API_BASE` when using `--summarize`

## GitHub Actions

- `/.github/workflows/nightly-pubmed.yml` runs every night at `02:30 UTC`.
- It commits refreshed outputs and deploys the `public/` folder to GitHub Pages.
- If `OPENAI_API_KEY` is configured in repository secrets, the workflow also runs with `--summarize` and includes tooltip summaries in the generated JSON/HTML.
- You can adjust timing, output folder, or query formula behavior in that workflow file.
- After pushing this repository, enable GitHub Pages in Settings → Pages with "Source: GitHub Actions".
- The workflow writes one table cell per query with a direct PubMed search link.

## Query formula

Current search query template is:

```text
"{query_term}"[All Fields] AND "{terpene_term}"[All Fields]
```

If you want a different PubMed syntax (for example `[Title/Abstract]`), edit `build_query()` in
`scripts/run_pubmed_nightly.py`.
