#!/usr/bin/env python3
"""Check RFC citations (the rfc-authoring skill).

For each RFC (default: every rfc/NNNN-*.md), checks that every in-text citation anchor
`](#ref-...)` has a matching `<a id="ref-..."></a>` entry under `## References`, that every
entry is cited, and that every entry carries a link. With `--crossref`, also resolves each
DOI through the Crossref API and prints its first author, year and title, to compare
against the citation by eye. Exits 1 on any structural problem or unresolved DOI.

Usage:
    python3 .agents/skills/rfc-authoring/check_citations.py [rfc/0007-*.md ...] [--crossref]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CITED = re.compile(r"\]\(#(ref-[a-z0-9]+)\)")
ENTRY = re.compile(r'<a id="(ref-[a-z0-9]+)"></a>(.*?)(?=\n- <a id=|\Z)', re.S)
DOI_LINK = "](https://doi.org/"


def _dois(entry: str) -> list[str]:
    """DOIs in an entry's Markdown links. Many DOIs contain parentheses, e.g.
    10.1016/S0012-821X(00)00350-2, so read to the link's balancing `)`, not the first."""
    found = []
    start = entry.find(DOI_LINK)
    while start >= 0:
        i = start + len(DOI_LINK)
        depth = 0
        for j in range(i, len(entry)):
            if entry[j] == "(":
                depth += 1
            elif entry[j] == ")":
                if depth == 0:
                    found.append(entry[i:j])
                    break
                depth -= 1
        start = entry.find(DOI_LINK, i)
    return found


def _crossref(doi: str) -> str | None:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(urllib.parse.unquote(doi), safe="/")
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            message = json.load(response)["message"]
    except Exception:  # noqa: BLE001 -- any failure means "didn't resolve"
        return None
    authors = message.get("author") or [{}]
    year = (message.get("issued", {}).get("date-parts") or [[None]])[0][0]
    return f"{authors[0].get('family')} {year}: {(message.get('title') or [''])[0][:90]}"


def check(path: Path, crossref: bool) -> int:
    text = path.read_text()
    body, _, references = text.partition("## References")
    cited = set(CITED.findall(body))
    entries = dict(ENTRY.findall(references))
    problems = 0
    for anchor in sorted(cited - entries.keys()):
        print(f"{path.name}: cited but no entry: {anchor}")
        problems += 1
    for anchor in sorted(entries.keys() - cited):
        print(f"{path.name}: entry never cited: {anchor}")
        problems += 1
    for anchor, entry in sorted(entries.items()):
        dois = _dois(entry)
        if not dois and "https://" not in entry:
            print(f"{path.name}: {anchor} has no link")
            problems += 1
        if crossref:
            for doi in dois:
                found = _crossref(doi)
                print(f"{path.name}: {anchor:24s} {found or 'UNRESOLVED ' + doi}")
                problems += found is None
    return problems


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="*", type=Path)
    parser.add_argument("--crossref", action="store_true", help="resolve each DOI through Crossref")
    args = parser.parse_args()
    files = args.files or sorted((ROOT / "rfc").glob("[0-9][0-9][0-9][0-9]-*.md"))
    problems = sum(check(path, args.crossref) for path in files)
    print(f"{len(files)} RFC(s) checked, {problems} problem(s)")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
