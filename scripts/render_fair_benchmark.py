#!/usr/bin/env python3
"""Render the self-contained Lavish board from saved benchmark evidence."""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results/fair-core")
    parser.add_argument("--output", type=Path, default=ROOT / ".lavish/sots-fair-bench/index.html")
    args = parser.parse_args()
    evidence = {phase: json.loads((args.results / f"{phase}.json").read_text())
                for phase in ("before", "after")}
    evidence["notes"] = json.loads((args.results / "notes.json").read_text())
    if evidence["before"]["config"] != evidence["after"]["config"]:
        raise ValueError("Incompatible phases")
    html = (ROOT / "scripts/fair_benchmark_board.html").read_text()
    # JSON is data, never HTML: prevent a path/string from closing the script.
    payload = json.dumps(evidence, allow_nan=False).replace("<", "\\u003c")
    html = html.replace("__EVIDENCE__", payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html)
    print(args.output)


if __name__ == "__main__":
    main()
