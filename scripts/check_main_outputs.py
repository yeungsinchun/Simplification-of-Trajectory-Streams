#!/usr/bin/env python3
"""Compare all CI simplification outputs byte for byte against a main binary.

Example: python3 scripts/check_main_outputs.py --baseline build/main-simplify
Both binaries write ignored data/<id>/simplify.txt; this check overwrites it.
"""

import argparse
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parent.parent
PAIRS = (
    ("coarse-e", "299", "1"),
    ("fine-e", "0.5", "300"),
    ("mid", "5", "50"),
    ("large-d", "1", "999"),
)


def run(binary: Path, ident: int, epsilon: str, delta: str) -> bytes:
    subprocess.run(
        [str(binary), "--in", str(ident), "--out", "-e", epsilon, "-d", delta],
        cwd=ROOT,
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    return (ROOT / "data" / str(ident) / "simplify.txt").read_bytes()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True,
                        help="Release binary built from main")
    parser.add_argument("--candidate", type=Path, default=Path("build/simplify"))
    args = parser.parse_args()
    baseline = args.baseline.resolve()
    candidate = args.candidate.resolve()
    failures = []
    for label, epsilon, delta in PAIRS:
        for ident in range(1, 11):
            expected = run(baseline, ident, epsilon, delta)
            actual = run(candidate, ident, epsilon, delta)
            if actual != expected:
                failures.append(f"{label} ID {ident}")
        print(f"{label}: 10 cases checked")
    if failures:
        print("Output mismatches: " + ", ".join(failures))
        return 1
    print("All 40 outputs match main byte for byte.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
