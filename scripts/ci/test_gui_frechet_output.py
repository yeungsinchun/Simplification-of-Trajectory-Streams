#!/usr/bin/env python3
"""Check the GUI's Frechet line against its original six-digit output."""

import json
import os
from pathlib import Path
import pty
import select
import subprocess
import time


ROOT = Path(__file__).resolve().parents[2]
DATASET = "11"


def main() -> None:
    trace = subprocess.run(
        [str(ROOT / "build/simplify"), "--in", DATASET, "--web-server", "--json-stream"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    distance = json.loads(trace.stdout.splitlines()[0])["expected_frechet"]
    assert format(distance, ".17g") != format(distance, ".6g")
    expected = f"Expected Frechet distance: {distance:.6g}"

    master, slave = pty.openpty()
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    process = subprocess.Popen(
        [str(ROOT / "build/simplify_with_gui"), "--in", DATASET],
        cwd=ROOT,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=slave,
        stderr=subprocess.DEVNULL,
    )
    os.close(slave)
    output = b""
    try:
        deadline = time.monotonic() + 30
        marker = b"Expected Frechet distance:"
        while not any(
            line.startswith(marker) and line.endswith(b"\n")
            for line in output.splitlines(keepends=True)
        ) and time.monotonic() < deadline:
            ready, _, _ = select.select([master], [], [], max(0, deadline - time.monotonic()))
            if ready:
                output += os.read(master, 1)
        lines = output.decode(errors="replace").splitlines()
        actual = next((line for line in lines if line.startswith("Expected Frechet distance:")), None)
        assert actual == expected, (actual, expected)
    finally:
        process.terminate()
        process.wait(timeout=5)
        os.close(master)


if __name__ == "__main__":
    main()
