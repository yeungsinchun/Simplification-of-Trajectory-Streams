#!/usr/bin/env python3
"""Small unit test for bench_report.py — runs locally without CI artifacts."""

import json
import tempfile
import pathlib
import subprocess
import sys

def _make_fake_doc(label="mid-e-d300-small", epsilon=5, delta=50, gated=10, mean_orig=4.5, mean_new=2.0, mean_ratio=0.444, overall_ok=True):
    # 10 cases, each with orig/new ms and phase ops (all gated, no floor)
    cases = []
    for i in range(11, 21):
        orig_ms = 4.0 + (i % 3) * 0.2
        new_ms = orig_ms * 0.44
        # simulate TIMER_MS: intersect, find_F, hull_Gi, boundary_P
        orig_ops = {"intersect": 2.0, "find_F": 1.2, "hull_Gi": 0.3, "boundary_P": 0.05, "total": 3.55}
        new_ops = {"intersect": 0.8, "find_F": 0.5, "hull_Gi": 0.3, "boundary_P": 0.05, "total": 1.65}
        cases.append({
            "id": i,
            "e": epsilon,
            "d": delta,
            "orig_ms": orig_ms,
            "new_ms": new_ms,
            "ratio": new_ms/orig_ms,
            "gated": True,
            "status": "OK",
            "orig_ops": orig_ops,
            "new_ops": new_ops,
            "ops": new_ops,
        })
    return {
        "label": label,
        "epsilon": epsilon,
        "delta": delta,
        "slowdown_limit": 1.2,
        "min_bench_ms": 0.0,
        "bench_runs": 5,
        "mean_limit": 1.05,
        "gated_count": gated,
        "mean_orig_ms": mean_orig,
        "mean_new_ms": mean_new,
        "mean_ratio": mean_ratio,
        "mean_ok": True,
        "overall_ok": overall_ok,
        "cases": cases,
    }

def test_basic():
    with tempfile.TemporaryDirectory() as td:
        td_path = pathlib.Path(td)
        # create two fake reports in nested layout like CI artifacts
        for label in ["mid-e-d300-small", "fine-e-d300-large"]:
            sub = td_path / f"benchmark-{label}"
            sub.mkdir(parents=True)
            doc = _make_fake_doc(label=label,
                                 epsilon=5 if "mid" in label else 0.5,
                                 delta=50 if "mid" in label else 200,
                                 gated=10,
                                 mean_orig=4.2 if "mid" in label else 34.0,
                                 mean_new=1.9 if "mid" in label else 15.8,
                                 mean_ratio=0.45 if "mid" in label else 0.46)
            (sub / "benchmark.json").write_text(json.dumps(doc, indent=2))
        out_md = td_path / "summary.md"
        out_html = td_path / "report.html"
        out_comment = td_path / "comment.md"
        # invoke bench_report.py
        repo_root = pathlib.Path(__file__).resolve().parents[2]
        script = repo_root / "scripts" / "ci" / "bench_report.py"
        if not script.exists():
            script = pathlib.Path("scripts/ci/bench_report.py")
        proc = subprocess.run(
            [sys.executable, str(script),
             "--reports-dir", str(td_path),
             "--output-md", str(out_md),
             "--output-html", str(out_html),
             "--output-comment", str(out_comment),
             "--run-url", "https://github.com/example/repo/actions/runs/12345",
             "--sha", "abc123def456"],
            capture_output=True, text=True
        )
        assert proc.returncode == 0, f"bench_report failed: {proc.stdout}\n{proc.stderr}"
        md = out_md.read_text()
        html = out_html.read_text()
        comment = out_comment.read_text()

        # Check that markdown contains overall table headers and known labels
        assert "Overall per configuration" in md
        assert "mid-e-d300-small" in md
        assert "fine-e-d300-large" in md
        assert "orig ms" in md.lower()
        assert "speedup" in md.lower()
        # Check per-phase sections
        assert "Per-phase" in md
        assert "intersect" in md
        assert "find_F" in md
        # Check HTML has tables and inline bar
        assert "<table" in html
        assert "bar-track" in html or "bar " in html
        assert "mid-e-d300-small" in html
        assert "https://github.com/example/repo/actions/runs/12345" in html
        # Check comment has marker and headline
        assert "<!-- bench-report -->" in comment
        assert "Benchmark Report" in comment
        assert "artifact" in comment
        print("test_basic passed")

        # Test legacy handling: one doc with only 'ops' (no orig_ops)
        legacy_path = td_path / "benchmark-legacy" / "benchmark.json"
        legacy_path.parent.mkdir(exist_ok=True)
        legacy_doc = _make_fake_doc(label="legacy-small")
        # strip orig_ops to simulate old artifact
        for c in legacy_doc["cases"]:
            c.pop("orig_ops", None)
        legacy_path.write_text(json.dumps(legacy_doc))
        proc2 = subprocess.run(
            [sys.executable, str(script), "--reports-dir", str(td_path), "--output-md", str(td_path/"summary2.md")],
            capture_output=True, text=True
        )
        assert proc2.returncode == 0, f"legacy handling failed: {proc2.stderr}"
        print("legacy handling passed")

        for c in legacy_doc["cases"][:3]:
            c["gated"] = False
            c["status"] = "SKIP_FAST"
        legacy_doc["gated_count"] = 7
        legacy_path.write_text(json.dumps(legacy_doc))
        partial_dir = td_path / "partial"
        partial_dir.mkdir()
        (partial_dir / "benchmark.json").write_text(json.dumps(legacy_doc))
        partial_md = td_path / "partial.md"
        partial_html = td_path / "partial.html"
        partial_comment = td_path / "partial-comment.md"
        proc3 = subprocess.run(
            [sys.executable, str(script), "--reports-dir", str(partial_dir),
             "--output-md", str(partial_md), "--output-html", str(partial_html),
             "--output-comment", str(partial_comment), "--expected-configs", "2"],
            capture_output=True, text=True
        )
        assert proc3.returncode == 0, proc3.stderr
        md = partial_md.read_text()
        html = partial_html.read_text()
        comment = partial_comment.read_text()
        assert "1/2 configurations reported" in md and "1/2 configurations reported" in comment
        assert "INCOMPLETE MATRIX (1/2 reported)" in html
        assert "1/2 configurations passed" in md and "1/2 passed" in html and "1/2 passed" in comment
        assert "7/10 IDs gated" in md and "7/10 IDs gated" in html and "7/10 IDs gated" in comment
        assert "<td>no</td>" in html and "<td>yes</td>" in html
        assert "| intersect | n/a | n/a |" in md
        assert "<td class=\"mono\">intersect</td><td>n/a</td><td>n/a</td>" in html
        assert "intersect: n/a →" in comment

        phase_dir = td_path / "phase-matrix"
        phase_dir.mkdir()
        missing_phase = _make_fake_doc(label="mid-e-d300-small")
        complete_phase = _make_fake_doc(label="mid-e-d300-large")
        for doc in (missing_phase, complete_phase):
            for case in doc["cases"]:
                case["new_ops"]["prepare_clip_polygon"] = 0.2
                case["orig_ops"]["prepare_clip_polygon"] = 0.4
                case["new_ops"]["total"] += 0.2
                case["orig_ops"]["total"] += 0.4
                case["new_ops"]["total"] += 0.1
        missing_phase["cases"][0]["orig_ops"].pop("prepare_clip_polygon")
        for doc in (missing_phase, complete_phase):
            sub = phase_dir / doc["label"]
            sub.mkdir()
            (sub / "benchmark.json").write_text(json.dumps(doc))
        phase_md = td_path / "phase.md"
        phase_html = td_path / "phase.html"
        phase_comment = td_path / "phase-comment.md"
        proc4 = subprocess.run(
            [sys.executable, str(script), "--reports-dir", str(phase_dir),
             "--output-md", str(phase_md), "--output-html", str(phase_html),
             "--output-comment", str(phase_comment), "--expected-configs", "2"],
            capture_output=True, text=True
        )
        assert proc4.returncode == 0, proc4.stderr
        md = phase_md.read_text()
        html = phase_html.read_text()
        comment = phase_comment.read_text()
        assert "| prepare_clip_polygon | n/a | n/a | 4.00 |" in md
        assert "| mid-e-d300-small |" in md and "| mid-e-d300-large |" in md
        assert "<td class=\"mono\">prepare_clip_polygon</td><td>n/a</td><td>n/a</td><td>4.00</td>" in html
        small_detail = html.split('<h3 class="mono">mid-e-d300-small ', 1)[1].split('<h3 class="mono">', 1)[0]
        large_detail = html.split('<h3 class="mono">mid-e-d300-large ', 1)[1].split('<h3 class="mono">', 1)[0]
        assert '<td class="mono">prepare_clip_polygon</td><td>n/a</td><td>n/a</td>' in small_detail
        assert '<td class="mono">prepare_clip_polygon</td><td>4.00</td>' in large_detail
        assert "prepare_clip_polygon: n/a → 4.00 ms (n/a, n/a →" in comment
        assert "| _other_ (total − leaves) | n/a | n/a | 2.00 |" in md
        assert '<td class="mono"><em>other</em></td><td>n/a</td><td>n/a</td><td>2.00</td>' in html

if __name__ == "__main__":
    test_basic()
    print("all tests passed")
