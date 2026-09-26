#!/usr/bin/env python3
"""
Benchmark report aggregation and rendering for CI.

Reads per-configuration benchmark.json artifacts (produced by benchmark.yml)
and produces:
  - Markdown job summary (overall per-config table + per-phase before/after)
  - Self-contained HTML report (tables + inline CSS bar charts)
  - Sticky PR comment snippet (headline + link)

The script is intended to run both in CI (after downloading artifacts) and
locally on saved logs:

  python scripts/ci/bench_report.py \
    --reports-dir bench-report \
    --output-md summary.md \
    --output-html report.html \
    --output-comment comment.md \
    --run-url https://github.com/org/repo/actions/runs/123

If --reports-dir contains downloaded artifacts in subdirectories
(benchmark-<label>/benchmark.json) it recurses. It also handles legacy
JSON where only `ops` (new) is present; newer artifacts contain `orig_ops`
and `new_ops`.

No external dependencies; only stdlib.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import glob
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Any, Optional

# Mapping for friendly phase names; raw timer names are kept if not matched.
# Parent timers that should not be counted as leaf phases.
PARENT_TIMERS = {"total", "get_longest_stab", "get_longest_stab_web", "simplify"}

def collect_reports(reports_dir: Path) -> List[Path]:
    """Find all benchmark.json files under reports_dir recursively."""
    # Direct file
    candidates: List[Path] = []
    if (reports_dir / "benchmark.json").is_file():
        candidates.append(reports_dir / "benchmark.json")
    # Recursive glob
    for p in reports_dir.rglob("benchmark.json"):
        if p not in candidates:
            candidates.append(p)
    # Also allow benchmark-*.json at top level
    for p in reports_dir.glob("benchmark-*.json"):
        if p not in candidates:
            candidates.append(p)
    return sorted(candidates)

def load_report(path: Path) -> Optional[Dict[str, Any]]:
    try:
        with path.open() as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"warn: failed to load {path}: {e}", file=sys.stderr)
        return None

def fmt_ms(v: Optional[float]) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "n/a"
    if math.isinf(v):
        return "inf"
    return f"{v:.3f}"

def fmt_ratio(v: Optional[float]) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "n/a"
    if math.isinf(v):
        return "inf"
    return f"{v:.2f}×"

def speedup_str(ratio: Optional[float]) -> str:
    if ratio is None or math.isnan(ratio) or math.isinf(ratio):
        return "n/a"
    return f"{ratio:.2f}×"

def classify_timer(name: str) -> str:
    n = name.lower()
    if n == "total" or n == "get_longest_stab" or n == "get_longest_stab_web":
        return "__parent__"
    if "intersect" in n or n == "clip":
        return "intersect"
    if "find_f" in n:
        return "find_F"
    if "hull" in n:
        return "Gi hull"
    if "boundary" in n:
        return "boundary_P"
    if "prep" in n:
        return "Gi prep"
    return name

def aggregate_configs(docs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Normalize and sort docs for reporting."""
    # Define tier order for deterministic sorting matching workflow
    tier_order = {
        "extra-coarse": 0,
        "coarse": 1,
        "mid": 2,
        "fine": 3,
        "extra-fine": 4,
    }
    def sort_key(d: Dict[str, Any]):
        label = d.get("label", "")
        # Extract tier prefix: split by "-" first token
        # label like "extra-coarse-e-d300-small"
        tier = label.split("-e-")[0] if "-e-" in label else label.split("-")[0]
        # Normalize extra-coarse vs extra-fine etc
        # Determine tier group
        for k in tier_order:
            if label.startswith(k):
                tier_rank = tier_order[k]
                break
        else:
            tier_rank = 99
        # delta numer: 300 vs 1000 -> d300 before d1000
        delta_group = 0 if "d300" in label else 1
        size_group = 0 if "small" in label else 1
        return (tier_rank, delta_group, size_group, label)
    docs_sorted = sorted(docs, key=sort_key)
    return docs_sorted

def compute_phase_totals(docs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Compute per-label and global phase aggregates from TIMER_MS ops.

    Returns dict with:
      per_label: list of {label, orig_totals, new_totals, orig_total_time, new_total_time, cases}
      global_orig, global_new dicts
      all_phases: sorted list of phase names encountered (leaf only)
    """
    per_label = []
    global_orig: Dict[str, float] = defaultdict(float)
    global_new: Dict[str, float] = defaultdict(float)
    all_phases_set = set()

    for doc in docs:
        label = doc.get("label", "unknown")
        cases = doc.get("cases", [])
        orig_totals: Dict[str, float] = defaultdict(float)
        new_totals: Dict[str, float] = defaultdict(float)
        # collect totals per case
        for case in cases:
            # Support both new schema (orig_ops/new_ops) and legacy (ops)
            orig_ops = case.get("orig_ops")
            new_ops = case.get("new_ops")
            # fallback legacy: case["ops"] is new_ops
            if new_ops is None:
                new_ops = case.get("ops", {})
            if orig_ops is None:
                orig_ops = {}
                # legacy has no orig_ops
                # try orig_ops field maybe empty
                # keep empty
                pass
            # orig_ops/new_ops could be dict or empty
            if isinstance(orig_ops, dict):
                for k, v in orig_ops.items():
                    try:
                        fv = float(v)
                    except: 
                        continue
                    orig_totals[k] += fv
                    all_phases_set.add(k)
            if isinstance(new_ops, dict):
                for k, v in new_ops.items():
                    try:
                        fv = float(v)
                    except:
                        continue
                    new_totals[k] += fv
                    all_phases_set.add(k)
        # Compute denominator for shares
        # Total time: use "total" timer if present, else sum of leaf timers
        # Leaves are all except parents
        leaf_orig_sum = sum(v for k, v in orig_totals.items() if k not in PARENT_TIMERS)
        leaf_new_sum = sum(v for k, v in new_totals.items() if k not in PARENT_TIMERS)
        orig_total_time = orig_totals.get("total", leaf_orig_sum) if orig_totals else leaf_orig_sum
        new_total_time = new_totals.get("total", leaf_new_sum) if new_totals else leaf_new_sum
        # If total exists but is 0, fallback to leaf
        if orig_total_time == 0 and leaf_orig_sum > 0:
            orig_total_time = leaf_orig_sum
        if new_total_time == 0 and leaf_new_sum > 0:
            new_total_time = leaf_new_sum
        # Accumulate global
        for k, v in orig_totals.items():
            global_orig[k] += v
        for k, v in new_totals.items():
            global_new[k] += v

        per_label.append({
            "label": label,
            "doc": doc,
            "orig_totals": dict(orig_totals),
            "new_totals": dict(new_totals),
            "orig_total_time": orig_total_time,
            "new_total_time": new_total_time,
            "leaf_orig_sum": leaf_orig_sum,
            "leaf_new_sum": leaf_new_sum,
        })

    # Sort global phases: put known phases first in logical order, then alphabetic
    canonical_order = ["intersect", "find_F", "hull_Gi", "boundary_P", "Gi prep", "Gi hull"]
    # normalize names for ordering: map raw to friendly? Keep raw order but prioritize intersect/find_F/hull
    priority = {
        "intersect": 0,
        "clip": 0,
        "find_F": 1,
        "hull_Gi": 2,
        "Gi hull": 2,
        "boundary_P": 3,
        "Gi prep": 4,
        "prepare_clip_polygon": 4,
    }
    def phase_sort_key(name: str) -> Tuple[int, str]:
        if name in priority:
            return (priority[name], name)
        # classify
        cls = classify_timer(name)
        if cls in priority:
            return (priority[cls], name)
        if name in PARENT_TIMERS:
            return (99, name)
        return (10, name)
    all_phases = sorted(all_phases_set, key=phase_sort_key)
    # Filter out parent timers from display list for leaf phases, but keep them for reference
    leaf_phases = [p for p in all_phases if p not in PARENT_TIMERS]
    return {
        "per_label": per_label,
        "global_orig": dict(global_orig),
        "global_new": dict(global_new),
        "all_phases": all_phases,
        "leaf_phases": leaf_phases,
    }

def render_markdown(docs: List[Dict[str, Any]], phase_info: Dict[str, Any], run_url: str, run_id: str) -> str:
    lines: List[str] = []
    lines.append("# Benchmark Report")
    lines.append("")
    if run_url:
        lines.append(f"Run: [{run_url}]({run_url})")
    elif run_id:
        lines.append(f"Run ID: {run_id}")
    lines.append("")
    # Headline numbers
    gated_docs = [d for d in docs if d.get("gated_count", 0) > 0]
    all_ratios = [d.get("mean_ratio") for d in gated_docs if isinstance(d.get("mean_ratio"), (int,float)) and not math.isnan(d.get("mean_ratio", float('nan')))]
    if all_ratios:
        min_ratio = min(all_ratios)
        max_ratio = max(all_ratios)
        # speedup is orig/new? Wait mean_ratio is new/orig (from benchmark.yml: mean_ratio = mean_new/mean_orig)
        # So speedup = orig/new = 1/mean_ratio
        # But earlier report says speedup 2.04-2.40x where new is faster => ratio <1, speedup >1.
        # We should display speedup = 1/ratio.
        speedups = [1/r if r and r !=0 else float('inf') for r in all_ratios]
        min_sp = min(speedups)
        max_sp = max(speedups)
        # mean speedup as 1 / mean(mean_ratio) ??? simpler avg speedup = avg(orig/new)
        mean_sp = sum(speedups)/len(speedups) if speedups else float('nan')
        lines.append(f"**Headline:** {len(gated_docs)}/{len(docs)} configurations gated (orig ≥1 ms).")
        lines.append(f"Speedups (orig/new): **{min_sp:.2f}–{max_sp:.2f}×** (mean {mean_sp:.2f}×) on gated mean (new/orig ratios {min(all_ratios):.3f}–{max(all_ratios):.3f}).")
        lines.append("")
        # Also show pass/fail counts
        pass_cnt = sum(1 for d in docs if d.get("overall_ok"))
        lines.append(f"Gating: {pass_cnt}/{len(docs)} configurations passed (mean ≤1.05× and per-ID ≤1.20×).")
        lines.append("")
    else:
        # No gated configs: show Sig totals?
        lines.append(f"{len(docs)} configurations (no gated IDs above floor).")
        lines.append("")

    # Overall per-configuration table
    lines.append("## Overall per configuration")
    lines.append("")
    lines.append("| Configuration | ε | δ | size | orig ms (gated mean) | new ms (gated mean) | speedup (orig/new) | gated IDs | status |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for doc in docs:
        label = doc.get("label", "n/a")
        eps = doc.get("epsilon", "")
        delta = doc.get("delta", "")
        # size derived from label
        size = "small" if "small" in label else "large" if "large" in label else doc.get("size", "")
        mean_orig = doc.get("mean_orig_ms")
        mean_new = doc.get("mean_new_ms")
        mean_ratio = doc.get("mean_ratio")
        gated = doc.get("gated_count", 0)
        overall_ok = doc.get("overall_ok", False)
        # speedup = orig/new
        if isinstance(mean_ratio, (int,float)) and mean_ratio and mean_ratio != 0:
            speedup = 1/mean_ratio
            speed_s = f"{speedup:.2f}×"
        else:
            speed_s = "n/a"
        # Status: if gated==0 => SKIP (shown as below floor)
        # else if overall_ok => PASS else FAIL
        if gated == 0:
            status = "SKIP_FLOOR"
        else:
            status = "PASS" if overall_ok else "FAIL"
            if not doc.get("mean_ok", True):
                status += " (mean)"
            # Check per-case fails
            fails = [c for c in doc.get("cases", []) if c.get("status") == "FAIL_SLOW"]
            if fails:
                status += f" ({len(fails)} slow)"
        # Format ms
        orig_s = fmt_ms(mean_orig) if gated else fmt_ms(mean_orig) if mean_orig else "Σ"
        new_s = fmt_ms(mean_new) if gated else fmt_ms(mean_new) if mean_new else "Σ"
        # For n/a with no gated, show sum? Already mean_orig is n/a when gated 0
        # In that case doc has mean_orig_ms None, but we could show total? Fallback to showing n/a
        lines.append(f"| {label} | {eps} | {delta:.6g} | {size} | {orig_s} | {new_s} | {speed_s} | {gated}/10 | {status} |")
    lines.append("")
    lines.append(f"_Thresholds: per-ID new ≤ orig×1.20, mean(new) ≤ mean(orig)×1.05 (gated IDs only)._")
    lines.append("")

    # Per-phase global table
    lines.append("## Per-phase before/after (global sum over IDs)")
    lines.append("")
    lines.append("_Aggregated from `TIMER_MS` lines. Each ID's phase time is one `--time` run; values summed over the 10 IDs per configuration then over all configurations._")
    lines.append("")
    global_orig = phase_info["global_orig"]
    global_new = phase_info["global_new"]
    leaf_phases = phase_info["leaf_phases"]
    all_phases = phase_info["all_phases"]
    # Compute global total times for shares (using total timer if present else sum leaves)
    leaf_orig_global = sum(global_orig.get(k, 0) for k in leaf_phases)
    leaf_new_global = sum(global_new.get(k, 0) for k in leaf_phases)
    global_orig_total = global_orig.get("total", leaf_orig_global) or leaf_orig_global or 1
    global_new_total = global_new.get("total", leaf_new_global) or leaf_new_global or 1
    # Build phase table: phase | orig ms | orig share | new ms | new share | speedup
    lines.append("| Phase | orig ms | orig share | new ms | new share | speedup (orig/new) |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    # Show leaf phases sorted
    for phase in leaf_phases:
        o_ms = global_orig.get(phase, 0)
        n_ms = global_new.get(phase, 0)
        o_share = (o_ms / global_orig_total * 100) if global_orig_total else 0
        n_share = (n_ms / global_new_total * 100) if global_new_total else 0
        # speedup orig/new
        sp = (o_ms / n_ms) if n_ms > 0 else float('inf') if o_ms>0 else 0
        sp_s = f"{sp:.2f}×" if sp != float('inf') and sp !=0 else "inf" if sp==float('inf') else "n/a"
        lines.append(f"| {phase} | {o_ms:.2f} | {o_share:.1f}% | {n_ms:.2f} | {n_share:.1f}% | {sp_s} |")
    # Other
    if global_orig_total > leaf_orig_global or global_new_total > leaf_new_global:
        o_other = max(0, global_orig_total - leaf_orig_global)
        n_other = max(0, global_new_total - leaf_new_global)
        o_share = (o_other / global_orig_total * 100) if global_orig_total else 0
        n_share = (n_other / global_new_total * 100) if global_new_total else 0
        sp = (o_other / n_other) if n_other>0 else float('inf') if o_other>0 else 0
        sp_s = f"{sp:.2f}×" if sp != float('inf') and sp!=0 else "inf" if sp==float('inf') else "n/a"
        lines.append(f"| _other_ (total − leaves) | {o_other:.2f} | {o_share:.1f}% | {n_other:.2f} | {n_share:.1f}% | {sp_s} |")
    # Totals
    # global_orig_total and new may be total timer, which includes overhead; show
    lines.append(f"| **total** | {global_orig_total:.2f} | 100% | {global_new_total:.2f} | 100% | {global_orig_total/(global_new_total) if global_new_total else 0:.2f}× |")
    lines.append("")
    # Also note overhead
    # Add per-configuration phase details as collapsible? Markdown simple table per config with limited phases
    # Provide compact per-config phase speedup heat? For markdown, add a summarized per-label phase share table
    # We'll add a short note and delegate details to HTML
    if phase_info["per_label"]:
        has_phase_data = any(pl["orig_totals"] or pl["new_totals"] for pl in phase_info["per_label"])
        if has_phase_data:
            lines.append("### Per-configuration phase speedup (orig/new)")
            lines.append("")
            # Header: config | intersect | find_F | hull | prep | boundary | ...
            # Determine top phases to show: intersect, find_F, hull_Gi, boundary_P
            # Use leaf_phases top 4-5
            top_phases = [p for p in leaf_phases if p not in PARENT_TIMERS][:6]
            if not top_phases:
                top_phases = leaf_phases[:6]
            header = "| Configuration | " + " | ".join(top_phases) + " |"
            lines.append(header)
            lines.append("|---|" + "---:|" * len(top_phases))
            for pl in phase_info["per_label"]:
                label = pl["label"]
                row = f"| {label} |"
                for ph in top_phases:
                    o = pl["orig_totals"].get(ph, 0)
                    n = pl["new_totals"].get(ph, 0)
                    if o == 0 and n == 0:
                        row += " n/a |"
                    else:
                        sp = (o / n) if n > 0 else float('inf')
                        row += f" {sp:.2f}× |"
                lines.append(row)
            lines.append("")
            lines.append("_Full per-phase before/after ms and shares are in the HTML report artifact._")
            lines.append("")

    lines.append("---")
    lines.append("_Generated by `scripts/ci/bench_report.py` from `benchmark.json` artifacts. Gated means use `orig_ms ≥1.0 ms`._")
    if run_url:
        lines.append(f"_Run: {run_url}_")
    lines.append("")
    return "\n".join(lines)

def render_html(docs: List[Dict[str, Any]], phase_info: Dict[str, Any], run_url: str, run_id: str, sha: str = "") -> str:
    # Compute globals for headline
    gated_docs = [d for d in docs if d.get("gated_count", 0) > 0]
    all_ratios = [d.get("mean_ratio") for d in gated_docs if isinstance(d.get("mean_ratio"), (int,float))]
    # speedups
    speedups = []
    for r in all_ratios:
        if r and r != 0:
            speedups.append(1/r)
    headline_speedup_range = f"{min(speedups):.2f}–{max(speedups):.2f}×" if speedups else "n/a"
    mean_sp = sum(speedups)/len(speedups) if speedups else 0
    pass_cnt = sum(1 for d in docs if d.get("overall_ok"))
    # Build HTML
    # Find max orig mean for bar scaling
    max_orig_mean = max([d.get("mean_orig_ms") or 0 for d in docs if d.get("mean_orig_ms")] + [1])
    max_new_mean = max([d.get("mean_new_ms") or 0 for d in docs if d.get("mean_new_ms")] + [1])
    max_mean = max(max_orig_mean, max_new_mean)
    # Phase totals for bar scaling
    global_orig = phase_info["global_orig"]
    global_new = phase_info["global_new"]
    leaf_phases = phase_info["leaf_phases"]
    # Colors for phases
    phase_colors = {
        "intersect": "#ff6b78",
        "clip": "#ff6b78",
        "find_F": "#4f9dff",
        "hull_Gi": "#3ddc97",
        "Gi hull": "#3ddc97",
        "boundary_P": "#f6c453",
        "Gi prep": "#b692ff",
        "prepare": "#b692ff",
    }
    # Helper to get color
    def phase_color(name: str) -> str:
        if name in phase_colors:
            return phase_colors[name]
        c = classify_timer(name)
        if c in phase_colors:
            return phase_colors[c]
        # hash to color
        hues = ["#22d3ee", "#8b5cf6", "#ec4899", "#14b8a6", "#f97316", "#64748b"]
        return hues[hash(name) % len(hues)]

    # Build per-label phase bars for HTML detail
    per_label = phase_info["per_label"]
    # Compute max phase time for bar scaling globally
    max_phase_time = max(
        max([v for k,v in global_orig.items() if k not in PARENT_TIMERS] + [1]),
        max([v for k,v in global_new.items() if k not in PARENT_TIMERS] + [1]),
    )
    # Also per-label max for config table bars: use max_orig_mean as before

    html_parts: List[str] = []
    html_parts.append("""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Benchmark Report</title>
<style>
  :root{--bg:#0f1115;--panel:#161a22;--line:#2a303d;--text:#e6e8ec;--dim:#919aab;--accent:#4f9dff;--green:#3ddc97;--amber:#f6c453;--danger:#ff6b78;--cyan:#22d3ee;--purple:#b692ff}
  *{box-sizing:border-box} body{margin:0;color:var(--text);background:var(--bg);font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Helvetica,Arial,sans-serif;line-height:1.5;padding:24px}
  a{color:var(--accent)} h1,h2,h3{color:var(--text)} .wrap{max-width:1200px;margin:0 auto}
  .card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:18px;margin:16px 0}
  table{width:100%;border-collapse:collapse;font-size:13px} th,td{padding:8px 10px;border-bottom:1px solid var(--line);text-align:right;white-space:nowrap} th:first-child,td:first-child{text-align:left;white-space:normal} th{color:var(--dim);font-size:11px;text-transform:uppercase;letter-spacing:.05em;background:#12161d}
  .bar{height:12px;border-radius:4px;display:inline-block;vertical-align:middle}
  .bar-orig{background:var(--danger)} .bar-new{background:var(--green)}
  .bar-track{height:12px;background:#0a0d12;border-radius:4px;overflow:hidden;display:inline-block;vertical-align:middle;width:90px}
  .seg{height:100%;display:inline-block}
  .pill{padding:2px 8px;border-radius:99px;font-size:11px;font-weight:600}
  .pill-pass{background:rgba(61,220,151,.12);color:var(--green);border:1px solid #315c50}
  .pill-fail{background:rgba(255,107,120,.12);color:var(--danger);border:1px solid #5c2f33}
  .pill-skip{background:rgba(246,196,83,.12);color:var(--amber);border:1px solid #5c4a2c}
  .muted{color:var(--dim)} .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
  .grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px} @media(max-width:800px){.grid2{grid-template-columns:1fr}}
  .headline{font-size:28px;letter-spacing:-.02em} .headline small{color:var(--dim);font-size:14px;font-weight:400}
  .legend{display:flex;gap:12px;flex-wrap:wrap;margin:8px 0} .legend span{display:inline-flex;align-items:center;gap:6px;font-size:12px;color:var(--dim)} .legend i{width:12px;height:12px;border-radius:3px;display:inline-block}
</style>
</head>
<body><div class="wrap">""")
    html_parts.append(f"""<h1>Benchmark Report</h1>
<p class="muted">""")
    if run_url:
        html_parts.append(f'Run: <a href="{run_url}">{run_url}</a>')
    else:
        html_parts.append(f'Run ID: {run_id or "local"}')
    if sha:
        html_parts.append(f' · commit: <span class="mono">{sha[:8]}</span>')
    if docs:
        eps_example = docs[0].get("epsilon", "")
        html_parts.append(f' · {len(docs)} configurations · {len([d for d in docs if d.get("gated_count")])} gated')
    html_parts.append("</p>")

    # Headline card
    html_parts.append('<div class="card">')
    if speedups:
        html_parts.append(f'<div class="headline">{headline_speedup_range} <small>gated speedup (orig/new) · mean {mean_sp:.2f}× · {pass_cnt}/{len(docs)} passed</small></div>')
        html_parts.append(f'<p class="muted">Gated mean ratios (new/orig): {min(all_ratios):.3f}–{max(all_ratios):.3f}. Per-ID limit 1.20×, mean limit 1.05×.</p>')
    else:
        html_parts.append(f'<div class="headline">n/a <small>no gated configs</small></div>')
    html_parts.append("</div>")

    # Legend
    html_parts.append("""<div class="legend">
  <span><i style="background:#ff6b78"></i> orig ms</span>
  <span><i style="background:#3ddc97"></i> new ms</span>
  <span><i style="background:#4f9dff"></i> find_F</span>
  <span><i style="background:#3ddc97"></i> Gi hull</span>
  <span><i style="background:#f6c453"></i> boundary_P</span>
  <span><i style="background:#b692ff"></i> Gi prep</span>
</div>""")

    # Overall table
    html_parts.append('<div class="card"><h2>Overall per configuration</h2><div style="overflow:auto"><table>')
    html_parts.append('<thead><tr><th>Configuration</th><th>ε</th><th>δ</th><th>orig ms<br><span class="muted">gated mean</span></th><th>new ms<br><span class="muted">gated mean</span></th><th>speedup<br><span class="muted">orig/new</span></th><th>gated</th><th>status</th><th>time bar<br><span class="muted">orig → new</span></th></tr></thead><tbody>')
    for doc in docs:
        label = doc.get("label","")
        eps = doc.get("epsilon","")
        delta = doc.get("delta","")
        mean_orig = doc.get("mean_orig_ms")
        mean_new = doc.get("mean_new_ms")
        mean_ratio = doc.get("mean_ratio")
        gated = doc.get("gated_count",0)
        overall_ok = doc.get("overall_ok", False)
        mean_ok = doc.get("mean_ok", True)
        # status pill
        if gated == 0:
            status_html = '<span class="pill pill-skip">SKIP</span>'
        elif overall_ok and mean_ok:
            status_html = '<span class="pill pill-pass">PASS</span>'
        else:
            status_html = '<span class="pill pill-fail">FAIL</span>'
        # speedup
        if isinstance(mean_ratio, (int,float)) and mean_ratio and mean_ratio != 0:
            speedup = 1/mean_ratio
            speed_s = f"{speedup:.2f}×"
            speedup_color = "color:#3ddc97" if speedup >1 else "color:#ff6b78" if speedup <1 else ""
        else:
            speed_s = "n/a"
            speedup_color = ""
        # bars
        orig_w = 0
        new_w = 0
        if isinstance(mean_orig, (int,float)) and isinstance(mean_new, (int,float)) and max_mean>0:
            orig_w = min(90, mean_orig / max_mean * 90)
            new_w = min(90, mean_new / max_mean * 90)
        elif gated==0:
            # show Σ totals maybe? but we have None -> 0
            orig_w = new_w = 0
        bar_html = f'<span class="bar-track"><span class="bar bar-orig" style="width:{orig_w:.1f}px"></span></span> <span class="muted">→</span> <span class="bar-track"><span class="bar bar-new" style="width:{new_w:.1f}px"></span></span>'
        # ms strings
        orig_s = f"{mean_orig:.3f}" if isinstance(mean_orig, (int,float)) else "n/a"
        new_s = f"{mean_new:.3f}" if isinstance(mean_new, (int,float)) else "n/a"
        delta_s = f"{delta:.6g}" if isinstance(delta, (int,float)) else str(delta)
        html_parts.append(f'<tr><td class="mono">{label}</td><td>{eps}</td><td>{delta_s}</td><td>{orig_s}</td><td>{new_s}</td><td style="{speedup_color}">{speed_s}</td><td>{gated}/10</td><td>{status_html}</td><td>{bar_html}</td></tr>')
    html_parts.append('</tbody></table></div>')
    html_parts.append('<p class="muted">Thresholds: per-ID new ≤ orig×1.20, mean(new) ≤ mean(orig)×1.05 (gated IDs only, orig ≥1.0 ms). Σ rows have no gated IDs.</p></div>')

    # Per-phase global
    html_parts.append('<div class="card"><h2>Per-phase before/after (global sum over IDs)</h2>')
    html_parts.append('<p class="muted">Aggregated from <span class="mono">TIMER_MS</span> lines. Each ID\'s phase time is one <span class="mono">--time</span> run; values summed over the 10 IDs per configuration then over all configurations. Share = phase / total.</p>')
    # Determine totals
    leaf_orig_global = sum(global_orig.get(k, 0) for k in leaf_phases)
    leaf_new_global = sum(global_new.get(k, 0) for k in leaf_phases)
    global_orig_total = global_orig.get("total", leaf_orig_global) or leaf_orig_global or 1
    global_new_total = global_new.get("total", leaf_new_global) or leaf_new_global or 1
    html_parts.append('<div style="overflow:auto"><table><thead><tr><th>Phase</th><th>orig ms</th><th>orig share</th><th>new ms</th><th>new share</th><th>speedup</th><th>share bar (orig → new)</th></tr></thead><tbody>')
    max_phase = max([global_orig.get(p,0) for p in leaf_phases] + [global_new.get(p,0) for p in leaf_phases] + [1])
    for phase in leaf_phases:
        o_ms = global_orig.get(phase, 0)
        n_ms = global_new.get(phase, 0)
        o_share = o_ms / global_orig_total * 100 if global_orig_total else 0
        n_share = n_ms / global_new_total * 100 if global_new_total else 0
        sp = (o_ms / n_ms) if n_ms else float('inf')
        sp_s = f"{sp:.2f}×" if sp != float('inf') else "∞"
        col = phase_color(phase)
        o_w = o_ms / max_phase * 80
        n_w = n_ms / max_phase * 80
        bar = f'<span class="bar" style="background:{col};width:{o_w:.1f}px"></span> → <span class="bar" style="background:{col};opacity:.55;width:{n_w:.1f}px"></span>'
        html_parts.append(f'<tr><td class="mono">{phase}</td><td>{o_ms:.2f}</td><td>{o_share:.1f}%</td><td>{n_ms:.2f}</td><td>{n_share:.1f}%</td><td>{sp_s}</td><td>{bar}</td></tr>')
    # other
    if global_orig_total > leaf_orig_global or global_new_total > leaf_new_global:
        o_other = max(0, global_orig_total - leaf_orig_global)
        n_other = max(0, global_new_total - leaf_new_global)
        o_share = o_other / global_orig_total * 100 if global_orig_total else 0
        n_share = n_other / global_new_total * 100 if global_new_total else 0
        sp = (o_other / n_other) if n_other else float('inf')
        sp_s = f"{sp:.2f}×" if sp != float('inf') else "∞"
        o_w = o_other / max_phase * 80 if max_phase else 0
        n_w = n_other / max_phase * 80 if max_phase else 0
        bar = f'<span class="bar" style="background:#64748b;width:{o_w:.1f}px"></span> → <span class="bar" style="background:#64748b;opacity:.55;width:{n_w:.1f}px"></span>'
        html_parts.append(f'<tr><td class="mono"><em>other</em></td><td>{o_other:.2f}</td><td>{o_share:.1f}%</td><td>{n_other:.2f}</td><td>{n_share:.1f}%</td><td>{sp_s}</td><td>{bar}</td></tr>')
    html_parts.append(f'<tr style="font-weight:700;background:#12161d"><td>total</td><td>{global_orig_total:.2f}</td><td>100%</td><td>{global_new_total:.2f}</td><td>100%</td><td>{global_orig_total/(global_new_total) if global_new_total else 0:.2f}×</td><td></td></tr>')
    html_parts.append('</tbody></table></div></div>')

    # Per-configuration phase details (expandable)
    if per_label and any(pl["orig_totals"] or pl["new_totals"] for pl in per_label):
        html_parts.append('<div class="card"><h2>Per-configuration phase details</h2>')
        html_parts.append('<p class="muted">Each configuration\'s phase times summed over its 10 IDs. Bars show orig (solid) → new (faded). Share = phase / total for that configuration.</p>')
        for pl in per_label:
            label = pl["label"]
            orig_totals = pl["orig_totals"]
            new_totals = pl["new_totals"]
            orig_total = pl["orig_total_time"] or 1
            new_total = pl["new_total_time"] or 1
            doc = pl["doc"]
            # Determine phases present in this label
            phases_here = sorted(set(list(orig_totals.keys()) + list(new_totals.keys())), key=lambda x: (0 if x not in PARENT_TIMERS else 99, x))
            phases_here = [p for p in phases_here if p not in PARENT_TIMERS]
            if not phases_here:
                continue
            html_parts.append(f'<h3 class="mono">{label} <span class="muted">ε={doc.get("epsilon")} δ={doc.get("delta",0):.6g}</span></h3>')
            html_parts.append('<div style="overflow:auto"><table><thead><tr><th>Phase</th><th>orig ms</th><th>orig %</th><th>new ms</th><th>new %</th><th>speedup</th><th>bar</th></tr></thead><tbody>')
            max_here = max(max(orig_totals.values(), default=1), max(new_totals.values(), default=1))
            for ph in phases_here:
                o_ms = orig_totals.get(ph, 0)
                n_ms = new_totals.get(ph, 0)
                o_share = o_ms / orig_total * 100 if orig_total else 0
                n_share = n_ms / new_total * 100 if new_total else 0
                sp = (o_ms / n_ms) if n_ms else (float('inf') if o_ms else 0)
                sp_s = f"{sp:.2f}×" if sp not in (float('inf'),0) else ("∞" if sp==float('inf') else "n/a")
                col = phase_color(ph)
                o_w = o_ms / max_here * 70 if max_here else 0
                n_w = n_ms / max_here * 70 if max_here else 0
                bar = f'<span class="bar" style="background:{col};width:{o_w:.1f}px"></span> → <span class="bar" style="background:{col};opacity:.5;width:{n_w:.1f}px"></span>'
                html_parts.append(f'<tr><td class="mono">{ph}</td><td>{o_ms:.2f}</td><td>{o_share:.1f}%</td><td>{n_ms:.2f}</td><td>{n_share:.1f}%</td><td>{sp_s}</td><td>{bar}</td></tr>')
            # Totals row
            html_parts.append(f'<tr style="font-weight:700"><td>total</td><td>{orig_total:.2f}</td><td>100%</td><td>{new_total:.2f}</td><td>100%</td><td>{orig_total/(new_total) if new_total else 0:.2f}×</td><td></td></tr>')
            html_parts.append('</tbody></table></div>')
        html_parts.append('</div>')

    # Detailed per-ID table collapsible? Could add but keep concise
    html_parts.append('<div class="card"><h2>Per-ID details</h2><p class="muted">Gated IDs contribute to mean gate; others are below 1 ms floor. Ratio = new/orig.</p>')
    for doc in docs:
        label = doc.get("label","")
        cases = doc.get("cases", [])
        if not cases:
            continue
        html_parts.append(f'<h3 class="mono">{label}</h3>')
        html_parts.append('<div style="overflow:auto"><table><thead><tr><th>ID</th><th>orig ms</th><th>new ms</th><th>ratio (new/orig)</th><th>speedup</th><th>gated</th><th>status</th></tr></thead><tbody>')
        for c in sorted(cases, key=lambda x: x.get("id")):
            cid = c.get("id")
            o_ms = c.get("orig_ms")
            n_ms = c.get("new_ms")
            ratio = c.get("ratio")
            gated = c.get("gated")
            status = c.get("status","")
            if isinstance(ratio, (int,float)) and ratio:
                speed = 1/ratio
                speed_s = f"{speed:.2f}×"
            else:
                speed_s = "n/a"
            ratio_s = f"{ratio:.3f}" if isinstance(ratio, (int,float)) else "n/a"
            pill = '<span class="pill pill-pass">OK</span>' if status=="OK" else '<span class="pill pill-fail">FAIL</span>' if "FAIL" in status else '<span class="pill pill-skip">SKIP</span>'
            if status=="SKIP_FLOOR":
                pill='<span class="pill pill-skip">SKIP</span>'
            html_parts.append(f'<tr><td>{cid}</td><td>{o_ms:.3f}</td><td>{n_ms:.3f}</td><td>{ratio_s}</td><td>{speed_s}</td><td>{"yes" if gated else "no"}</td><td>{pill}</td></tr>')
        html_parts.append('</tbody></table></div>')
    html_parts.append('</div>')

    html_parts.append("""
<div class="card muted" style="font-size:12px">
Generated by <span class="mono">scripts/ci/bench_report.py</span> from benchmark artifacts. Thresholds: per-ID new ≤ orig×1.20, mean(new) ≤ mean(orig)×1.05 (gated IDs only). TIMER_MS phase times are from one <span class="mono">--time</span> run per ID and summed per configuration.
</div>""")
    html_parts.append("</div></body></html>")
    return "\n".join(html_parts)

def render_comment(docs: List[Dict[str, Any]], phase_info: Dict[str, Any], run_url: str) -> str:
    lines: List[str] = []
    lines.append("<!-- bench-report -->")
    lines.append("## Benchmark Report")
    if run_url:
        lines.append(f"Run: {run_url}")
    # Headline speedups
    gated_docs = [d for d in docs if d.get("gated_count",0)>0]
    ratios = [d.get("mean_ratio") for d in gated_docs if isinstance(d.get("mean_ratio"), (int,float))]
    if ratios:
        speedups = [1/r for r in ratios if r]
        lines.append(f"**Speedup (orig/new) on gated means: {min(speedups):.2f}–{max(speedups):.2f}× (mean {sum(speedups)/len(speedups):.2f}×)** — {len(gated_docs)}/{len(docs)} configurations gated.")
    # Pass/fail
    pass_cnt = sum(1 for d in docs if d.get("overall_ok"))
    lines.append(f"Gating: **{pass_cnt}/{len(docs)} passed** (per-ID ≤1.20×, mean ≤1.05×).")
    lines.append("")
    # Compact overall table (markdown)
    lines.append("| Configuration | orig ms | new ms | speedup | status |")
    lines.append("|---|---:|---:|---:|---|")
    for doc in docs:
        label = doc.get("label","")
        mean_ratio = doc.get("mean_ratio")
        mean_orig = doc.get("mean_orig_ms")
        mean_new = doc.get("mean_new_ms")
        gated = doc.get("gated_count",0)
        overall_ok = doc.get("overall_ok", False)
        if isinstance(mean_ratio, (int,float)) and mean_ratio:
            sp = 1/mean_ratio
            sp_s = f"{sp:.2f}×"
        else:
            sp_s = "n/a"
        orig_s = fmt_ms(mean_orig) if gated else "Σ"
        new_s = fmt_ms(mean_new) if gated else "Σ"
        status = "✅" if (gated==0 or overall_ok) else "❌"
        # shorten label for comment: remove extra prefix? Keep full
        lines.append(f"| {label} | {orig_s} | {new_s} | {sp_s} | {status} |")
    lines.append("")
    # Phase headline if available
    leaf_phases = phase_info.get("leaf_phases", [])
    global_orig = phase_info.get("global_orig", {})
    global_new = phase_info.get("global_new", {})
    if leaf_phases and global_orig and global_new:
        leaf_orig_sum = sum(global_orig.get(p,0) for p in leaf_phases)
        leaf_new_sum = sum(global_new.get(p,0) for p in leaf_phases)
        total_o = global_orig.get("total", leaf_orig_sum) or leaf_orig_sum
        total_n = global_new.get("total", leaf_new_sum) or leaf_new_sum
        lines.append("**Per-phase (global sum):**")
        lines.append("")
        # show top phases
        for ph in leaf_phases[:4]:
            o = global_orig.get(ph,0)
            n = global_new.get(ph,0)
            sp = o/n if n else float('inf')
            sp_s = f"{sp:.2f}×" if sp != float('inf') else "∞"
            lines.append(f"- {ph}: {o:.1f} → {n:.1f} ms ({sp_s}, {o/total_o*100:.1f}% → {n/total_n*100:.1f}% share)")
        lines.append("")
    lines.append("_Full tables and inline charts are in the `benchmark-report` artifact (`report.html`). See job summary for per-configuration and per-phase tables._")
    if run_url:
        lines.append(f"_Artifact: {run_url} (download `benchmark-report`)_")
    lines.append("")
    return "\n".join(lines)

def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate benchmark.json artifacts into report")
    parser.add_argument("--reports-dir", type=Path, required=True, help="Directory containing benchmark.json artifacts (recursively searched)")
    parser.add_argument("--output-md", type=Path, help="Write markdown summary to file (also stdout if not given)")
    parser.add_argument("--output-html", type=Path, help="Write self-contained HTML report")
    parser.add_argument("--output-comment", type=Path, help="Write sticky PR comment markdown")
    parser.add_argument("--run-url", type=str, default=os.environ.get("GITHUB_RUN_URL", ""), help="URL to link in report (e.g. https://github.com/org/repo/actions/runs/123)")
    parser.add_argument("--run-id", type=str, default=os.environ.get("GITHUB_RUN_ID", ""), help="Run ID if URL not available")
    parser.add_argument("--sha", type=str, default=os.environ.get("GITHUB_SHA", ""), help="Commit SHA for header")
    parser.add_argument("--check", action="store_true", help="Exit 1 if any configuration failed gating")
    args = parser.parse_args()

    reports_dir: Path = args.reports_dir
    if not reports_dir.exists():
        print(f"error: reports dir {reports_dir} does not exist", file=sys.stderr)
        return 2
    paths = collect_reports(reports_dir)
    if not paths:
        print(f"warn: no benchmark.json found under {reports_dir}", file=sys.stderr)
        # Still produce empty report placeholders
        docs: List[Dict[str, Any]] = []
        phase_info = {"per_label": [], "global_orig": {}, "global_new": {}, "all_phases": [], "leaf_phases": []}
    else:
        docs = []
        for p in paths:
            d = load_report(p)
            if d is not None:
                docs.append(d)
        if not docs:
            print("warn: no docs loaded", file=sys.stderr)
            phase_info = {"per_label": [], "global_orig": {}, "global_new": {}, "all_phases": [], "leaf_phases": []}
        else:
            docs = aggregate_configs(docs)
            phase_info = compute_phase_totals(docs)

    # Auto-derive run_url if not given but env has repo and run_id
    run_url = args.run_url
    if not run_url and args.run_id:
        repo = os.environ.get("GITHUB_REPOSITORY", "")
        if repo:
            run_url = f"https://github.com/{repo}/actions/runs/{args.run_id}"

    md = render_markdown(docs, phase_info, run_url, args.run_id)
    html = render_html(docs, phase_info, run_url, args.run_id, args.sha)
    comment = render_comment(docs, phase_info, run_url)

    if args.output_md:
        args.output_md.parent.mkdir(parents=True, exist_ok=True)
        args.output_md.write_text(md)
        print(f"wrote markdown to {args.output_md}")
    else:
        print(md)

    if args.output_html:
        args.output_html.parent.mkdir(parents=True, exist_ok=True)
        args.output_html.write_text(html)
        print(f"wrote html to {args.output_html}")

    if args.output_comment:
        args.output_comment.parent.mkdir(parents=True, exist_ok=True)
        args.output_comment.write_text(comment)
        print(f"wrote comment to {args.output_comment}")

    # Also print to stdout if files given? already printed md if no output_md

    if args.check:
        failed = [d for d in docs if not d.get("overall_ok", False) and d.get("gated_count",0)>0]
        if failed:
            print(f"check failed: {len(failed)} configurations failed gating", file=sys.stderr)
            return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
