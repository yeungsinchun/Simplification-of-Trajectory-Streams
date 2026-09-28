---
name: sots-perf-report
description: How to write a performance report for a Simplification-of-Trajectory-Streams change (PR description, review board, or scout report). Use whenever a change claims a speedup or touches the hot path of get_longest_stab / find_F / intersect. Covers the required sections, how to measure each one, and the reconciliations reviewers always ask for.
---

# SOTS performance report

A speedup claim is not reviewable from the headline number alone. Every perf
report for this repo must let a reviewer answer five questions without asking:
**how much faster, where did the time go before and after, what each change
contributed, what the percentages are relative to, and is the output
unchanged.** The sections below are required, in this order.

## 1. Headline: before/after per CI setting

- Use the 20 `(ε, δ, size)` settings from `.github/workflows/README.md`:
  5 ε tiers (299, 30, 5, 0.5, 0.1) ×2 δ formulas `300/(1+ε)` & `1000/(1+ε)` ×2 sizes small ~100 pts (IDs 11–20) / large ~1000 pts (IDs 21–30), derived via `scripts/derive_benchmark_data.py`.
- Metric: `SIMPLIFY_CORE_MS` from an **untimed** Release build (no `--time`),
  per setting (10 IDs: 11–20 small or 21–30 large). Interleave base and head runs (ABAB…) and report medians, then the
  CI-style gated mean over all 10 IDs per setting as defined in `.github/workflows/README.md`.
- Table columns: setting · base ms · head ms · speedup.

## 2. Output identity

- State `N / 200` byte-identical outputs (10 IDs × 20 settings: IDs 11–20 small, 21–30 large), via
  `./build/simplify <id> -e <ε> -d <δ>` then `cmp data/<id>/simplify.txt`
  across builds. If output changes, say so first, not last.

## 3. Phase time share, BEFORE and AFTER (required chart)

Never show only the head's breakdown. Show paired bars per setting (base vs
head), with bar length proportional to core time, and per-phase
`base ms → head ms (base% → head%)` for the hot phases.

- Phases: boundary P · Gi hull · Gi prep · find_F · clip (F ∩ Gi) · caches ·
  loop/other. If a phase only exists on one side (e.g. a cache), show it as 0
  on the other side.
- If base does work *inside* another call (e.g. per-candidate Gi dedup/CCW
  inside `intersect`/`clip`), give it its own nested scope and subtract it
  from the parent, so the categories match across revisions.
- Instrument with a local steady-clock accumulator (fixed enum slots, no map
  lookups). `timer.h`'s `TIMER()` does a `std::map<std::string>` lookup per
  call and distorts hot phases; use `--time` only for coarse attribution.
- Remove the clock overhead: scale phase ms so they sum to the untimed
  `SIMPLIFY_CORE_MS`, and say so. Note that nested scopes inflate overhead
  more on the revision that has more of them.
- Call out share shifts that come purely from other phases shrinking (for
  example, clip's share rising while its absolute ms fell).

## 4. Per-optimization contribution (one-feature-off ablation)

- Disable one optimization at a time in a temporary build; report the median
  penalty on a small ID (ID1 mid) and a large one (ID10 fine-e). Show negative
  results (for example, a prune that costs more than it saves).
- Include operation counts that explain the size (stream steps vs candidate
  visits, clips performed vs pruned, cache hits).
- Ablations overlap; say they must not be summed.

## 5. Reconciliations reviewers will ask for

- **Share vs ablation.** A phase can be ~1–3% of head runtime and still be worth
  +30% when un-hoisted. Show the arithmetic: extra calls × (penalty / extra
  calls) = per-call cost; per-call cost × head calls = head share. Check the
  result against the clocked phase share.
- **Baseline of every %.** Say explicitly that percentages are of core
  time (`SIMPLIFY_CORE_MS`, the timed `get_longest_stab` loop), not process
  wall time. Give wall vs core for one small and one large input, since
  startup and I/O dominate small inputs.
- **Ablation vs real diff.** If an ablation rebuilds the *new* richer object,
  note that the base's true cost differed.

## 6. Headroom and "can it be removed?"

- For each piece of setup work, state why it exists and whether it is a no-op
  for this input shape, with evidence. Example: Gi is a translated
  `CGAL::convex_hull_2` template, so dedup and CCW are no-ops. If feasible,
  prototype the removal locally and report identity and timing.
- End by naming where the remaining time is (the largest head phases). An A/B
  swing larger than the phase you touched is code-layout noise; say so rather
  than claiming it.

## 7. Complexity

- State time and space complexity before and after, with symbols defined
  (n points, p anchors, v live-region vertices, g Gi vertices). Say whether the
  change is asymptotic or constant-factor.

## Hygiene

- Record the host, compiler, build type, commit SHAs, run counts, and
  aggregation (median vs mean) next to every table.
- Keep instrumentation local (scratch copies such as `git archive <sha> | tar -x`);
  never commit profiling scopes.
