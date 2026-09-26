# Project agent memory

This file is the project's committed home for project-intrinsic agent knowledge: build, test, release, architecture, and sharp-edge notes that should travel with the code.

- Add durable project-specific notes here as they are discovered through real work.

## CI gates and timing

- Correctness and Benchmark workflows share a 20-entry `(ε, δ, size)` matrix (5 ε tiers ×2 δ formulas 300/(1+ε) & 1000/(1+ε) ×2 sizes small ~100 pts IDs 11..20 / large ~1000 pts IDs 21..30, derived via `scripts/derive_benchmark_data.py`); bars and artifact layout are documented in `.github/workflows/README.md`. Their `orig` checkout is the PR base (`github.event.pull_request.base.sha`); `base_sha` is not a context field and silently falls back to the previous push.
- Benchmark adds an aggregated report: `bench-report` job collects all `benchmark-<label>/benchmark.json` (now with `orig_ops`/`new_ops` from `TIMER_MS`), renders a per-configuration overall table and per-phase before/after table with time shares via `scripts/ci/bench_report.py` (usable locally on saved logs), writes the same tables to `GITHUB_STEP_SUMMARY` and to a self-contained `benchmark-report/report.html` artifact, and posts a single sticky PR comment (`<!-- bench-report -->`). See `.github/workflows/README.md` for the invocation.
- Opt-in phase timers: `./build/simplify <id> --time` writes a human TIMING SUMMARY plus machine `TIMER_MS <name> <ms> <calls>` lines on stderr (`timer.h`, wired from `simplify.cpp`).

## Maintaining this file

Keep this file for knowledge useful to almost every future agent session in this project.
Do not repeat what the codebase already shows; point to the authoritative file or command instead.
Prefer rewriting or pruning existing entries over appending new ones.
When updating this file, preserve this bar for all agents and keep entries concise.
