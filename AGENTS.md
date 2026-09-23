# Project agent memory

This file is the project's committed home for project-intrinsic agent knowledge: build, test, release, architecture, and sharp-edge notes that should travel with the code.

- Add durable project-specific notes here as they are discovered through real work.

## CI gates and timing

- Correctness and Benchmark workflows share a 20-entry `(ε, δ, size)` matrix (5 ε tiers ×2 δ formulas 300/(1+ε) & 1000/(1+ε) ×2 sizes small ~100 pts IDs 11..20 / large ~1000 pts IDs 21..30, derived via `scripts/derive_benchmark_data.py`); bars and artifact layout are documented in `.github/workflows/README.md`. Their `orig` checkout is the PR base (`github.event.pull_request.base.sha`); `base_sha` is not a context field and silently falls back to the previous push.
- Opt-in phase timers: `./build/simplify <id> --time` writes a human TIMING SUMMARY plus machine `TIMER_MS <name> <ms> <calls>` lines on stderr (`timer.h`, wired from `simplify.cpp`).
- Local `SIMPLIFY_CORE_MS` means can swing under host load. Interleave main and candidate runs, and inspect per-ID minima alongside the five-run means used by CI.

## Maintaining this file

Keep this file for knowledge useful to almost every future agent session in this project.
Do not repeat what the codebase already shows; point to the authoritative file or command instead.
Prefer rewriting or pruning existing entries over appending new ones.
When updating this file, preserve this bar for all agents and keep entries concise.
