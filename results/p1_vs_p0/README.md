# P1 vs P0 performance evidence (IDs 1-50)

Compares **P0 tip** (wedge prune only, commit `f4123f3`) against **P0+P1**
(Gi dedup once per stab step + skip Sutherland-Hodgman `to_ccw`).

Not a main-vs-branch CI gate: this isolates the incremental cost of P1.

## Method

- Taxi trajectories 1-50 (T-Drive normalized)
- Median of N timed runs of `SIMPLIFY_CORE_MS` per ID
- Ratio = opt_ms / base_ms (`<1` means P1 faster)

## Results

### e=0.5, d=200 (median of 3)

| Metric | Value |
|---|---|
| Median ratio | 0.938 |
| Mean / geomean | 0.935 / 0.927 |
| Faster / slower / flat | 64% / 24% / 12% |
| Total core time | 4387 ms -> 4023 ms (**-8.3%**, saved 364 ms) |

### e=299, d=1 CI params (median of 5)

| Metric | Value |
|---|---|
| Median ratio | 0.959 |
| Mean / geomean | 0.962 / 0.958 |
| Faster / slower / flat | 68% / 18% / 14% |
| Total core time | 701 ms -> 654 ms (**-6.8%**, saved 47 ms) |

Plots and CSVs live next to this file under `results/p1_vs_p0/` and
`results/p1_vs_p0_ci5/`.

## Reading the variance

P1 is a small constant-factor win on the intersect path (one Gi dedup per
stream step instead of per alive candidate; skip orientation fix). Absolute
savings are a few percent, so short trajectories and thermal noise can flip
individual ratios above 1.0. Aggregate wall time and median ratio still favor
P1 on both parameter settings.
