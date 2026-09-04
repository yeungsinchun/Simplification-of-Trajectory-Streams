---
name: algo-pr-description
description: >-
  Write clear algorithm-change PR descriptions for this trajectory-simplification
  repo, including after no-mistakes opens a PR. Use when creating or editing pull
  requests, drafting PR bodies, driving no-mistakes axi run/respond, or when the
  user asks for PR description quality. Requires a short numbered stab-loop
  walkthrough, not an intent dump.
disable-model-invocation: false
---

# Algorithm PR descriptions

When making a PR (including via no-mistakes or `gh` / `gh-axi`), describe the
**algorithm change** clearly. Do not paste a raw `--intent` paragraph as the
main readable story.

## Required shape (human-facing PR body)

Lead with a short explanation a reviewer can follow without opening the diff:

1. **What is being changed** - one or two sentences (e.g. which candidates die
   early, which geometry step is skipped).
2. **Context in the stab loop** - numbered steps for the relevant prefix logic,
   for example:
   1. `P` = boundary anchors of the δ-disk around `p0` (candidates for the prefix).
   2. Each index `i` keeps a stab region `S[i]`.
   3. At each new stream point, build `Gi = conv(G)` and try `F(S[i], P[i]) ∩ Gi`.
   4. If that intersection is empty, candidate `i` dies (`dead[i] = true`).
3. **What this PR adds/changes inside that loop** - how the new path relates to
   those steps (what it skips, what it reuses, what it must not change).
4. **Correctness stance** - Frechet-safe / identical outputs / conservative
   prune / etc., in one short paragraph.
5. **Performance** - honest summary: typical speedup, and whether some IDs or
   parameter regimes can still look flat or slightly slower (noise, miss-path
   cost, CI 1.5x bar).
6. **Test plan** - Frechet, points, perf gates (IDs 1-10, CI `e`/`d` if used).

Optional collapsed `<details>` for long logs are fine; they must not replace the
readable algo section above.

## Do not

- Lead with a wall of intent / acceptance-criteria prose.
- Say only "enable X" without saying **what X prunes or skips** (candidates in
  `P`, stream points, polygons, …).
- Claim uniform speedup if measurements showed mixed ratios.

## Verbatim rule from the project owner

When making PR, must describe the algo change clearly as in the numbered stab-loop walkthrough above (what is pruned, the 4 steps, how the new path fits), not a huge block of intent text.

## Instructing no-mistakes (important)

no-mistakes **hardcodes** PR body assembly (see its PR step docs):

- Prepends `## Intent` from `--intent`
- Lets the PR agent write only a short `## What Changed` (1-3 bullets)
- Appends `## Risk Assessment`, `## Testing`, `## Pipeline`

There is **no** `.no-mistakes.yaml` knob to replace that template. To get
descriptions like PR #6 anyway:

### 1) Pass a structured `--intent` (feeds `## Intent`)

When starting `no-mistakes axi run --intent "..."`, write the intent as the
numbered algo walkthrough (the Required shape sections 1-3), not a single dense
paragraph of acceptance criteria. Keep constraints short at the end
(e.g. "No Frechet/points/perf regression on IDs 1-10 at e=299 d=1.").

### 2) Rewrite the PR body after the PR step (required)

As soon as no-mistakes reports a PR URL (or on `checks-passed`), rewrite the
forge PR description with `gh-axi pr edit <n> --body-file ...` (or `gh pr edit`)
into the **Required shape** above.

Preserve any `<!-- no-mistakes-pipeline-attestation:v1 ... -->` HTML comment if
present. You may keep a short `## Pipeline` note plus attestation; do not leave
the default Intent-dump body as the only description.

### 3) Trusted repo review guidance

`.no-mistakes.yaml` on the **default branch** may carry `review.path_instructions`
so Review also expects algo clarity. That steers review findings; it does not
change the PR body template. PR readability still depends on steps 1-2.

## Example lead (wedge prune)

```markdown
## Summary

Early-out for **alive candidates** in `P` (boundary anchors for the current
prefix), not input stream points.

For each prefix starting at `p0`:

1. `P` = boundary anchors of the δ-disk grid around `p0`.
2. Each index `i` keeps a stab region `S[i]`.
3. At each new stream point, build `Gi = conv(G)` and try `F(S[i], P[i]) ∩ Gi`.
4. If that intersection is empty, candidate `i` dies.

This PR runs `wedge_gi_disjoint` before step 3's expensive work: if `Gi` is
robustly outside the tangent cone of `S[i]` at `P[i]`, mark the candidate dead
without building `F` or clipping. Same outcome as an empty intersection.
Tangents are shared with `find_F` on a miss so we do not pay the scan twice.
```
