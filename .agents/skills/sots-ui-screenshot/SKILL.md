---
name: sots-ui-screenshot
description: Capture before/after screenshots at laptop and mobile for the web viewer. Use whenever a PR changes the web viewer's UI, layout, styling or behavior.
---

# SOTS UI Screenshot

Applies when a PR changes the web viewer's UI, layout, styling or behavior.

## Requirement

Capture before and after screenshots at both sizes:

- Laptop 1280x800
- Mobile 390x844

Use the same trajectory, URL, viewport and scroll per pair. One pair per size — total 4 images (before/after × laptop/mobile).

A UI PR with only one size is incomplete.

## Capture

Use the isolated automation browser for viewing only. Never log in, never read Keychain/cookies/credentials or any stored credentials. The browser is for viewing only — do not reuse a logged-in profile or decrypt stored credentials.

## Attach

Upload via `gh --attach` (gh 2.100.0, up to 50 files per invocation) so images render inline in the PR description as `https://github.com/user-attachments/assets/...` URLs. Keep the bare URL or `![...](url)` markdown that `gh` inserts verbatim — do not rewrite to `releases/download/...`.

```bash
gh pr edit <number> --attach /path/to/laptop-before.png --attach /path/to/laptop-after.png --attach /path/to/mobile-before.png --attach /path/to/mobile-after.png
# or per-comment:
gh pr comment <number> --attach /path/to/image.png --body "context"
```

If `gh --attach` is unavailable, fall back to web UI drag-and-drop, but `gh --attach` is preferred.
