---
name: github-pr-inline-video
description: Attach videos to GitHub PRs so they play inline. Use when adding video evidence to a PR description or comment — GitHub only renders videos inline from user-attachments/assets uploads.
---

# GitHub PR Inline Video

## Invariant

GitHub plays video inline (embedded player) in a PR only when the video is uploaded as a file attachment in the PR description or a PR comment. GitHub then serves it from `https://github.com/user-attachments/assets/...` and renders an inline player.

The following do **not** play inline:

- Raw HTML `<video>` tags — stripped/sanitized by GitHub markdown.
- Markdown image syntax `![](...)` pointing to external URLs — including release assets such as `https://github.com/OWNER/REPO/releases/download/.../*.mp4` — rendered as a plain link or download, not an inline player.

Only the `user-attachments/assets` attachment flow renders inline.

## Correct flow

### Preferred scripted: `gh --attach` (verified gh 2.100.0)

`gh` 2.100.0 adds first-class attachment upload that produces real `user-attachments/assets` URLs without browser automation. Prefer this over manual drag-and-drop when `gh` is available:

```bash
# Append a video asset to the existing PR description (keeps current body):
gh pr edit <number> --attach /path/to/video.mp4

# Or create a new PR comment with an attachment:
gh pr comment <number> --attach /path/to/video.mp4 --body "context text"
gh issue comment <number> --attach /path/to/video.mp4 --body "context text"
```

Up to 50 files per invocation (`--attach` may be repeated). `gh` appends the `https://github.com/user-attachments/assets/<uuid>` markdown to the PR body/comment verbatim — keep it as-is (do not rewrite to `releases/download/...`). The existing body is preserved; `--attach` does not overwrite it.

**Video alt-text limitation (verified):** do not append `#alt text` to a video attachment (`video.mp4#something` fails with `cannot set alt text on video`). Upload the bare `.mp4`.

### Fallback: web UI drag-and-drop (always works)

If `gh --attach` is unavailable, edit the PR description or add a new comment, drag-and-drop the `.mp4` into the markdown editor (or click to attach), wait for upload, and keep the markdown GitHub inserts verbatim, e.g.:

```markdown
https://github.com/user-attachments/assets/<uuid>
```

or

```markdown
![description](https://github.com/user-attachments/assets/<uuid>)
```

Either form renders as an inline playable video. Do not rewrite it to a `releases/download/...` URL.

### What not to do

- Raw `gh api` upload to an assumed `POST /repos/{owner}/{repo}/attachments` or `POST https://github.com/upload/policies/assets` with a PAT does **not** work in this repo (verified: returns `404 Not Found` or `"Oh no"` HTML). Do not invent or assume such endpoints.
- Never use a `releases/download/...` URL for PR-inline video — it renders as a plain link, not an inline player.

### Prohibited: credential / Keychain / cookie reuse

Do **not** run `security find-generic-password -s "Chrome Safe Storage" -w`, do **not** read the macOS Keychain, `Chrome Safe Storage`, browser Cookies / `Cookies` SQLite, or any stored credential, do **not** decrypt a browser profile, and do **not** log in via copied cookies as another user. The `user-attachments/assets` upload must be via `gh --attach` or manual browser drag-and-drop as the logged-in PR author — not via stolen credentials.

Once you have the verified `user-attachments/assets` URL, add it to the PR (if you obtained it via a separate upload):

```bash
gh pr comment <number> --body "![description](https://github.com/user-attachments/assets/<uuid>)"
# or via the API:
gh api repos/OWNER/REPO/issues/<number>/comments -f body="![description](https://github.com/user-attachments/assets/<uuid>)"
```

Never use a `releases/download/...` URL for PR-inline video.

## Evidence and motivating case

Project evidence: `docs/pr-evidence.md` documents this rule and the full attachment mechanics — refer to it for details.

Motivating case: PR #28 video comment https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/pull/28#issuecomment-5894279277 published `sots-pr28-viewer-video.mp4` as a release asset (`releases/download/...`) and referenced it with both a plain link and `![](...releases/download/...mp4)`. Neither rendered inline. The same file would have played inline if attached directly to the PR as a `user-attachments/assets` upload.
