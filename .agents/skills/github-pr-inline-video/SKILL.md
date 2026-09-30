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

### Preferred: web UI drag-and-drop (always works)

1. Edit the PR description or add a new comment.
2. Drag-and-drop the `.mp4` into the markdown editor (or click to attach) and wait for upload to complete.
3. Keep the markdown GitHub inserts verbatim, e.g.:

   ```markdown
   https://github.com/user-attachments/assets/<uuid>
   ```

   or

   ```markdown
   ![description](https://github.com/user-attachments/assets/<uuid>)
   ```

   Either form renders as an inline playable video. Do not rewrite it to a `releases/download/...` URL.

### Scripted alternative (verify before use)

There is **no stable single-command `gh` attachment upload**. To automate, you must upload via the GitHub API to obtain a `user-attachments/assets` URL and then reference that URL in a PR comment. Any `gh api` upload invocation must be verified against the current GitHub API documentation before use — do not invent or assume an endpoint (e.g., do not assume `POST /repos/{owner}/{repo}/attachments` exists without verifying it). The web UI above is the guaranteed path.

Once you have the verified `user-attachments/assets` URL, add it to the PR:

```bash
gh pr comment <number> --body "![description](https://github.com/user-attachments/assets/<uuid>)"
# or via the API:
gh api repos/OWNER/REPO/issues/<number>/comments -f body="![description](https://github.com/user-attachments/assets/<uuid>)"
```

Never use a `releases/download/...` URL for PR-inline video.

## Evidence and motivating case

Project evidence: `docs/pr-evidence.md` documents this rule and the full attachment mechanics — refer to it for details.

Motivating case: PR #28 video comment https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/pull/28#issuecomment-5894279277 published `sots-pr28-viewer-video.mp4` as a release asset (`releases/download/...`) and referenced it with both a plain link and `![](...releases/download/...mp4)`. Neither rendered inline. The same file would have played inline if attached directly to the PR as a `user-attachments/assets` upload.
