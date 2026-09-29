# PR video evidence — inline playback rule

Videos render and play inline inside a GitHub PR only when uploaded as
file attachments in the PR description or a PR comment. GitHub then serves
the file from `https://github.com/user-attachments/assets/...` and renders
an inline player.

## Why external URLs and `<video>` fail

- Raw HTML `<video>` tags are sanitized/stripped by GitHub markdown and do
  not play inline.
- Standard markdown image/video syntax pointing to external URLs — including
  release assets such as
  `https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/releases/download/...`
  — is blocked or shown as a plain link/download. GitHub does not embed
  release-asset `.mp4` links as an inline player, even when the URL ends in
  `.mp4` or is wrapped in `![](...)`.
- Only the `user-attachments/assets` flow is treated as an inline video
  attachment.

## Correct attachment flow

Upload the `.mp4` file as an attachment in the PR description or a PR
comment so GitHub creates a `user-attachments/assets` URL. Do not use
release assets for PR-inline video.

### Web UI (simplest)

1. Edit the PR description or add a new comment.
2. Drag-and-drop the `.mp4` file (or click the attachment area) into the
   markdown editor and wait for the upload to finish.
3. GitHub inserts markdown like:

   ```markdown
   https://github.com/user-attachments/assets/<uuid>
   ```

   or

   ```markdown
   ![viewer walkthrough](https://github.com/user-attachments/assets/<uuid>)
   ```

   Either form renders as an inline playable video — keep the markdown
   GitHub generated, do not rewrite it to a release-asset URL.

### Scripted with `gh` (when automation is needed)

There is no stable `gh` subcommand that uploads an attachment in one step;
the attachment must be created via the GitHub attachments API and then
referenced in a comment. The pattern is:

```bash
# 1. Upload the video as an attachment (creates a user-attachments URL).
#    The endpoint is the same one the web UI uses; gh wraps it:
gh api --method POST \
  -H "Accept: application/vnd.github+json" \
  -H "Content-Type: video/mp4" \
  --input sots-pr28-viewer-video.mp4 \
  https://api.github.com/repos/yeungsinchun/Simplification-of-Trajectory-Streams/attachments
# Response includes the https://github.com/user-attachments/assets/... URL

# 2. Reference that URL in a PR comment — this is what renders inline:
gh pr comment 28 --body "![viewer walkthrough](https://github.com/user-attachments/assets/<uuid>)"
# or via the API:
gh api repos/yeungsinchun/Simplification-of-Trajectory-Streams/issues/28/comments \
  -f body="![viewer walkthrough](https://github.com/user-attachments/assets/<uuid>)"
```

If the API upload step is unavailable in your `gh` version, use the web
drag-and-drop flow above and then use `gh pr comment` / `gh api` only for
the markdown reference. Never substitute a release-asset URL such as
`https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/releases/download/sots-pr28-viewer-video-20250930/sots-pr28-viewer-video.mp4` —
that URL does not embed inline.

## Motivating case — PR #28

PR #28 (https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/pull/28)
added viewer trace-diet and render-cache work. Its video evidence
(`sots-pr28-viewer-video.mp4`, 45 s, 1280×720) was published as a release
asset at
`https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/releases/download/sots-pr28-viewer-video-20250930/sots-pr28-viewer-video.mp4`
and referenced in
https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/pull/28#issuecomment-5894279277
with both a plain link and `![](...release...mp4)`. Neither rendered as an
inline player. The same file would have played inline if uploaded as a PR
attachment (`user-attachments/assets`).

Next time: record the video as before, but attach the `.mp4` directly to
the PR description or a PR comment (drag-and-drop or the `gh api` attachment
upload above) instead of creating a release asset. Keep the
`user-attachments/assets` markdown GitHub generates.
