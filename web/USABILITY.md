# Visualizer usability findings

Running log of website usability issues found while checking desktop and mobile. Each item is fixed in its own change when practical.

## Fixed

### Top-bar spinners were not where the numbers appear

**Where:** header params bar, `Computed Fréchet distance` and `Simplification time`.

**Problem:** While those metrics were loading, the spinner was rendered *before* the label. The numeric value later appeared *after* the label. Users looking at the empty value slot saw no progress indicator, and swapping spinner-for-number also shifted neighboring params.

The Fréchet metric was also omitted until simplification finished, so the whole params row jumped when that chip appeared.

**Fix:** The spinner now lives in a reserved value slot after the label. Both metrics stay visible for the life of a loaded (or loading) trace, so the slot is replaced in place instead of the chip appearing later.

**Evidence:**

- Desktop (1280px): `web/usability/desktop-params-metrics.png`
- Mobile (390px): `web/usability/mobile-params-metrics.png`
- Markup fixture used for those screenshots: `web/usability/params-metric-states.html`

### Mobile loaded header hid Fréchet distance and simplification time

**Where:** mobile trace screen (`max-width: 720px`) after `body.trace-loaded-mobile`.

**Problem:** The compact loaded header hid all of `.headerBody`, and `#paramsBar` lives inside that node. Computed Fréchet distance and Simplification time were in the DOM with in-slot spinners, but users only saw Back and the title.

**Fix:** The loaded header still hides file controls and start instructions, but keeps `#paramsBar` as a second row under Back / title. On mobile that row shows only Computed Fréchet distance and Simplification time, so the metrics stay on screen while they load and after values arrive.

**Evidence:**

- Mobile (390px): `web/usability/mobile-loaded-params.png`
- Markup fixture used for that screenshot: `web/usability/mobile-loaded-params.html`

### Desktop drop-hint heading was misspelled

**Where:** empty-state `#dropHint` panel, `.desktop-instructions h2` (desktop only; the panel is hidden below 720px).

**Problem:** The heading was written as `Intruction`. The stylesheet uppercases it (`text-transform: uppercase`), so the first thing a desktop visitor saw in the instructions card was `INTRUCTION`. Mobile start-help and the sidebar already used `Instructions`.

**Fix:** The heading is now `Instructions`, so the card reads `INSTRUCTIONS`.

**Evidence:**

- Desktop (1280px): `web/usability/desktop-drop-hint-heading.png`
- Markup fixture used for that screenshot: `web/usability/desktop-drop-hint-heading.html`

### Mobile start-screen δ value was clipped

**Where:** 390px start form, `#deltaInput` in the two-column ε / δ row.

**Problem:** Each param label was a horizontal flex row (symbol + `input[type=number]`) with `min-width: 0`. Chrome's spinner arrows kept their width, so the default δ value `500` was clipped. ε `0.9` often still fit because it is narrower.

**Fix:** On viewports max-width 720px each param stacks the symbol above a full-column number field and hides the unused spinner arrows (`appearance: textfield`). Desktop stays an inline 62px field.

**Evidence:**

- Mobile live start screen (390px): `web/usability/mobile-delta-field-live.png`
- Before/after markup fixture (390px): `web/usability/mobile-delta-field.png`
- Markup fixture used for that screenshot: `web/usability/mobile-delta-field.html`

### Canvas stayed blank while a trace was computing

**Where:** `#traceLoadingHud` over `#canvasContainer`, after Load Trace and before stream geometry arrives.

**Problem:** The HUD existed in markup but CSS forced `#traceLoadingHud { display: none !important; }`, and viewer logic never toggled it. Load Trace hid `#dropHint` immediately, so desktop and the mobile loaded canvas were empty until the first header/geometry message. Progress lived only in the header (`Computing trace…` / metric slots), not where the trajectory would appear.

**Fix:** The HUD is a centered overlay on the canvas. It is shown from `showTraceLoading()` and hidden when geometry arrives (`initTraceUI`) or the load fails/finishes. A failed load also restores the drop hint so desktop is not left on a blank canvas.

**Evidence:**

- Desktop live loading (1280px): `web/usability/desktop-canvas-loading-hud-live.png`
- Mobile live loading (390px): `web/usability/mobile-canvas-loading-hud-live.png`
- Before/after markup fixture: `web/usability/canvas-loading-hud.html`
- Fixture screenshot (1280px): `web/usability/canvas-loading-hud.png`

### Hidden native select caused 390px start-screen overflow

**Where:** 390px start form, `#traceSelect.visually-hidden` in `.preloaded-row`.

**Problem:** Mobile hides the native `<select>` with `.visually-hidden` and shows a custom Preloaded trigger. `.upload-form select { width: 100% }` beat `.visually-hidden { width: 1px }`, so the control was about 390px wide starting at x=26. `documentElement.scrollWidth` was 416, and the start screen could pan sideways.

**Fix:** `.visually-hidden` box metrics now use `!important`, and the mobile form width rule skips `.visually-hidden` selects. The native control stays in the accessibility tree at 1px. Desktop (`min-width: 721px`) still shows the native select.

**Evidence:**

- Mobile live start screen (390px): `web/usability/mobile-select-overflow-live.png`
- Before/after markup fixture (390px): `web/usability/mobile-select-overflow.png`
- Markup fixture used for that screenshot: `web/usability/mobile-select-overflow.html`

### Mobile playback dock appeared before canvas geometry

**Where:** 390px loaded-trace screen, `#mobileTransport`, from Load Trace until the first stream header/geometry message.

**Problem:** `showTraceLoading()` entered the mobile loaded layout immediately (`body.trace-loaded-mobile`), and that class showed the bottom playback dock. Buttons were disabled, but Play still used the accent color over an empty canvas. Desktop `#playbackBar` already waited for `initTraceUI`; mobile did not.

**Fix:** Playback chrome (mobile dock, desktop bar, and the extra canvas bottom padding) is hidden from `showTraceLoading()` and only shown in `initTraceUI`, when stream geometry exists. A failed or cleared load hides it again.

**Evidence:**

- Mobile live loading (390px): `web/usability/mobile-playback-dock-loading-live.png`
- Before/after markup fixture (390px): `web/usability/mobile-playback-dock-loading.png`
- Markup fixture used for that screenshot: `web/usability/mobile-playback-dock-loading.html`

### Mobile status sidebar appeared before canvas geometry

**Where:** 390px loaded-trace screen, `#sidebar` Status panel, from Load Trace until the first stream header/geometry message.

**Problem:** `enterMobileTraceLayout()` showed the sidebar immediately and filled it with bootstrap values (`p #0`, `v_i #1`, `step 1`). The canvas HUD was honest, but the panel under it looked like real playback state. Desktop already kept `#sidebar` at `display: none` until `initTraceUI`.

**Fix:** Status chrome stays hidden from `showTraceLoading()` (and a reload clears any previous inline `display: flex`). Bootstrap numbers are not written until stream geometry exists. On mobile the canvas also fills the leftover viewport so hiding the panel does not leave a 52vh stub with empty space below. The sidebar returns in `initTraceUI` with `playback-ready`.

**Evidence:**

- Mobile live loading (390px): `web/usability/mobile-sidebar-loading-live.png`
- Before/after markup fixture (390px): `web/usability/mobile-sidebar-loading.png`
- Markup fixture used for that screenshot: `web/usability/mobile-sidebar-loading.html`

### Mobile playback dock covered Status and blocked scrolling

**Where:** 390px loaded-trace screen after stream geometry arrives (`body.trace-loaded-mobile.playback-ready`), `#sidebar` Status panel and `#mobileTransport`.

**Problem:** The fixed playback dock sat on the Status numbers. `#canvasWrap` used 76px of bottom padding meant for the canvas, which opened a gap between the map and Status instead of clearing the dock. `#app` / `main` stayed viewport-sized with `min-height: 0`, so `scrollHeight` matched the 390x844 viewport and Layers / View / Instructions were clipped with no way to reach them.

**Fix:** The loaded mobile page grows and scrolls. Dock clearance moved to `#sidebar` padding. The empty canvas gap is gone, so Status (p #0, v_i #1, step, alive) sits above the dock at rest. Scrolling to the end keeps Instructions above the dock. Hiding Status (`mobile-panel-closed`) grows the canvas to the leftover viewport above the dock.

**Evidence:**

- Mobile live playback-ready (390px): `web/usability/mobile-dock-sidebar-live.png`
- Mobile live scrolled to Instructions (390px): `web/usability/mobile-dock-sidebar-scrolled-live.png`
- Mobile live with Status hidden (390px): `web/usability/mobile-dock-sidebar-panel-closed-live.png`
- Before/after markup fixture: `web/usability/mobile-dock-sidebar.html`
- Fixture screenshot: `web/usability/mobile-dock-sidebar.png`

### Mobile Instructions listed keyboard shortcuts

**Where:** 390px loaded-trace sidebar, `#sidebar` Instructions card (`body.trace-loaded-mobile.playback-ready`).

**Problem:** The card told phone users to press arrows, Shift, Space, and C / X. Those keys are not available, and playback is the bottom dock (Step, Candidate, Play). Desktop still needs the keyboard table.

**Fix:** At max-width 720px the card describes the dock buttons and touch gestures (Step, Candidate, Play, drag, pinch, Hide / Controls). Desktop (`min-width: 721px`) still shows the keyboard shortcut table.

**Evidence:**

- Mobile live Instructions (390px): `web/usability/mobile-instructions-shortcuts-live.png`
- Before/after markup fixture: `web/usability/mobile-instructions-shortcuts.html`
- Fixture screenshot: `web/usability/mobile-instructions-shortcuts.png`

### Mobile dock had no simplified-segment control

**Where:** 390px loaded-trace screen, `#mobileTransport` after `body.trace-loaded-mobile.playback-ready`.

**Problem:** Desktop can jump previous / next simplified segment with Shift+arrows and the playback-bar prefix buttons. The phone dock only had Step, Candidate, and Play, so that jump had no touch equivalent.

**Fix:** The dock keeps the five existing controls on the first row and adds Segment « / » on a second row (half-width each, so labels stay readable at 390px). Those buttons click the same prefix controls as desktop. Sidebar padding and the closed-panel Controls offset grew with the taller dock. The mobile Instructions table now names Segment.

**Evidence:**

- Mobile live playback-ready dock (390px): `web/usability/mobile-dock-segment-live.png`
- Mobile live Instructions (390px): `web/usability/mobile-dock-segment-instructions-live.png`
- Before/after markup fixture: `web/usability/mobile-dock-segment.html`
- Fixture screenshot: `web/usability/mobile-dock-segment.png`

### Mobile dock had no playback speed control

**Where:** 390px loaded-trace screen, `#mobileTransport` after `body.trace-loaded-mobile.playback-ready`.

**Problem:** Desktop `#playbackBar` has 0.25×-4× speed presets. That bar is `display: none` below 720px, so the phone dock could Play / Pause but auto-advance stayed at 1×.

**Fix:** The dock keeps Step / Candidate / Play and Segment on the first two rows and adds the same five speed presets on a third row. Tapping a rate marks it active on both the dock and the desktop bar. The playback-ready canvas is 38vh so Status still sits above the taller dock. Sidebar padding and the closed-panel Controls offset grew with the extra row. The mobile Instructions table now names Speed.

**Evidence:**

- Mobile live playback-ready dock (390px): `web/usability/mobile-dock-speed-live.png`
- Mobile live Instructions (390px): `web/usability/mobile-dock-speed-instructions-live.png`
- Before/after markup fixture: `web/usability/mobile-dock-speed.html`
- Fixture screenshot: `web/usability/mobile-dock-speed.png`

### Play did not stop at the last frame

**Where:** desktop `#playBtn` and 390px `#mobilePlayBtn`, after a loaded trace reaches the last prefix and last step.

**Problem:** Auto-advance called `goToStep` one past the end. On the last frame that is a no-op, but the play timer kept scheduling ticks and Play stayed on Pause. Tapping 4× only shortened the leftover delay. Clicking Play on the last step also restarted even when later candidates on that step were still unplayed.

**Fix:** Playback asks whether a next candidate, step, or segment exists. If none does, the timer clears and both Play controls return to Play on the last frame. Play / Space only restarts from the first segment when nothing is left to advance.

**Evidence:**

- Desktop live stopped-at-end (1280px): `web/usability/desktop-playback-end-live.png`
- Mobile live stopped-at-end (390px): `web/usability/mobile-playback-end-live.png`
- Before/after markup fixture: `web/usability/playback-end-stop.html`
- Fixture screenshot: `web/usability/playback-end-stop.png`

### Desktop layers formula overflowed the sidebar

**Where:** desktop `#sidebar` Layers list, especially at 900px where the panel is 280px.

**Problem:** Layer labels are a nowrap flex row (checkbox + swatch + MathJax). The single \(S_i[p] = \mathrm{conv}(G_{v_i}) \cap F(S_{i-1}[p], p)\) equation was about 231px after the checkbox, so it painted past the 280px sidebar and 13px past the 900px window. `#sidebar` only set `overflow-y: auto`, so the overflow was visible instead of clipped.

**Fix:** The equation splits at the equals sign so it can wrap. Layer labels wrap, stay `max-width: 100%`, and clip leftover MathJax assistive boxes. The sidebar uses `overflow-x: hidden`. At 900px the visible formula ends at 882px (inside the panel). 390px open layers already fit.

**Evidence:**

- Desktop live 900px: `web/usability/desktop-layers-overflow-live.png`
- Desktop live 1280px: `web/usability/desktop-layers-overflow-1280-live.png`
- Mobile live open layers (390px): `web/usability/mobile-layers-open-live.png`
- Before/after markup fixture: `web/usability/desktop-layers-overflow.html`
- Fixture screenshot: `web/usability/desktop-layers-overflow.png`

### Desktop loading hid Fréchet distance and simplification time

**Where:** desktop header `#paramsBar`, from Load Trace until the first stream header/geometry message (`min-width: 721px`).

**Problem:** `renderParamsBarPreview()` returned immediately on desktop, so the params bar stayed empty while the canvas HUD and `Computing trace…` ran. Computed Fréchet distance and Simplification time only appeared in `initTraceUI`, which grew the 1280px header from 48px to 110px and the 900px header from 81px to 144px. Mobile already reserved those two chips with in-slot spinners.

**Fix:** Loading always writes the two metric chips into `#paramsBar`. Desktop CSS puts a non-empty params bar on its own full-width row so the chips start at the same left edge they keep after values arrive. Mobile 390px behavior is unchanged.

**Evidence:**

- Desktop live loading (1280px): `web/usability/desktop-loading-params-live.png`
- Desktop live loading (900px): `web/usability/desktop-loading-params-900-live.png`
- Mobile live loading regression (390px): `web/usability/mobile-loading-params-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-loading-params.html`
- Fixture screenshot: `web/usability/desktop-loading-params.png`

### Desktop Load Trace button collapsed while computing

**Where:** desktop header `#loadBtn`, from Load Trace until the stream finishes (`min-width: 721px`).

**Problem:** `setLoadButtonBusy(true)` replaced the `Load Trace` label with a 14px spinner. The control shrank from about 82px to 36px, so `Computing trace…` and the rest of the header row jumped left. Mobile already stretches `#loadBtn` to the full form column, so the collapse was a desktop layout shift.

**Fix:** The idle `Load Trace` label stays in the button and keeps its width (`flex-shrink: 0`). The spinner is centered in that reserved slot (`visibility: hidden` on the label, `aria-busy` on the control). The busy button is not dimmed to 0.35 so the in-slot spinner stays readable. Idle and busy both measure 82px at 1280 and 900.

**Evidence:**

- Desktop live loading (1280px): `web/usability/desktop-load-btn-width-live.png`
- Desktop live loading (900px): `web/usability/desktop-load-btn-width-900-live.png`
- Mobile live loading regression (390px): `web/usability/mobile-load-btn-width-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-load-btn-width.html`
- Fixture screenshot: `web/usability/desktop-load-btn-width.png`

### Desktop loading grew the header when the remaining param chips arrived

**Where:** desktop header `#paramsBar`, from Load Trace until stream geometry arrives (`min-width: 721px`).

**Problem:** Loading reserved only Computed Fréchet distance and Simplification time. Those two chips filled one params row (header about 79px at 1280, 113px at 900). When geometry arrived, ε, δ, grid length, disk radius, Fréchet bound, stream counts, and ratio wrapped a second row and grew the header by about 31px (110px at 1280, 144px at 900). The two metric chips stayed put, but the canvas jumped.

**Fix:** Desktop loading writes the same secondary chips as the loaded bar. ε and δ use the form values immediately. The other chips keep ellipsis slots, including Actual Fréchet distance and ratio, so the params row is already two lines. Filling numbers later keeps header height at 110px at 1280. Mobile 390px still shows only the two metric chips.

**Evidence:**

- Desktop live loading (1280px): `web/usability/desktop-loading-params-wrap-live.png`
- Desktop live loaded (1280px): `web/usability/desktop-loading-params-wrap-loaded-live.png`
- Desktop live loading (900px): `web/usability/desktop-loading-params-wrap-900-live.png`
- Desktop live loaded (900px): `web/usability/desktop-loading-params-wrap-900-loaded-live.png`
- Mobile live loading regression (390px): `web/usability/mobile-loading-params-wrap-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-loading-params-wrap.html`
- Fixture screenshot: `web/usability/desktop-loading-params-wrap.png`

### Desktop Computing status grew the first header row

**Where:** desktop header `#uploadStatus` and `.upload-form` ε / δ labels, at 900px, from Load Trace until geometry arrives (`min-width: 721px`).

**Problem:** `Computing trace…` plus its spinner wrapped to two lines (30px). That squeeze also stacked ε / δ under their symbols (35px). The controls row was 41px while loading and 26px after `✓ Loaded Trace 1`, so the header jumped from 158px to 144px. Params stayed 46px; the leftover was the status wrap, not the chips. 1280px already stayed 110px.

**Fix:** `#uploadStatus` stays on one line in a reserved 126×26px slot. Desktop ε / δ labels stay an inline nowrap flex pair, and the upload / filename controls do not wrap. Between 721px and 899px the form wraps as a whole so those reserved widths do not clip. Loading and loaded both measure 144px at 900px and 110px at 1280. Mobile 390px is unchanged at 117px.

**Evidence:**

- Desktop live loading (900px): `web/usability/desktop-loading-status-wrap-900-live.png`
- Desktop live loaded (900px): `web/usability/desktop-loading-status-wrap-900-loaded-live.png`
- Desktop live loading regression (1280px): `web/usability/desktop-loading-status-wrap-1280-live.png`
- Desktop live loaded regression (1280px): `web/usability/desktop-loading-status-wrap-1280-loaded-live.png`
- Mobile live loading regression (390px): `web/usability/mobile-loading-status-wrap-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-loading-status-wrap.html`
- Fixture screenshot: `web/usability/desktop-loading-status-wrap.png`

### Desktop 721px params grew a third row when numbers filled

**Where:** desktop header `#paramsBar`, between 721px and 899px, from Load Trace until stream geometry arrives.

**Problem:** Loading already reserved the secondary chips with ellipsis placeholders. At 721px those chips wrapped to two params rows (header about 184px). When geometry filled `159.0990`-style numbers, chips widened enough to wrap a third row and the header jumped to about 213px. The form row stayed at 67px; only the params bar grew. 900px and 1280px already stayed on two params rows.

**Fix:** Between 721px and 899px, a non-empty `#paramsBar` keeps `min-height: 76px` (three wrap rows). Ellipsis still paints two rows inside that slot; filled numbers use the third row already reserved. Loading and loaded both measure 214px at 721px. 900px stays 150px, 1280px stays 116px, and mobile 390px stays 117px without the reserve.

**Evidence:**

- Desktop live loading (721px): `web/usability/desktop-721-params-wrap-live.png`
- Desktop live loaded (721px): `web/usability/desktop-721-params-wrap-loaded-live.png`
- Desktop live regression (900px): `web/usability/desktop-721-params-wrap-900-regression-live.png`
- Desktop live regression (1280px): `web/usability/desktop-721-params-wrap-1280-regression-live.png`
- Mobile live regression (390px): `web/usability/desktop-721-params-wrap-390-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-721-params-wrap.html`
- Fixture screenshot: `web/usability/desktop-721-params-wrap.png`

### First-visit visitors had no guided intro

**Where:** start screen on desktop and mobile (`#uiTour`), before any trace is loaded.

**Problem:** Novices saw ε / δ, Baseline, and Load Trace with no plain-language path through the first load. Existing Instructions copy assumed the visitor already knew what to do.

**Fix:** A four-step popover tour starts once on first visit (`localStorage` key `simplify-viewer-tour-v1`). It explains the product, choosing a trajectory, what ε / δ mean, and Load Trace, with a spotlight on the live controls. Skip / Done remembers completion. A header `?` control relaunches the tour. ε / δ inputs also expose the same plain-language `title` tooltips.

**Evidence:**

- Desktop live step 2 (1280px): `web/usability/first-visit-tour-1280-live.png`
- Mobile live step 2 (390px): `web/usability/first-visit-tour-390-live.png`
- Before/after markup fixture: `web/usability/first-visit-tour.html`
- Fixture screenshot: `web/usability/first-visit-tour.png`

### Desktop playback bar wrapped Speed onto a second row

**Where:** desktop `#playbackBar`, especially 721-1024px after stream geometry arrives.

**Problem:** Long captions (`pts for current segment`, `segments in simplified curve`, `boundary anchors`) plus full-size inputs made the bar wider than the canvas. Speed presets wrapped under the nav groups, so bar height jumped from about 61px (1280) to about 91px (1024/900). Resizing the window shifted the canvas.

**Fix:** Captions are plain `Step` / `Segment` / `Candidate` (with longer `title` tooltips). Desktop keeps `flex-wrap: nowrap`. Between 721px and 1100px the bar uses compact button/input/speed sizes, and below 820px the captions hide. Height stays one row: 61px at 1280, 53px at 1024/900, 48px at 721, with no horizontal scroll.

**Evidence:**

- Desktop live (1280px): `web/usability/desktop-playback-bar-wrap-1280-live.png`
- Desktop live (1024px): `web/usability/desktop-playback-bar-wrap-1024-live.png`
- Desktop live (900px): `web/usability/desktop-playback-bar-wrap-900-live.png`
- Desktop live (721px): `web/usability/desktop-playback-bar-wrap-721-live.png`
- Mobile regression (390px, desktop bar hidden): `web/usability/desktop-playback-bar-wrap-390-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-playback-bar-wrap.html`
- Fixture screenshot: `web/usability/desktop-playback-bar-wrap.png`

### Layer toggles used paper math without plain language

**Where:** desktop and mobile `#layersSection` Simplify accordion, after a trace is loaded.

**Problem:** Several layer labels were pure notation (`conv(G_{v_i})`, `F(S_{i-1}[p], p)`, the long `S_i[p] = …` equation). Visitors who had not read the paper could not tell what each toggle drew. The long `S_i[p]` equation also left a MathJax assistive MathML tree whose unclipped box could sit about 10px past a 900px window (visible formula and `scrollWidth` already stayed inside).

**Fix:** Every math layer keeps a short symbol plus a plain-language gloss and a longer `title` tooltip. `S_i[p]` is labeled `(current search region)` instead of the full intersection formula. MathJax `renderActions.assistiveMml` is cleared so assistive copies are not injected, with CSS clip kept as a fallback. At 900px the farthest layer label ends at 877px, assistive count is 0, and `documentElement.scrollWidth` stays 900.

**Evidence:**

- Desktop live (900px): `web/usability/desktop-layer-glosses-900-live.png`
- Desktop live (1280px): `web/usability/desktop-layer-glosses-1280-live.png`
- Mobile regression (390px open layers): `web/usability/desktop-layer-glosses-390-regression-live.png`
- Before/after markup fixture: `web/usability/desktop-layer-glosses.html`
- Fixture screenshot: `web/usability/desktop-layer-glosses.png`

### Compare / Baseline selection wrapped the header and Load Trace

**Where:** desktop header `.header-baseline` (and mobile start form), when selecting DOTS / DP / SQUISH before Load Trace.

**Problem:** Compare controls lived on the same form row as Load Trace. Choosing algorithms revealed LSSD / PED ε / Ratio fields that wrapped onto extra rows, moved Load Trace, and at wider desktop widths could make `header { flex-wrap: wrap }` put the title on its own line so the whole header jumped (for example 54px → 142px at 1280 with all three selected). Acronyms also had no plain-language help.

**Fix:** Compare sits on its own full-width strip under Load Trace. Desktop keeps that strip `nowrap` with horizontal scroll if needed, and the header itself is `nowrap` with a shrinking `.file-controls` so the title stays beside the controls. Labels read `Compare:` with short param names (`limit` / `match` / `keep`) plus tooltips; Status uses `start point` / `current point` / `open` instead of bare `alive`. Selecting none → all three keeps header height at 88px (1280), 122px (900), and 123px (721) with Load Trace unmoved.

**Evidence:**

- Desktop live idle (1280px): `web/usability/desktop-compare-bar-wrap-1280-idle-live.png`
- Desktop live all three (1280px): `web/usability/desktop-compare-bar-wrap-1280-live.png`
- Desktop live all three (900px): `web/usability/desktop-compare-bar-wrap-900-live.png`
- Desktop live all three (721px): `web/usability/desktop-compare-bar-wrap-721-live.png`
- Mobile live all three (390px): `web/usability/desktop-compare-bar-wrap-390-live.png`
- Status glosses after load (1280px): `web/usability/desktop-compare-bar-wrap-status-live.png`
- Before/after markup fixture: `web/usability/desktop-compare-bar-wrap.html`
- Fixture screenshot: `web/usability/desktop-compare-bar-wrap.png`

### Cloud Run CI/CD blocked on JSON service-account keys

**Where:** `.github/workflows/deploy.yml` and repository secrets.

**Problem:** The first deploy workflow expected `GCP_SA_KEY`, but the GCP project enforces `constraints/iam.disableServiceAccountKeyCreation`, so a JSON key could not be created or stored. Pushes to `main` could not publish.

**Fix:** Configured Workload Identity Federation for `yeungsinchun/Simplification-of-Trajectory-Streams` (pool `github-actions`, provider `github`, SA `github-actions-deploy@…`) and switched `deploy.yml` to OIDC auth with `id-token: write`. No repository secret is required; optional variables can override project / region / provider / SA.

### After Load Trace, Step / Segment / Candidate were unexplained

**Where:** desktop `#playbackBar` and mobile `#mobileTransport`, after the first successful Load Trace.

**Problem:** The first-visit tour ended at Load Trace. Once the green path appeared, novices faced Step / Segment / Candidate with only short captions and no plain-language walkthrough of how to replay the algorithm. Layers still said Baseline while the header said Compare.

**Fix:** A three-step playback tour opens once after the first successful load (`localStorage` key `simplify-viewer-playback-tour-v1`), spotlighting the playback chrome and explaining Step, then Segment / Candidate / Play. The header `?` stays available after load and relaunches this guide while a trace is loaded; otherwise it relaunches the start tour. Layers accordion and Results empty / hint copy now say Compare instead of Baseline.

**Evidence:**

- Desktop live step 1 (1280px): `web/usability/playback-tour-1280-live.png`
- Desktop live Step (1280px): `web/usability/playback-tour-1280-step-live.png`
- Desktop live Segment/Candidate (1280px): `web/usability/playback-tour-1280-segment-live.png`
- Desktop ? relaunch after load (1280px): `web/usability/playback-tour-1280-relaunch-live.png`
- Mobile live step 1 (390px): `web/usability/playback-tour-390-live.png`
- Mobile ? visible after skip (390px): `web/usability/playback-tour-390-help-live.png`
- Before/after markup fixture: `web/usability/playback-tour.html`
- Fixture screenshot: `web/usability/playback-tour.png`

### Results tab and Dead candidates stayed opaque after load

**Where:** left-edge `#resultsPanelOpen`, Layers `#toggle-dead-candidates`, and the post-load playback tour.

**Problem:** After Load Trace the Results tab appeared only when the compare API finished, and the playback tour never mentioned it. Novices who finished Step / Segment / Candidate guidance still did not know Results holds scores or optional Compare runs. Layers also labeled ruled-out points as Dead candidates with no gloss.

**Fix:** Results chrome is shown as soon as the trace is ready (`showResultsPanel` before the playback tour; `clearCompare({ hideChrome: false })` while compare data reloads). Playback tour step 5/5 spotlights Results and explains optional Compare. On mobile the Results tab sits above the playback dock instead of under it. The layer toggle reads Rejected candidates with a plain-language tooltip; the Results button title / aria-label name scores and Compare.

**Evidence:**

- Before/after markup fixture: `web/usability/results-tour.html`
- Fixture screenshot: `web/usability/results-tour.png`
- Desktop live Results step (1280px): `web/usability/results-tour-1280-live.png`
- Mobile live Results step (390px): `web/usability/results-tour-390-live.png`

### Layers stayed paper-jargon and were skipped by the post-load tour

**Where:** sidebar `#layersSection` / `#fitBtn`, and the post-load playback tour.

**Problem:** After Load Trace, layer rows still said reachability / δ-ball / vertex / boundary anchors, and the playback tour jumped from Segment / Candidate to Results. Novices could finish both tours without learning that Layers control map overlays or that Fit to data resets the view.

**Fix:** Layer glosses and tooltips use plain wording (search circle, allowed area, anchor points). Playback tour inserts a Layers step (4/5) that opens the mobile Layers accordion, spotlights the Simplify toggles, and mentions Fit to data. Results / Compare remains the final step. Results empty copy no longer implies Compare is required before scores appear.

**Evidence:**

- Before/after markup fixture: `web/usability/layers-tour.html`
- Fixture screenshot: `web/usability/layers-tour.png`
- Desktop live Layers step (1280px): `web/usability/layers-tour-1280-live.png`
- Mobile live Layers step (390px): `web/usability/layers-tour-390-live.png`

### Params bar and Results still used Fréchet / |stream| jargon

**Where:** desktop `#paramsBar` chips after Load Trace, mobile Match error / Time metrics, and Results compare table.

**Problem:** After the tours, novices still faced `Computed Fréchet distance`, `len_grid`, `R (disk radius)`, `a-priori Fréchet bound`, `|stream|`, `|simplified|`, `ratio`, and a Results row labeled `Frechet`. Those labels assume paper vocabulary and MathJax typesetting.

**Fix:** Params chips use short plain labels with tooltips (`Match error`, `Time`, `grid step`, `radius`, `error budget`, `trace error`, `original`, `kept`, `kept %`). Results uses `Match error` / `Time (ms)` and a plain footer gloss. Desktop chips no longer need MathJax. Desktop `header` may wrap again so `#paramsBar` keeps `flex: 1 0 100%` on its own row (file-controls still shrinks beside the title); that stops the bar from being crushed into a ~160px side column at 721px after the Compare nowrap change.

**Evidence:**

- Before/after markup fixture: `web/usability/desktop-params-glosses.html`
- Fixture screenshot: `web/usability/desktop-params-glosses.png`
- Desktop live (1280px): `web/usability/desktop-params-glosses-1280-live.png`
- Desktop live (900px): `web/usability/desktop-params-glosses-900-live.png`
- Mobile live metrics (390px): `web/usability/desktop-params-glosses-390-live.png`

### ε/δ and Status still used opaque shorthand for novices

**Where:** start-screen `#epsilonInput` / `#deltaInput`, drop hint / mobile start help, Status `#statusGrid`, Layers `#toggle-stream`.

**Problem:** After Match error / tour work, novices who skipped or forgot the tour still saw bare `ε` / `δ` with no visible meaning, Status rows labeled `step` / `open`, Layers `Full stream`, and start copy that said “set ε/δ” without explaining the words.

**Fix:** Inputs show short glosses (`ε match`, `δ grid`) with aria-labels; below 900px desktop the gloss text hides so Load Trace stays beside the fields (tooltips + start copy still explain). Drop hint and mobile instructions say `ε match` / `δ grid`. Status uses `path step` / `candidates`. Layers says `Original path`. Invalid-input alert uses the same plain wording.

**Evidence:**

- Before/after markup fixture: `web/usability/eps-delta-glosses.html`
- Fixture screenshot: `web/usability/eps-delta-glosses.png`
- Desktop live (1280px): `web/usability/eps-delta-glosses-1280-live.png`
- Desktop live (900px): `web/usability/eps-delta-glosses-900-live.png`
- Desktop live (721px): `web/usability/eps-delta-glosses-721-live.png`
- Desktop live (820px, glosses hidden): `web/usability/eps-delta-glosses-820-live.png`
- Mobile live (390px): `web/usability/eps-delta-glosses-390-live.png`
- Desktop loaded Status/Layers (1280px): `web/usability/eps-delta-glosses-1280-loaded-live.png`
- Mobile loaded Status/Layers (390px): `web/usability/eps-delta-glosses-390-loaded-live.png`

### Post-load Compare had no visible Run next to the header pills

**Where:** desktop `#headerBaseline` Compare strip after Load Trace; `#baselineStatus` / `#baselineRunBtn` only inside the closed Results panel.

**Problem:** After a trace loaded, selecting DOTS / DP / SQUISH in the header only wrote “Click Run compare…” into the Results panel status node. With Results closed, novices saw param fields appear but no next action, and tooltips still said “baseline”.

**Fix:** A compact header `Run` button (`#headerBaselineRunBtn`) appears once Results is available, stays on the same nowrap Compare strip, and shares busy/disabled state with Results `Run compare`. Start and playback tour copy mention Compare + Run; pill titles drop “baseline” jargon; failure copy says “Compare run failed”.

**Evidence:**

- Before/after fixture: `web/usability/desktop-compare-header-run.html`
- Fixture screenshot (900px): `web/usability/desktop-compare-header-run.png`
- Desktop live after load with DOTS + Run (1280px): `web/usability/desktop-compare-header-run-1280-live.png`
- Desktop live narrow (900px): `web/usability/desktop-compare-header-run-900-live.png`

### Compare params and Layers hint still used thresh / Results-only copy

**Where:** `#headerBaseline` param labels, Results Compare fields, `#baselineLayerHint`, and mobile start Instructions.

**Problem:** After the header `Run` control landed, Layers Compare still said “Run a compare from the Results panel…”, so novices who used the header path got contradictory guidance. Header fields also kept paper-ish `DOTS thresh` / `DP ε` while other chrome already used plain match/limit wording. Mobile start help mentioned Compare before Load Trace but never mentioned Run or Results.

**Fix:** Header/Results fields read `DOTS limit`, `DP match` / `match error`, and `SQUISH keep`. Layers hint and error strings point at header `Run` or Results `Run compare`. Mobile Instructions put Compare + Run after Load Trace and mention Results for scores. The 900px fixture strip stays one row with no page overflow.

**Evidence:**

- Before/after fixture: `web/usability/compare-param-glosses.html`
- Fixture screenshot (900px): `web/usability/compare-param-glosses.png`

### Upload button and mobile Compare path used opaque / unreachable wording

**Where:** `#uploadBtn`, `#dropHint` start copy, mobile `.mobile-start-help`, and the post-load Results tour step.

**Problem:** The primary upload control said `Upload original.txt`, which assumes repository sample filenames. Mobile Instructions told novices to pick Compare and tap `Run` after load, but `body.trace-loaded-mobile` hides `.file-controls` (including the header Compare strip and `Run`). The Results tour also led with “header Run” even though phones only keep `Run compare` inside Results.

**Fix:** The button reads `Upload trajectory` with a plain-text format tooltip; drop hint and start-tour copy match. Mobile Instructions send users to Results → Run compare after Load Trace. The Results tour puts Results / Run compare first and mentions header Run only for wider screens. The 900px fixture strip stays one row with no page overflow.

**Evidence:**

- Before/after fixture: `web/usability/upload-trajectory-label.html`
- Fixture screenshot (900px): `web/usability/upload-trajectory-label.png`
- Desktop fixture (1280px): `web/usability/upload-trajectory-label-1280.png`
- Mobile fixture (390px): `web/usability/upload-trajectory-label-390.png`
- Desktop live start (1280px): `web/usability/upload-trajectory-label-1280-live.png`
- Mobile live start (390px): `web/usability/upload-trajectory-label-390-live.png`

### Layers Compare hint still led with unreachable header Run on phones

**Where:** `#baselineLayerHint` in Layers → Compare, and `#baselineStatus` after selecting DOTS / DP / SQUISH.

**Problem:** Iteration 17 fixed mobile Instructions and the Results tour, but Layers still said “press Run beside Compare…”, and selection status in Results used the same desktop-first wording. After Load Trace, `body.trace-loaded-mobile` hides `.file-controls`, so header Run is unreachable on phones.

**Fix:** Layers hint and selection status lead with Results → Run compare, and mention header Run only as a wider-screen alternative. Empty-selection status says “press Run compare”. Fixture and live 390px evidence confirm file-controls are hidden while the new hint stays on-screen with no page overflow (`scrollWidth` 390).

**Evidence:**

- Before/after fixture: `web/usability/layers-compare-hint-mobile.html`
- Fixture screenshot (900px): `web/usability/layers-compare-hint-mobile.png`
- Fixture screenshot (390px): `web/usability/layers-compare-hint-mobile-390.png`
- Fixture screenshot (1280px): `web/usability/layers-compare-hint-mobile-1280.png`
- Mobile live loaded Layers (390px): `web/usability/layers-compare-hint-mobile-390-live.png`
- Desktop live loaded Layers (1280px): `web/usability/layers-compare-hint-mobile-1280-live.png`

## Still open

No open layout, novice-copy, or Cloud Run CI/CD items from this pass.

### Match / Match error still said the original (resolved)

**Where:** Match form/chip titles, Match error Scores footer and tips, match limit / kept % tips, Gray path overlay tip, Accuracy tour, and playback-tour intro.

**Problem:** Map overlays and Step already used `Gray path`, but Match / Match error / tour tips still said `the original` / `gray original`. Skip-tour users mapping Match error to the Gray path overlay could not tell those tips meant that same path.

**Fix:** Tips say `Gray path` / `original points` (count metric kept). Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/match-error-gray-path.html`
- Fixture screenshot (900px): `web/usability/match-error-gray-path.png`

### Step still said original points (resolved)

**Where:** Status Step tooltip, desktop/mobile Step titles and aria-labels, Instructions gloss, playback tour Step copy, and Status point-index title.

**Problem:** Map overlays already named the input overlay `Gray path`, but Step still said `original points` / `original path points`, and Status indices said `original trajectory`. Skip-tour users mapping Step and start/current indices to Gray path could not tell those controls meant that same path.

**Fix:** Copy says `Gray path points` / `Gray path`. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/step-gray-path-points.html`
- Fixture screenshot (900px): `web/usability/step-gray-path-points.png`

### Match / Grid / Scores still said simplified path (resolved)

**Where:** Match / Grid form titles, Load title/aria-label, Results This run / Match error / Kept points / Time tips, Map overlays This run / Gray path tips, desktop `#paramsBar` chips, Accuracy tour copy.

**Problem:** Green path so far / Full green path / next green-path point already shared green-path vocabulary, but Match / Grid / Load / Scores / Gray path tips still said `simplified path` or `simplification`. Skip-tour users mapping those controls to the green path could not tell the tips meant that same path.

**Fix:** Tips say `green path` / `this run` / `each path` / `build the green path`. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/green-path-not-simplified.html`
- Fixture screenshot (900px): `web/usability/green-path-not-simplified.png`

### Map overlays / circle radius said next simplified point (resolved)

**Where:** Map overlays `#toggle-F-Si` / `#toggle-S` tooltips and the desktop `#paramsBar` circle radius chip tooltip.

**Problem:** Green path so far / Full green path / Next-point zone already shared green-path and next-point vocabulary, but those tips still said `next simplified point`. Skip-tour users mapping Next-point zone and circle radius to the green path could not tell the tip meant the next point on that green path.

**Fix:** Tips read `next green-path point`. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/next-green-path-point.html`
- Fixture screenshot (900px): `web/usability/next-green-path-point.png`

### Map overlays said Still-allowed area (resolved)

**Where:** Map overlays `#toggle-F` / `#toggle-F-Si` labels/tooltips and `#toggle-S` Next-point zone tooltip.

**Problem:** Match / Match error / Match limit already shared Match vocabulary, but the blue/cyan overlays still said `Still-allowed area` / `Still-allowed (this option)`. Skip-tour users could not tell those areas are the residual region within the Match limit that Next-point zone is cut from.

**Fix:** Labels read `Match-safe area` and `Match-safe (this option)`. Related tooltips and the Next-point zone tip reuse those names. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/match-safe-area.html`
- Fixture screenshot (900px): `web/usability/match-safe-area.png`

### Map overlays said Path so far (resolved)

**Where:** Map overlays `#toggle-simplified` label/tooltip, Gray path / Full green path tooltips, and playback Map overlays tour copy.

**Problem:** Gray path and Full green path already led with color, but the mid-playback green overlay still said `Path so far`. Skip-tour users mapping Map overlays to the gray/green vocabulary could not tell that toggle was the same green path.

**Fix:** Label reads `Green path so far`. Related tooltips and the Map overlays tour use the same name. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/green-path-so-far.html`
- Fixture screenshot (900px): `web/usability/green-path-so-far.png`

### Map overlays said Original path (resolved)

**Where:** Map overlays `#toggle-stream` label/tooltip and playback Map overlays tour copy.

**Problem:** Path so far / Full green path and the playback tour already used green / gray wording, but the input overlay toggle still said `Original path`. Skip-tour users mapping Map overlays to the tour’s “gray original” could not tell which toggle was that gray path.

**Fix:** Label reads `Gray path` with a tooltip that names Path so far and Full green path. Tour copy says gray path. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/gray-path.html`
- Fixture screenshot (900px): `web/usability/gray-path.png`

### Map overlays said Search circle (start/current point) (resolved)

**Where:** Map overlays `#toggle-ball-p0` / `#toggle-ball-pi` labels/tooltips; desktop `#paramsBar` `circle radius` chip tooltip; playback Map overlays tour copy.

**Problem:** Status already led with `start point` / `current point`, but the purple/pink toggles still said `Search circle (start point)` / `Search circle (current point)`. Skip-tour users opening Map overlays first saw algorithm-first naming instead of Status vocabulary, and the circle radius chip tip still said Search circle.

**Fix:** Labels read `Start-point circle` and `Current-point circle`. Tooltips and the circle radius chip tip reuse those names; the Map overlays tour says start-point / current-point circles. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/start-point-circle.html`
- Fixture screenshot (900px): `web/usability/start-point-circle.png`

### Map overlays said Full result path (resolved)

**Where:** Map overlays `#toggle-final-simplify` label/tooltip.

**Problem:** Segment, Path so far, and This run already shared green-path wording, but the finished-path toggle still said `Full result path`. Skip-tour users could not tell it was the same green path shown mid-playback by Path so far / Segment.

**Fix:** Label reads `Full green path` with a tooltip that names Path so far and Segment. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/full-green-path.html`
- Fixture screenshot (900px): `web/usability/full-green-path.png`

### Map overlays said Next landing zone (resolved)

**Where:** Map overlays `#toggle-S` label/tooltip, Still-allowed tooltip, and playback Map overlays tour copy.

**Problem:** After Option / Options near current / Option markers shared next-point vocabulary, the purple overlap toggle still said `Next landing zone`. Skip-tour users could not tell that zone is the same next-point idea as Option on the playback bar.

**Fix:** Label reads `Next-point zone` with a tooltip that names Options near current, Still-allowed area, and Option on the playback bar. Still-allowed and tour copy use the same name. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/next-point-zone.html`
- Fixture screenshot (900px): `web/usability/next-point-zone.png`

### Map overlays said Paths & search (resolved)

**Where:** Map overlays accordion `#accordionSimplify` summary and playback Map overlays tour copy.

**Problem:** Results Scores already labeled the green-path column `This run`, but the Map overlays accordion for the same run still said `Paths & search`. Skip-tour users mapping Scores to overlay toggles saw two names for one run.

**Fix:** Accordion summary reads `This run` with a Scores-aligned tooltip. Tour copy uses the same name and mentions the Scores column. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/this-run-overlays-accordion.html`
- Fixture screenshot (900px): `web/usability/this-run-overlays-accordion.png`

### Playback and Status said Candidate (resolved)

**Where:** Status Option row, desktop `#playbackBar` Option caption, mobile Option dock buttons, Instructions / tour / Map overlays tooltips that named Candidate.

**Problem:** After Option markers / Options near current landed, Status, playback, and Instructions still said `Candidate`. Skip-tour users could not tell that the playback control cycles the same next-point options as those overlays.

**Fix:** Visible labels read `Option` (narrow desktop short `Opt`). Status, Instructions, tour, and overlay tooltips use the same Option vocabulary. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/option-playback-label.html`
- Fixture screenshot (900px): `web/usability/option-playback-label.png`

### Map overlays said Candidate markers (resolved)

**Where:** Map overlays Paths & search toggles `#toggle-P` / `#toggle-dead-candidates` / `#toggle-F-Si`.

**Problem:** After Options near current landed, orange/red marker toggles still said `Candidate markers` / `Rejected candidates`, and the cyan toggle said `Still-allowed (this candidate)`. Skip-tour users could not map those markers to Options near current, or tell that Candidate on the playback bar cycles the same options.

**Fix:** Labels read `Option markers`, `Rejected options`, and `Still-allowed (this option)`. Tooltips link Options near current and Candidate playback. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/option-markers.html`
- Fixture screenshot (900px): `web/usability/option-markers.png`

### Map overlays said Candidate region (resolved)

**Where:** Map overlays Paths & search toggle `#toggle-Gi`, `#toggle-S` tooltip, and playback Map overlays tour copy.

**Problem:** After Next landing zone and Still-allowed area landed, the yellow overlay still said `Candidate region`. Skip-tour users could not tell it is the set of next-point options near the Status current point, and the tour still said “candidate regions”.

**Fix:** Label reads `Options near current` with a Status-linked tooltip. Next landing zone tooltip and the Map overlays tour use the same capitalized name. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/options-near-current.html`
- Fixture screenshot (900px): `web/usability/options-near-current.png`

### Map overlays said Allowed area so far (resolved)

**Where:** Map overlays Paths & search toggles `#toggle-F` / `#toggle-F-Si` and `#toggle-S` tooltip.

**Problem:** After Path so far and Next landing zone landed, the blue/cyan overlays still said `Allowed area so far` / `Allowed area (this candidate)`. Skip-tour users could confuse “so far” with Path so far, and the Next landing zone tooltip already said “still-allowed” without a matching toggle name.

**Fix:** Labels read `Still-allowed area` and `Still-allowed (this candidate)`. Tooltips lead with Still-allowed / Match / Candidate / Next landing zone vocabulary. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/still-allowed-area.html`
- Fixture screenshot (900px): `web/usability/still-allowed-area.png`

### Map overlays said Current search region (resolved)

**Where:** Map overlays Paths & search toggle `#toggle-S`.

**Problem:** After Candidate region / Allowed area / Search circle labels landed, the purple overlap toggle still said `Current search region`. Skip-tour users could not tell that the zone is where the next simplified point can land.

**Fix:** Label reads `Next landing zone`. Tooltip leads with that purpose and mentions the candidate / still-allowed overlap. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/next-landing-zone.html`
- Fixture screenshot (900px): `web/usability/next-landing-zone.png`

### Params chip said search radius instead of Search circle (resolved)

**Where:** desktop `#paramsBar` `search radius` chip; Map overlays Search circle (start / current) tooltips.

**Problem:** After Match / Grid / grid cell / saved Match landed, the radius chip still said `search radius` while Map overlays said `Search circle`. Circle tooltips also claimed the size was “from Grid”, so skip-tour users could not map the chip number to the purple/pink overlays.

**Fix:** Chip reads `circle radius` with a Search circle tooltip. Overlay titles say the size matches the circle radius chip. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/circle-radius-search-circle.html`
- Fixture screenshot (900px): `web/usability/circle-radius-search-circle.png`

### Params chips still said cell size / file match (resolved)

**Where:** desktop `#paramsBar` secondary chips after Load.

**Problem:** After Match / Grid / match limit landed, chips still said `cell size` and `file match`. Skip-tour users could not tell that the first is a Grid-derived length, or that the second is a saved Match error (distinct from the live Match error chip).

**Fix:** Chips read `grid cell` and `saved Match` with Grid- and Match-aligned tooltips. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/saved-match-grid-cell.html`
- Fixture screenshot (900px): `web/usability/saved-match-grid-cell.png`

### Params chips still said max error / orig. points (resolved)

**Where:** desktop `#paramsBar` secondary chips; SQUISH Compare field titles (header + Results).

**Problem:** After Match error / kept points landed, chips still said `max error` beside Match error and abbreviated `orig. points` beside `kept points`. SQUISH titles still said `keep ratio`. Skip-tour users could not tell the Match limit from the live Match error, or map the count chips to one vocabulary.

**Fix:** Chips read `match limit` and `original points` with Match-aligned tooltips. SQUISH titles say `keep percent`. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/match-limit-original-points.html`
- Fixture screenshot (900px): `web/usability/match-limit-original-points.png`

### Scores Metric and DP/DOTS tooltips stayed paper-jargon (resolved)

**Where:** Results Scores table first column; Compare pill / overlay titles and DP match field titles (header + Results).

**Problem:** The Scores table headed the first column `Metric`, and Compare tooltips still said `Douglas-Peucker`, `point-to-edge`, and `streaming simplifier`. Skip-tour users mapping Scores to Compare could not tell what those algorithms do from long-press titles alone.

**Fix:** Scores uses `Score` to match the section name. DOTS / DP titles say stream / classic path shortener with a plain match-limit gloss; DP fields say match error. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/compare-score-plain-tooltips.html`
- Fixture screenshot (900px): `web/usability/compare-score-plain-tooltips.png`

### Step / Segment still said “simplified piece” (resolved)

**Where:** Instructions gloss, Status Segment/Step tooltips, desktop playback captions/titles, mobile Step/Segment controls, Map overlays search-circle / Allowed area tooltips, and the playback tour Step / Segment steps.

**Problem:** After Status and playback already shared Segment / Step / Candidate labels, the glosses still said “simplified piece” / “pieces of the simplified path”. Skip-tour users could not tell that Step walks inside the current Segment, or that Segment pieces are the green path.

**Fix:** Step reads “within the current Segment”; Segment / tour / tooltips say “pieces of the green path” (and Map overlays use current/previous Segment). Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/step-segment-green-path.html`
- Fixture screenshot (900px): `web/usability/step-segment-green-path.png`

### Wide desktop form still led with Greek ε / δ (resolved)

**Where:** start-screen `#epsilonInput` / `#deltaInput` labels at ≥900px (and the shared base form CSS used on every viewport).

**Problem:** Narrow desktop and phones already showed plain `Match` / `Grid`, and help / tour / chips used the same words, but wide desktop still painted Greek-first `ε match` / `δ grid`. Skip-tour users on a typical laptop saw a different vocabulary than Instructions.

**Fix:** Base CSS now hides `.param-symbol` and capitalizes Match / Grid on every viewport (ε / δ remain in titles). Narrow-desktop still shrinks the number inputs so Load stays one row. Fixture `scrollWidth` stays within 1280.

**Evidence:**

- Before/after fixture: `web/usability/wide-match-grid-labels.html`
- Fixture screenshot (1280px): `web/usability/wide-match-grid-labels.png`

### Start help / tour / alerts still led with Greek ε / δ (resolved)

**Where:** mobile and desktop start Instructions, empty-canvas drop hint, Accuracy tour step, invalid Match/Grid alert, server timeout / invalid-εδ JSON, DP Compare tooltips, search-circle layer tooltips; mobile form labels at max-width 720px.

**Problem:** After loaded params chips and 721–899px fields preferred plain `Match` / `Grid`, start help, the first-visit Accuracy tour, validation alerts, and server errors still led with Greek `ε match` / `δ grid`. Phones also kept Greek symbols above the stacked number fields, so skip-tour users saw a different vocabulary than the chips.

**Fix:** Instructions, tour, alerts, and server copy now lead with `Match` / `Grid`. Form labels on every viewport (including phones and wide desktop) hide Greek symbols and capitalize the glosses. ε / δ remain only in tooltips. Fixture `scrollWidth` stays within 390.

**Evidence:**

- Before/after fixture: `web/usability/match-grid-help-copy.html`
- Fixture screenshot (390px): `web/usability/match-grid-help-copy.png`

### Narrow desktop hid Match / Grid glosses (resolved)

**Where:** start-screen `#epsilonInput` / `#deltaInput` labels at 721–899px; loaded `#paramsBar` chips for ε / δ.

**Problem:** To keep Load on one row, the 721–899px rule hid `.param-gloss`, leaving bare Greek `ε` / `δ`. Skip-tour users (and anyone who forgot the tour) could not tell what the fields meant without hovering tooltips. Loaded params chips also led with `ε match` / `δ grid`.

**Fix:** At 721–899px the form keeps capitalized Match / Grid glosses (now the default on every viewport) and slightly narrower number inputs so Load stays on the same row (`scrollWidth` 820 / 900). Loaded params chips read `Match` / `Grid` with ε / δ only in tooltips.

**Evidence:**

- Before/after fixture: `web/usability/narrow-match-grid-labels.html`
- Fixture screenshot (820px): `web/usability/narrow-match-grid-labels.png`

### Server / Match error failures still said binary jargon (resolved)

**Where:** `web/server.py` API error payloads (`type:error` stream messages, Compare `baseline_error`, Frechet / upload / missing-trajectory JSON); Match error chip in `viewer.js` `renderParamsBar`.

**Problem:** After client load/stream copy was plain, server-authored failures still said `Binary execution failed`, `Trace N not found`, `Simplify failed`, `lssd must be positive`, and raw binary paths / stderr. The Match error chip also showed bare `failed`. Novices who hit a timeout, missing Compare binary, or Match error outage saw developer wording.

**Fix:** Server responses now use plain trajectory / Compare / Match error guidance (technical detail stays in server logs). The client sanitizes residual technical `error` / `baseline_error` strings, and the Match error chip shows `unavailable` with a short tooltip. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/server-error-plain-copy.html`
- Fixture screenshot (900px): `web/usability/server-error-plain-copy.png`

### Load/stream failures still said JSON / prefix jargon (resolved)

**Where:** `loadTraceStream` / `loadTraceText` failure paths in `viewer.js` (status text and alerts when a response or fallback file cannot be read).

**Problem:** Rare load failures still said `Invalid stream JSON`, `Received prefix before header`, `Received done before header`, `Could not parse JSON` plus raw parser text, and `Parsing NKB JSON…`. Novices who hit a bad server response or unreadable fallback file saw developer wording instead of what to try next.

**Fix:** Those paths now say plain trajectory copy (`Could not read trajectory data…`, `Trajectory data arrived out of order…`, `Reading trajectory…`, and a plain-text upload/pick alert). Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/load-error-plain-copy.html`
- Fixture screenshot (900px): `web/usability/load-error-plain-copy.png`

### Narrow desktop hid Step/Segment/Candidate captions (resolved)

**Where:** desktop `#playbackBar` `.pb-caption` at 721-820px.

**Problem:** To keep the playback bar one row, full captions were `display: none` below 820px, so Step / Segment / Candidate became three identical `← 0/0 →` groups. Skip-tour users (and anyone who forgot the tour) could not tell the controls apart without hovering tooltips.

**Fix:** At 721-820px the bar keeps short `Step` / `Seg` / `Cand` labels instead of hiding captions. Each nav group and input also has an aria-label matching the Instructions glosses. Fixture stays one row with `scrollWidth` within 820.

**Evidence:**

- Before/after fixture: `web/usability/playback-short-captions.html`
- Fixture screenshot (820px): `web/usability/playback-short-captions.png`

### Speed chips lost meaning when the label hid; upload format showed quotes (resolved)

**Where:** desktop `#playbackBar` `.speed-preset` buttons (especially 721-1100px where `.pb-speed-label` is `display: none`); empty-canvas `#dropHint .dropHint-format`.

**Problem:** Narrow desktop hides the Speed caption so the playback bar stays one row, but the `0.25×`–`4×` chips had no `title` / `aria-label`, so skip-tour users saw bare multipliers. The upload format line also rendered `"x y"` with quotes, which novices could copy into trajectory files even though samples are bare `x y`.

**Fix:** Each speed chip (desktop and mobile) has a Playback speed title and aria-label; the speed group is labeled too. The format line reads `N (point count)` then `x y` without quotes. Fixture Speed row stays one line with `scrollWidth` within 900.

**Evidence:**

- Before/after fixture: `web/usability/speed-titles-upload-format.html`
- Fixture screenshot (900px): `web/usability/speed-titles-upload-format.png`

### Status lacked a map/playback gloss; empty Load said upload (resolved)

**Where:** sidebar `#statusGloss`; empty `#loadBtn` status in `viewer.js`; default `#baselineLayerHint` (preloaded path).

**Problem:** Skip-tour users opening Status first saw bare start point / current point / Segment numbers with no link to the colored map markers or the playback bar. Pressing Load with nothing chosen said “select or upload”, which is unreachable on phones where Upload is `display: none`. The preloaded Compare map-overlay hint still ended with an Upload aside even when Compare was already available.

**Fix:** Status leads with a plain gloss mapping start/current to map markers and Segment/Step/Candidate to playback. Empty Load says “Please choose a trajectory first”. The default Compare hint keeps Results → Run compare and drops the Upload aside (upload-blocked copy still explains uploads when needed). Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/status-gloss-map-playback.html`
- Fixture screenshot (900px): `web/usability/status-gloss-map-playback.png`

### Mobile start help still mentioned Upload (resolved)

**Where:** header `.mobile-start-help` Instructions (max-width 720px); start-tour Load step and playback Results tour copy in `viewer.js`.

**Problem:** On phones, `#uploadBtn` / `.or-divider` are `display: none`, so Upload is unreachable. Start Instructions still said “Uploaded files show scores only” and over-qualified Compare as if Upload were an option. Tour Load / Results steps repeated that Upload dead-end language for every visitor.

**Fix:** Mobile start help points at Compare above, then Results → Run compare, with no Upload wording. Start-tour Load and Results tour keep the preloaded Compare path but drop the “Uploaded files show scores only” aside (desktop Upload remains explained only in the Choose trajectory step). Fixture `scrollWidth` stays within 390.

**Evidence:**

- Before/after fixture: `web/usability/mobile-start-no-upload-copy.html`
- Fixture screenshot (390px): `web/usability/mobile-start-no-upload-copy.png`

### Compare overlays / Scores headers lacked glosses (resolved)

**Where:** Map overlays → Compare `#baselineLayerToggles` after Run compare; Results Scores `#compareMetricsHead` algorithm columns.

**Problem:** Results Compare pills already showed `stream` / `classic` / `keep %`, but after a Compare run the Map overlays toggles and Scores table headers were bare `DOTS` / `DP` / `SQUISH`. Skip-tour users mapping dashed map paths or score columns to the pills lost the plain gloss.

**Fix:** Overlay toggles show the same short glosses with map-path tooltips. Scores headers reuse those glosses under the acronym (accordion summary stays short `DOTS / DP / SQUISH`). Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/compare-overlay-glosses.html`
- Fixture screenshot (900px): `web/usability/compare-overlay-glosses.png`

### Outer Layers section still used GIS jargon (resolved)

**Where:** sidebar `#layersSection` h2 / `#mobileLayersToggle`, playback tour step, mobile Instructions / start help after Load.

**Problem:** Iteration 32 renamed the inner Simplify accordion to `Paths & search`, but the outer section heading and mobile toggle still said `Layers`. Tour, Instructions, and start help used the same GIS word, so skip-tour users saw two names for one sidebar block that also holds Compare.

**Fix:** Outer section and mobile toggle read `Map overlays`. Tour, Hide/Controls shortcut, and mobile start help use the same wording and mention Paths & search / Compare. Upload Compare hint drops the leftover “layers” suffix. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/map-overlays-label.html`
- Fixture screenshot (900px): `web/usability/map-overlays-label.png`

### Layers still said segment start / Anchor points (resolved)

**Where:** Layers `#toggle-ball-p0` / `#toggle-P` / `#toggle-F-Si` and related layer tooltips after Load.

**Problem:** Status and the map already used `start point` / `current point`, but Layers still said `Search circle (segment start)` and `Anchor points`, and tooltips led with paper symbols (`Paper: Bp.`, `Paper: P.`). Skip-tour users opening Layers first could not map toggles to Status, and long-press titles on phones dumped math notation.

**Fix:** Layers read `Search circle (start point)`, `Allowed area (this candidate)`, and `Candidate markers` (paired with Rejected candidates). Tooltips use the same start point / current point / δ grid wording and drop `Paper:` math. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/layers-start-point-labels.html`
- Fixture screenshot (900px): `web/usability/layers-start-point-labels.png`

### Status Candidates did not match playback Candidate (resolved)

**Where:** sidebar Status `#statusGrid` Candidates row after Load.

**Problem:** Segment and Step already used `N / total` like the playback bar, but Status still showed `Candidates 3 still open`. Skip-tour users mapping Status to Candidate on the dock could not tell which option was selected or how many options existed.

**Fix:** Status labels the row `Candidate` and shows the same `current / total` as the playback Candidate control (still-open + just-rejected cycle pool). The still-open count moves into the tooltip. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/status-candidate-total.html`
- Fixture screenshot (900px): `web/usability/status-candidate-total.png`

### Status omitted Segment and used # point indices (resolved)

**Where:** sidebar Status `#statusGrid` / `#statusIndices` after Load.

**Problem:** Status showed Step as a bare number with no total and never listed Segment, while the playback bar already used `N / total` for both. Point indices also used a cryptic `#` prefix (`#0`, `#12`). Skip-tour users opening Status first could not map it to Segment / Step on the dock.

**Fix:** Status lists Segment and Step as `current / total` with the same tooltips as the playback captions. Point indices are plain numbers with a “Point number on the original trajectory” tooltip. Mobile Status reserves height for the third row so the panel does not jump.

**Evidence:**

- Before/after fixture: `web/usability/status-segment-step-totals.html`
- Fixture screenshot (900px): `web/usability/status-segment-step-totals.png`

### Results still said Simplified points (resolved)

**Where:** Results Scores first metric row after Load; Kept % tooltip; header `kept %` chip tooltip.

**Problem:** Header chips already said `kept points` / `kept %`, but Results still labeled the count row `Simplified points` and Kept % tooltips still led with that phrase. Novices mapping Scores to the header saw two names for the same keep count.

**Fix:** Scores row is `Kept points` with a tooltip that points at the header chip. Kept % tooltips (Scores + header chip) say “Kept points as a percent…”. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/scores-kept-points.html`
- Fixture screenshot (900px): `web/usability/scores-kept-points.png`

### DOTS K unit and opaque Loaded id (resolved)

**Where:** header / Results DOTS Compare field (`limit` + unit `K`); preloaded load success `#uploadStatus`; desktop `#paramsBar` count chips `original` / `kept`.

**Problem:** Preloaded success still said `✓ Loaded 51` (internal id) while Upload already said `✓ Loaded`. DOTS used `limit` with unit `K` titled only “Thousand”, so novices could not tell the field was a distance budget entered in thousands. Params chips `original` / `kept` looked like bare adjectives next to numbers.

**Fix:** Preloaded success matches Upload as `✓ Loaded` (header title already names the trajectory). DOTS reads `budget` / `distance budget` with unit `×1k` and tooltips that spell out thousands. Count chips are `orig. points` / `kept points`. Fixture Compare strip stays one row with `scrollWidth` within 900.

**Evidence:**

- Before/after fixture: `web/usability/dots-budget-point-chips.html`
- Fixture screenshot (900px): `web/usability/dots-budget-point-chips.png`

### Results / Layers still said Simplify (resolved)

**Where:** Results Scores first data column; Layers accordion summary; path toggles under that accordion; desktop `#paramsBar` `search r` chip.

**Problem:** After Compare glosses and Match error work, novices still saw a Scores column labeled `Simplify` next to DOTS / DP / SQUISH with no plain meaning, a Layers group also named `Simplify`, and path toggles `Simplified so far` / `Final simplified` with no tooltips. Params kept abbreviated `search r`. Skip-tour users could not tell which column was this run, or when to turn on the full green path.

**Fix:** Scores column is `This run` with a green-path tooltip. Layers accordion is `Paths & search`. Path toggles are `Path so far` / `Full result path` with tooltips. Params chip is `search radius`. User-facing “Simplify scores” copy is just “scores”. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/this-run-path-labels.html`
- Fixture screenshot (900px): `web/usability/this-run-path-labels.png`

### Params chips still said grid step / error budget / recorded (resolved)

**Where:** desktop `#paramsBar` secondary chips after Load; Results Scores `Time (ms)` row.

**Problem:** Match error / kept % were plain, but secondary chips still said `grid step`, `radius`, `error budget`, and `recorded`. Novices could not tell the orange file-saved value from the live Match error chip. Results also said `Time (ms)` while the header chip said `Time`.

**Fix:** Chips read `cell size`, `search r`, `max error`, and `file match` with tooltips. Results Scores uses `Time` with a milliseconds tooltip, matching the header.

**Evidence:**

- Before/after fixture: `web/usability/params-chip-plain-labels.html`
- Fixture screenshot (900px): `web/usability/params-chip-plain-labels.png`

### Status and playback Candidate wording disagreed (resolved)

**Where:** sidebar Status `path step` / `candidates` rows after Load; desktop `#playbackBar` Step / Segment / Candidate titles; mobile `#mobileTransport` aria-labels; Layers Anchor points tooltip.

**Problem:** Instructions and the playback tour already explained Step / Segment / Candidate in plain language, but Status still said `path step` and showed `candidates 3 / 12` (alive over anchor count). That ratio looks like “3 of 12 candidates” and does not match the dock. Candidate tooltips still said “search boundary”.

**Fix:** Status labels are `Step` / `Candidates` with `N still open`. Playback captions, button titles, and mobile aria-labels reuse the Instructions glosses (walk original points / jump simplified pieces / cycle next-point options). Anchor points drop “search boundary”. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/status-playback-gloss.html`
- Fixture screenshot (900px): `web/usability/status-playback-gloss.png`

### Post-load Instructions were shortcut-only (resolved)

**Where:** sidebar Instructions after Load (`#playbackInstructionsGloss`, desktop / mobile shortcut tables); View `#fitBtn`; Layers tour copy.

**Problem:** After Load, Instructions listed only keyboard / dock shortcuts. Novices who skipped the playback tour still saw Step / Segment / Candidate with no meaning. View also said `Fit to data`, which reads like a data action rather than resetting the map.

**Fix:** Instructions lead with a plain gloss for Step / Segment / Candidate / Play / Speed / Fit view (same ideas as the playback tour). Shortcut tables keep the keys underneath. The View button reads `Fit view`.

**Evidence:**

- Before/after fixture: `web/usability/playback-instructions-gloss.html`
- Fixture screenshot (900px): `web/usability/playback-instructions-gloss.png`

### Results Scores still said Compression and led with Fréchet (resolved)

**Where:** Results Scores table keep-share row; Match error footer / tooltips; upload success status; SQUISH keep validation.

**Problem:** The keep-share metric was labeled `Compression` even though the value is points kept (same idea as header `kept %`), so a novice reading `12.4%` could think only 12% was removed. Match error help still led with “discrete Fréchet distance”. Upload success said `✓ Generated`, and the SQUISH keep error used interval notation `(0, 100]`.

**Fix:** Scores row is `Kept %` with a plain tooltip. Match error footer / tooltips lead with drift-from-original wording and demote Fréchet to an optional technical link. Upload success is `✓ Loaded`. SQUISH validation says keep % must be greater than 0 and at most 100. Fixture `scrollWidth` stays within 900.

**Evidence:**

- Before/after fixture: `web/usability/scores-kept-pct.html`
- Fixture screenshot (900px): `web/usability/scores-kept-pct.png`

### Compare algorithm names and canvas “cur” stayed opaque (resolved)

**Where:** Results / header Compare pills (`DOTS` / `DP` / `SQUISH`); canvas overlay label beside the current path point; Status candidates tooltip; desktop keyboard shortcut for C/X; canvas loading HUD.

**Problem:** Compare pills were acronym-only. Desktop users could hover titles, but phones cannot, so novices who opened Results still did not know what DOTS / DP / SQUISH meant. The map still said `cur` while Status said `current point`. Shortcut copy still said `open candidates`, and the HUD said `Loading details…`.

**Fix:** Results pills show short visible glosses (`stream` / `classic` / `keep %`); the header strip hides those glosses so the Compare row stays one line. Canvas label is `current`. Status / Layers tooltips say “still being considered”. Shortcut text drops `open`. HUD reads `Loading trajectory…`.

**Evidence:**

- Before/after fixture: `web/usability/compare-algo-glosses.html`
- Fixture screenshot (900px): `web/usability/compare-algo-glosses.png`

### Mobile loaded header kept the long product name (resolved)

**Where:** mobile `#appTitle` / `header h1` after `body.trace-loaded-mobile` (Back + title + `?`).

**Problem:** After Load, phones hide the trajectory picker, but the header still showed the long product name `Trajectory Simplification Visualizer` with `white-space: nowrap`. Beside Back and `?` that title needed more width than a 390px row, so the name clipped or risked horizontal scroll, and novices could not see which trajectory was loaded.

**Fix:** Brand shortens to `Trajectory Simplifier`. On load / loading the header title becomes the selected trajectory label (or upload filename) with ellipsis (`min-width: 0`, `text-overflow: ellipsis`); failure or clear restores the brand. Fixture header width stays 390 with no page scroll.

**Evidence:**

- Before/after fixture: `web/usability/mobile-header-trajectory-title.html`
- Fixture screenshot (390px): `web/usability/mobile-header-trajectory-title.png`

### Status / canvas / Layers still led with paper notation (resolved)

**Where:** desktop `#statusIndices` labels; canvas `p` / `vN` overlays; Simplify layer toggle rows; loading / picker point-count copy.

**Problem:** Mobile Status already said `start point` / `current point`, but desktop still prefixed MathJax `p` / `vᵢ`, the map drew `p` / `v42`, and Layers kept symbols first (`Bp`, `Si[p]`, …). Novices could not match Status wording to the map, and loading still said `Loading pts…`.

**Fix:** Status uses plain `start point` / `current point` on every viewport. Canvas overlays say `start` / `current` in the same colors. Layer rows lead with plain glosses (paper symbols only in tooltips). Progress and picker counts say `points` / `Loading points…`; load failures say `Could not load trajectory`.

**Evidence:**

- Before/after fixture: `web/usability/status-canvas-plain-labels.html`
- Fixture screenshot (900px): `web/usability/status-canvas-plain-labels.png`

### Picker / load errors still said “trace” (resolved)

**Where:** mobile `#tracePicker` divider and unlabeled items; desktop `#traceSelect` fallback option labels; NDJSON load error strings; invalid JSON-upload `alert`; Results panel `aria-label`.

**Problem:** Start screens already said trajectory, but the picker divider still read `Other traces`, unlabeled ids fell back to `Trace N`, and failure alerts mentioned `simplify --web-server trace` / `Trace stream…`. Novices who skipped jargon elsewhere still hit mixed vocabulary on pick and error paths.

**Fix:** Divider is `Other trajectories`; unlabeled items / select options use `Trajectory N`. Load failures say plain “Loading failed” / “No data received…” / “Loading stopped…”. Invalid upload alert points at the plain-text format or preloaded list. Results `aria-label` is `Results scores and Compare`.

**Evidence:**

- Before/after fixture: `web/usability/picker-trajectory-wording.html`
- Fixture screenshot (900px): `web/usability/picker-trajectory-wording.png`

### Load Trace CTA still mixed “trace” with trajectory wording (resolved)

**Where:** `#loadBtn`, start Instructions / drop hint / tour, Compare status strings, mobile Back aria-label, and the params `trace error` chip.

**Problem:** Picker and Upload already said trajectory, but the primary CTA still read `Load Trace`, help copy repeated that label, Back said “trace selection”, and the params bar kept a `trace error` chip. Novices who skipped jargon elsewhere still hit mixed vocabulary on the same screen.

**Fix:** The button reads `Load` (shorter than `Load Trace`, with title / aria-label `Load trajectory and run simplification`). Instructions, tour, Compare status, and Results empty copy say `Load`. Back is “trajectory selection”. The params chip is `recorded` with the same Fréchet tooltip. Fixture row stays nowrap with no page overflow.

**Evidence:**

- Before/after fixture: `web/usability/load-button-trajectory.html`
- Fixture screenshot (900px): `web/usability/load-button-trajectory.png`

### Desktop empty-canvas Instructions showed playback keys before load (resolved)

**Where:** `#dropHint .desktop-instructions` on the empty canvas (desktop only); start picker `#preloadedLabel` / `#tracePicker` heading; loading status strings.

**Problem:** Before any trajectory loaded, desktop Instructions listed ←/→ Step, Segment, Candidate, and Space Play shortcuts. Novices who skipped the tour saw expert replay keys instead of how to start. The picker still said `Select trace…` / `Choose a preloaded trace` while Upload already said `trajectory`, and status used `Computing trace…` / `Loaded Trace N`.

**Fix:** Desktop empty Instructions now mirror the mobile start path (pick / upload → ε match / δ grid → Load Trace → Results / Compare), and point to the sidebar for keyboard shortcuts after load. Picker and user-facing copy say `trajectory`; status reads `Computing…` / `✓ Loaded <id>`.

**Evidence:**

- Before/after fixture: `web/usability/desktop-start-instructions.html`
- Fixture screenshot (900px): `web/usability/desktop-start-instructions.png`

### Results metrics heading said Compare while scores appear without it (resolved)

**Where:** Results panel `#resultsPanel` second section label; empty-canvas `#dropHint` headline.

**Problem:** After Load Trace, Simplify Match error / Time / points already fill the metrics table with no Compare algorithm selected. The section was still labeled `Compare`, so novices who only wanted scores thought they had to pick DOTS / DP / SQUISH first. The empty canvas headline also led with `Upload a trajectory…`, which underplayed the primary preloaded path (and on phones Upload is hidden).

**Fix:** The metrics block is labeled `Scores` with a tooltip that Compare adds columns when run. The drop hint reads `Choose a trajectory to get started` and lists preloaded before Upload.

**Evidence:**

- Before/after fixture: `web/usability/results-scores-label.html`
- Fixture screenshot (900px): `web/usability/results-scores-label.png`

### Compare looked available after Upload trajectory (resolved)

**Where:** header / Results Compare pills, `#baselineStatus`, `#baselineLayerHint`, start Instructions, and the load / Results tour steps after choosing Upload trajectory.

**Problem:** Compare only runs against preloaded traces (`/api/trace/<id>/compare`). After Upload trajectory, pills stayed clickable and status still said “Load a preloaded trace…”, which reads like the upload did not count even though Simplify scores already appear in Results.

**Fix:** Choosing an upload clears Compare selection, disables the DOTS / DP / SQUISH pills with a plain tooltip, and sets status / Layers hint / Results empty copy to “Compare needs a preloaded trace. Uploaded files show Simplify scores only.” Tour and mobile Instructions say the same. Picking a preloaded trace re-enables the pills.

**Evidence:**

- Before/after fixture: `web/usability/compare-upload-preloaded-only.html`
- Fixture screenshot (900px): `web/usability/compare-upload-preloaded-only.png`

### Cloud Run Deploy landed on main (resolved)

**Where:** `.github/workflows/deploy.yml` on `main`, GitHub Actions Deploy workflow, Cloud Run service `simplify-viewer`.

**Problem:** Deploy existed only on the usability branch, so GitHub listed only Benchmark and Correctness until a main merge.

**Fix:** [PR #11](https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/pull/11) merged to `main`. The push registered Deploy and published revision `simplify-viewer-00011-lzj` via Workload Identity Federation ([run 33982295165](https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/actions/runs/33982295165)). Live URL: https://simplify-viewer-522405269791.asia-east2.run.app (HTTP 200).

