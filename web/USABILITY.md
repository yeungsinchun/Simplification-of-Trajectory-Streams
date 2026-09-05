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

## Still open

No open layout or novice-copy items from this pass. First live Cloud Run publish still needs a successful `main` merge of [PR #11](https://github.com/yeungsinchun/Simplification-of-Trajectory-Streams/pull/11) (or a `workflow_dispatch` once Deploy is on `main`). Until then GitHub only lists Benchmark and Correctness.

