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

## Still open

A 1280/900/390 loading pass now keeps the Load Trace control at its idle width. Remaining items for a later pass:

- After geometry arrives, the rest of the desktop param chips still wrap a second params row (header 79px to 110px at 1280, 113px to 144px at 900). The two metric chips stay in place.
- At 900px the wrapped \(S_i[p]\) second line can still paint about 10px past the window (MathJax `conv(G_{v_i})` box). A 1280-only pass misses it.
- Desktop playback-bar height still grows from 91px (1280) to about 108px (1024/900) because the bar wraps; treat wrap as a problem only if a control is clipped or untappable.