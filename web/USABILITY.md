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

## Still open

These were noticed while checking desktop and mobile and are not fixed here:

- On mobile, the playback dock appears as soon as Load Trace starts, while the canvas still has no geometry. The buttons are disabled, but Play still looks tappable over the empty canvas.