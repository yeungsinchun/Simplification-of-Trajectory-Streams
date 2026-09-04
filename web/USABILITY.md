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

## Still open

These were noticed while checking desktop and mobile and are not fixed here:

- Canvas `#traceLoadingHud` (`Loading details…`) is forced off with `display: none !important`, so after Load Trace the drop hint disappears and the canvas stays blank until geometry arrives. Progress is only in the header (upload status / metric slots).
- On the 390px start screen the δ number field can hide its value (default 500) because the input is `min-width: 0` in a two-column grid next to the spinner arrows.