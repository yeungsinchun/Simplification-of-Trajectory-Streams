// ===========================================================================
//  Trajectory Simplification Visualizer — viewer logic
// ===========================================================================
//
// Consumes the JSON trace produced by `simplify --web-server` and renders an
// interactive step-through of the paper's streaming construction: for each
// "prefix" (one call to get_longest_stab) it shows the delta-disk candidate
// anchors P, and for each step within a prefix it shows the point v_i being
// consumed, the delta-disk hull G_i, each surviving candidate's free-space
// wedge F(S,p) fed into intersect(), the resulting stab region S, and the
// best (buffer) segment chosen so far.

(() => {
  "use strict";

  // -------------------------------------------------------------------------
  //  State
  // -------------------------------------------------------------------------

  const state = {
    trace: null,        // parsed JSON
    preview: null,      // lightweight original trajectory shown while the full trace loads
    prefixIdx: 0,
    stepIdx: 0,          // -1 means "before any step" (only P anchors, no consumption yet)
    candidateIdx: 0,     // which candidate to show F/S for (cycles through alive and dead)
    cam: { x: 0, y: 0, scale: 1 }, // world->screen: screen = (world - cam) * scale + center
    dragging: false,
    dragStart: null,
    camStart: null,
    playing: false,
    playTimer: null,
    computedFrechet: null,
    computingFrechet: false,
    currentStartPoint: null, // shared between renderStatus() and render() for the "p" label
  };

  // -------------------------------------------------------------------------
  //  DOM refs
  // -------------------------------------------------------------------------

  const el = (id) => document.getElementById(id);
  const canvas = el("canvas");
  const ctx = canvas.getContext("2d");
  const canvasWrap = el("canvasWrap");
  const canvasContainer = el("canvasContainer");
  const dropHint = el("dropHint");
  const paramsBar = el("paramsBar");
  const playbackBar = el("playbackBar");
  const statusGrid = el("statusGrid");
  const statusIndices = el("statusIndices");
  const traceLoadingHud = el("traceLoadingHud");
  const traceLoadingText = el("traceLoadingText");

  const uploadBtn = el("uploadBtn");
  const trajectoryInput = el("trajectoryInput");
  const trajectoryFileName = el("trajectoryFileName");
  const epsilonInput = el("epsilonInput");
  const deltaInput = el("deltaInput");
  const uploadStatus = el("uploadStatus");
  const traceSelect = el("traceSelect");
  const loadBtn = el("loadBtn");

  // Mobile off-canvas controls
  const headerToggle = el("headerToggle");
  const headerBody = el("headerBody");
  const playbackToggle = el("playbackToggle");
  const playbackBarEl = el("playbackBar");
  const mobileBackBtn = el("mobileBackBtn");
  const mobileLayersToggle = el("mobileLayersToggle");
  const layersSection = el("layersSection");
  const mobilePanelClose = el("mobilePanelClose");
  const mobilePanelOpen = el("mobilePanelOpen");
  const mobileTransport = el("mobileTransport");
  const mobileStepBackBtn = el("mobileStepBackBtn");
  const mobileCandidateBackBtn = el("mobileCandidateBackBtn");
  const mobilePlayBtn = el("mobilePlayBtn");
  const mobileCandidateForwardBtn = el("mobileCandidateForwardBtn");
  const mobileStepForwardBtn = el("mobileStepForwardBtn");

  const stepInput = el("stepInput");
  const segmentInput = el("segmentInput");
  const playBtn = el("playBtn");
  const frechetResult = el("frechetResult");

  // Current loaded file/trace state
  let currentFile = null;
  let currentTraceId = null;

  function setLoadButtonBusy(isBusy) {
    loadBtn.disabled = isBusy;
    loadBtn.setAttribute("aria-busy", String(isBusy));
    loadBtn.innerHTML = isBusy
      ? '<span class="button-spinner" aria-hidden="true"></span><span class="visually-hidden">Loading</span>'
      : "Load Trace";
  }

  // Speed presets
  const speedPresets = document.querySelectorAll(".speed-preset");
  let currentSpeedMultiplier = 1; // Default 1x
  const toggles = {};
  for (const key of ["stream", "simplified", "ball-p0", "ball-pi", "Gi", "F", "F-Si", "S", "P", "dead-candidates"]) {
    toggles[key] = el(`toggle-${key}`);
  }

  // -------------------------------------------------------------------------
  //  Play button state helper (defined early)
  // -------------------------------------------------------------------------

  function updatePlayButton() {
    const playIcon = playBtn.querySelector('.play-icon');
    const pauseIcon = playBtn.querySelector('.pause-icon');
    const label = playBtn.querySelector('.play-label');
    if (state.playing) {
      playIcon.style.display = 'none';
      pauseIcon.style.display = 'block';
      if (label) label.textContent = 'Pause';
    } else {
      playIcon.style.display = 'block';
      pauseIcon.style.display = 'none';
      if (label) label.textContent = 'Play';
    }
    // Mirror to the mobile FAB
    const fabPlay = playbackToggle && playbackToggle.querySelector('.play-fab-icon');
    const fabPause = playbackToggle && playbackToggle.querySelector('.pause-fab-icon');
    if (fabPlay && fabPause) {
      fabPlay.style.display = state.playing ? 'none' : 'block';
      fabPause.style.display = state.playing ? 'block' : 'none';
    }
    const mobilePlay = mobilePlayBtn && mobilePlayBtn.querySelector('.mobile-play-icon');
    const mobilePause = mobilePlayBtn && mobilePlayBtn.querySelector('.mobile-pause-icon');
    const mobileLabel = mobilePlayBtn && mobilePlayBtn.querySelector('.mobile-play-label');
    if (mobilePlay && mobilePause) {
      mobilePlay.style.display = state.playing ? 'none' : 'block';
      mobilePause.style.display = state.playing ? 'block' : 'none';
      mobilePlayBtn.setAttribute('aria-label', state.playing ? 'Pause' : 'Play');
      if (mobileLabel) mobileLabel.textContent = state.playing ? 'Pause' : 'Play';
    }
  }
  
  // Initialize play button state
  updatePlayButton();

  // -------------------------------------------------------------------------
  //  File loading
  // -------------------------------------------------------------------------

  async function apiErrorMessage(resp) {
    const fallback = `Server error ${resp.status}`;
    try {
      const text = await resp.text();
      if (!text) return fallback;
      try {
        const data = JSON.parse(text);
        return data.error || fallback;
      } catch {
        return text.length <= 200 ? text : fallback;
      }
    } catch {
      return fallback;
    }
  }

  function trajectoryPointsFromText(text) {
    const lines = text.trim().split(/\r?\n/);
    const count = Number.parseInt(lines[0], 10);
    if (!Number.isFinite(count) || count < 2) return [];
    const points = [];
    for (let i = 1; i < lines.length && points.length < count; i++) {
      const [x, y] = lines[i].trim().split(/\s+/).map(Number);
      if (Number.isFinite(x) && Number.isFinite(y)) points.push([x, y]);
    }
    return points;
  }

  function previewBBox(points, epsilon, delta, yOffset) {
    let minCoord = Infinity;
    let maxCoord = -Infinity;
    for (const [x, y] of points) {
      minCoord = Math.min(minCoord, x, y);
      maxCoord = Math.max(maxCoord, x, y);
    }
    const gridReach = (1 + epsilon) * delta;
    const padding = gridReach + Math.max(1, gridReach * 1e-6);
    return [
      minCoord - padding,
      minCoord - padding + yOffset,
      maxCoord + padding,
      maxCoord + padding + yOffset,
    ];
  }

  function showTracePreview(points) {
    if (!Array.isArray(points) || points.length < 2 || state.trace) return;

    const isSampleTrace = currentTraceId !== null
      && [51, 52, 53].includes(Number.parseInt(currentTraceId, 10));
    const yOffset = isSampleTrace ? 500 : 0;
    const previewPoints = isSampleTrace
      ? points.map(([x, y]) => [x, y + yOffset])
      : points;
    const epsilon = Number.parseFloat(epsilonInput.value);
    const delta = Number.parseFloat(deltaInput.value);

    state.preview = {
      points: previewPoints,
      bbox: previewBBox(points, epsilon, delta, yOffset),
    };
    dropHint.style.display = "none";
    el("sidebar").style.display = "flex";
    el("resizer").style.display = "block";
    canvas.classList.add("has-trace");
    document.body.classList.add("trace-loaded-mobile", "trace-preview-loading");
    document.body.classList.remove("trace-preview-error");
    traceLoadingHud.classList.add("visible");
    traceLoadingHud.classList.remove("error");
    traceLoadingText.textContent = "Loading details…";
    fitToBBox(state.preview.bbox);
    render();
  }

  function finishTracePreview() {
    state.preview = null;
    document.body.classList.remove("trace-preview-loading", "trace-preview-error");
    traceLoadingHud.classList.remove("visible", "error");
  }

  function failTracePreview(message) {
    if (!state.preview) return;
    document.body.classList.remove("trace-preview-loading");
    document.body.classList.add("trace-preview-error");
    traceLoadingHud.classList.add("visible", "error");
    traceLoadingText.textContent = message;
  }

  async function showPreloadedTracePreview(traceId) {
    try {
      const resp = await fetch(`/api/trace/${traceId}/original`);
      if (!resp.ok) return;
      const data = await resp.json();
      showTracePreview(data.points);
    } catch (err) {
      console.warn("Could not load trajectory preview:", err);
    }
  }

  function loadTraceText(text, name) {
    let parsed;
    try {
      const sizeKB = (text.length / 1024).toFixed(0);
      console.log(`[Client] Starting JSON.parse of ${sizeKB}KB string...`);
      uploadStatus.textContent = `Parsing ${sizeKB}KB JSON...`;
      uploadStatus.style.color = "#e8c547";
      
      const parseStart = performance.now();
      parsed = JSON.parse(text);
      const parseTime = ((performance.now() - parseStart) / 1000).toFixed(2);
      console.log(`[Client] JSON.parse completed in ${parseTime}s`);
    } catch (e) {
      alert("Could not parse JSON: " + e.message);
      uploadStatus.textContent = `Error: ${e.message}`;
      uploadStatus.style.color = "#ff5f6d";
      return;
    }
    if (!parsed || !Array.isArray(parsed.prefixes)) {
      alert("This does not look like a simplify --web-server trace (missing 'prefixes').");
      uploadStatus.textContent = "Invalid trace format";
      uploadStatus.style.color = "#ff5f6d";
      return;
    }
    stopPlaying();
    
    // For sample traces (51, 52, 53), shift all y-coordinates up by 500
    const isSampleTrace = (currentTraceId !== null)
      && (parseInt(currentTraceId, 10) === 51
       || parseInt(currentTraceId, 10) === 52
       || parseInt(currentTraceId, 10) === 53);
    if (isSampleTrace) {
      const Y_OFFSET = 500;
      // Shift bbox
      if (parsed.bbox) {
        parsed.bbox[1] += Y_OFFSET; // ymin
        parsed.bbox[3] += Y_OFFSET; // ymax
      }
      // Shift stream (full input trajectory)
      if (parsed.stream) {
        parsed.stream.forEach(p => p[1] += Y_OFFSET);
      }
      // Shift all points in prefixes
      if (parsed.prefixes) {
        parsed.prefixes.forEach(pfx => {
          // Shift P (original points)
          if (pfx.P) pfx.P.forEach(p => p[1] += Y_OFFSET);
          // Shift p0 (start anchor)
          if (pfx.p0) pfx.p0[1] += Y_OFFSET;
          // Shift S (simplified so far)
          if (pfx.S) pfx.S.forEach(p => p[1] += Y_OFFSET);
          // Shift output (the simplified segment for this prefix)
          if (pfx.output) {
            pfx.output.forEach(p => p[1] += Y_OFFSET);
          }
          // Shift all step data
          if (pfx.steps) {
            pfx.steps.forEach(step => {
              // Shift pi (current point)
              if (step.pi) step.pi[1] += Y_OFFSET;
              // Shift Gi (convex hull)
              if (step.Gi) step.Gi.forEach(p => p[1] += Y_OFFSET);
              // Shift candidate data
              if (step.candidates) {
                step.candidates.forEach(cand => {
                  if (cand.F) cand.F.forEach(p => p[1] += Y_OFFSET);
                  if (cand.S) cand.S.forEach(p => p[1] += Y_OFFSET);
                  if (cand.rays) cand.rays.forEach(ray => ray.forEach(p => p[1] += Y_OFFSET));
                });
              }
            });
          }
        });
      }
    }
    
    finishTracePreview();
    state.trace = parsed;
    state.prefixIdx = 0;
    state.stepIdx = parsed.prefixes.length && parsed.prefixes[0].steps.length ? 0 : -1;
    state.computedFrechet = null;
    state.computingFrechet = false;
    dropHint.style.display = "none";
    el("sidebar").style.display = "flex";
    el("resizer").style.display = "block";
    canvas.classList.add("has-trace");
    document.body.classList.add("trace-loaded-mobile");
    if (playbackBar) playbackBar.classList.add("visible");
    if (playbackToggle) playbackToggle.classList.add("visible");
    renderParamsBar();
    setupSliders();
    fitToData();
    render();
    // Update button states - disable step back if we're at step 0 (i=1, the first step)
    const stepBackBtn = el("stepBackBtn");
    stepBackBtn.disabled = (state.prefixIdx === 0 && state.stepIdx <= 0);
    // Typeset MathJax in sidebar
    if (window.MathJax) {
      MathJax.typesetPromise([el("sidebar")]).catch(() => {});
    }
    // Automatically start Fréchet computation
    startFrechetComputation();
  }

  // Choose trajectory file and auto-generate trace
  uploadBtn.addEventListener("click", () => trajectoryInput.click());
  trajectoryInput.addEventListener("change", async (e) => {
    const f = e.target.files && e.target.files[0];
    trajectoryFileName.textContent = f ? f.name : "no file chosen";
    
    if (!f) return;
    
    currentFile = f;
    currentTraceId = null;
    traceSelect.value = "";
    syncPreloadedTrigger();
    uploadStatus.textContent = "File selected. Click 'Load Trace' to generate.";
    uploadStatus.style.color = "#8b93a3";
  });

  // Load button handler - generates trace from current file or preloaded trace
  loadBtn.addEventListener("click", async () => {
    const eps = parseFloat(epsilonInput.value);
    const delta = parseFloat(deltaInput.value);
    if (isNaN(eps) || eps <= 0 || isNaN(delta) || delta <= 0) {
      alert("Invalid epsilon or delta");
      return;
    }

    // Collapse header drawer once loading starts so the canvas gets more room.
    if (headerBody && headerBody.classList.contains("open")) closeHeader();

    // Check if we have a file or a selected trace
    if (currentFile) {
      // Generate from uploaded file
      uploadStatus.textContent = "";
      setLoadButtonBusy(true);

      try {
        const formData = new FormData();
        formData.append("file", currentFile);
        formData.append("epsilon", eps);
        formData.append("delta", delta);

        console.log(`[Client] Starting trace upload with eps=${eps}, delta=${delta}`);
        const startTime = performance.now();

        const traceRequest = fetch("/api/trace", { method: "POST", body: formData });
        const originalText = await currentFile.text();
        showTracePreview(trajectoryPointsFromText(originalText));

        const resp = await traceRequest;
        if (!resp.ok) {
          throw new Error(await apiErrorMessage(resp));
        }
        const traceText = await resp.text();
        
        const elapsed = ((performance.now() - startTime) / 1000).toFixed(2);
        console.log(`[Client] Trace upload completed in ${elapsed}s`);
        
        loadTraceText(traceText, `${currentFile.name} (ε=${eps}, δ=${delta})`);
        uploadStatus.textContent = "✓ Generated";
        uploadStatus.style.color = "#3ddc97";
      } catch (err) {
        uploadStatus.textContent = `Error: ${err.message}`;
        uploadStatus.style.color = "#ff5f6d";
        failTracePreview("Could not load details");
      } finally {
        setLoadButtonBusy(false);
      }
    } else if (currentTraceId) {
      // Generate from preloaded trace
      uploadStatus.textContent = "";
      setLoadButtonBusy(true);

      try {
        console.log(`[Client] Starting trace load (ID=${currentTraceId}) with eps=${eps}, delta=${delta}`);
        const startTime = performance.now();

        const traceRequest = fetch(`/api/trace/${currentTraceId}?epsilon=${eps}&delta=${delta}`);
        await showPreloadedTracePreview(currentTraceId);

        const resp = await traceRequest;
        if (!resp.ok) {
          throw new Error(await apiErrorMessage(resp));
        }
        const traceText = await resp.text();
        
        const elapsed = ((performance.now() - startTime) / 1000).toFixed(2);
        console.log(`[Client] Trace load completed in ${elapsed}s`);
        
        loadTraceText(traceText, `Trace ${currentTraceId} (ε=${eps}, δ=${delta})`);
        uploadStatus.textContent = `✓ Loaded Trace ${currentTraceId}`;
        uploadStatus.style.color = "#3ddc97";
      } catch (err) {
        uploadStatus.textContent = `Error: ${err.message}`;
        uploadStatus.style.color = "#ff5f6d";
        failTracePreview("Could not load details");
      } finally {
        setLoadButtonBusy(false);
      }
    } else {
      uploadStatus.textContent = "Please select a file or preloaded trace first";
      uploadStatus.style.color = "#ff5f6d";
    }
  });


  // -------------------------------------------------------------------------
  //  Preloaded trace dropdown
  // -------------------------------------------------------------------------

  // Trace metadata keyed by id — populated on load, used by the change handler.
  const traceMetaById = {};

  // Cached ordered trace list (id + label + n_points) used by the picker overlay.
  let tracesList = [];

  const preloadedTrigger = el("preloadedTrigger");
  const preloadedLabel = el("preloadedLabel");
  const tracePicker = el("tracePicker");
  const tracePickerBackdrop = el("tracePickerBackdrop");
  const tracePickerClose = el("tracePickerClose");
  const tracePickerList = el("tracePickerList");

  function openTracePicker() {
    if (!tracePicker) return;
    tracePicker.classList.add("open");
    tracePickerBackdrop.classList.add("open");
    tracePicker.setAttribute("aria-hidden", "false");
    document.body.style.overflow = "hidden";
  }
  function closeTracePicker() {
    if (!tracePicker) return;
    tracePicker.classList.remove("open");
    tracePickerBackdrop.classList.remove("open");
    tracePicker.setAttribute("aria-hidden", "true");
    document.body.style.overflow = "";
  }

  function renderTracePicker() {
    if (!tracePickerList) return;
    tracePickerList.innerHTML = "";
    if (!tracesList.length) {
      const empty = document.createElement("div");
      empty.className = "trace-picker-empty";
      empty.textContent = "No preloaded traces available.";
      tracePickerList.appendChild(empty);
      return;
    }
    let addedDivider = false;
    tracesList.forEach(t => {
      const id = t.id !== undefined ? t.id : t;
      const n = t.n_points;
      const label = t.label || `Trace ${id}`;
      const isSample = (id === 51 || id === 52 || id === 53);
      if (!isSample && !addedDivider) {
        const sep = document.createElement("div");
        sep.className = "trace-picker-divider";
        sep.textContent = "Other traces";
        tracePickerList.appendChild(sep);
        addedDivider = true;
      }
      const btn = document.createElement("button");
      btn.type = "button";
      btn.className = "trace-picker-item";
      btn.dataset.traceId = id;
      if (String(traceSelect.value) === String(id)) btn.classList.add("selected");
      const labelEl = document.createElement("span");
      labelEl.textContent = label;
      btn.appendChild(labelEl);
      if (n != null) {
        const meta = document.createElement("span");
        meta.className = "trace-picker-item-meta";
        meta.textContent = `${n.toLocaleString()} pts`;
        btn.appendChild(meta);
      }
      btn.addEventListener("click", () => {
        traceSelect.value = id;
        traceSelect.dispatchEvent(new Event("change"));
        closeTracePicker();
      });
      tracePickerList.appendChild(btn);
    });
  }

  function syncPreloadedTrigger() {
    if (!preloadedTrigger || !preloadedLabel) return;
    // If a file is currently chosen, show its name in the trigger.
    if (currentFile && currentFile.name) {
      preloadedLabel.textContent = currentFile.name;
      preloadedTrigger.classList.add("has-value");
      return;
    }
    const opt = traceSelect.options[traceSelect.selectedIndex];
    const value = traceSelect.value;
    if (value && opt && opt.textContent) {
      preloadedLabel.textContent = opt.textContent;
      preloadedTrigger.classList.add("has-value");
    } else {
      preloadedLabel.textContent = "Select trace…";
      preloadedTrigger.classList.remove("has-value");
    }
  }

  if (preloadedTrigger) preloadedTrigger.addEventListener("click", openTracePicker);
  if (tracePickerClose) tracePickerClose.addEventListener("click", closeTracePicker);
  if (tracePickerBackdrop) tracePickerBackdrop.addEventListener("click", closeTracePicker);
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      if (tracePicker && tracePicker.classList.contains("open")) closeTracePicker();
      if (playbackBarEl && playbackBarEl.classList.contains("open")) closePlayback();
      if (headerBody && headerBody.classList.contains("open")) closeHeader();
    }
  });

  // ---- Mobile header toggle (collapse file controls) ----
  function openHeader() {
    if (!headerBody) return;
    headerBody.classList.add("open");
    headerToggle.classList.add("active");
    headerToggle.setAttribute("aria-expanded", "true");
  }
  function closeHeader() {
    if (!headerBody) return;
    headerBody.classList.remove("open");
    headerToggle.classList.remove("active");
    headerToggle.setAttribute("aria-expanded", "false");
  }
  if (headerToggle) {
    headerToggle.addEventListener("click", () => {
      if (headerBody.classList.contains("open")) closeHeader();
      else openHeader();
    });
  }

  // ---- Mobile playback FAB + bottom-sheet ----
  function openPlayback() {
    if (!playbackBarEl) return;
    playbackBarEl.classList.add("open");
    playbackToggle.classList.add("active");
    playbackToggle.setAttribute("aria-expanded", "true");
    document.body.classList.add("playback-open");
  }
  function closePlayback() {
    if (!playbackBarEl) return;
    playbackBarEl.classList.remove("open");
    playbackToggle.classList.remove("active");
    playbackToggle.setAttribute("aria-expanded", "false");
    document.body.classList.remove("playback-open");
  }
  if (playbackToggle) {
    playbackToggle.addEventListener("click", () => {
      if (playbackBarEl.classList.contains("open")) closePlayback();
      else openPlayback();
    });
  }
  if (mobileBackBtn) {
    mobileBackBtn.addEventListener("click", () => window.location.reload());
  }
  if (mobileLayersToggle && layersSection) {
    mobileLayersToggle.addEventListener("click", () => {
      const isOpen = layersSection.classList.toggle("layers-open");
      mobileLayersToggle.setAttribute("aria-expanded", String(isOpen));
    });
  }
  if (mobilePanelClose) {
    mobilePanelClose.addEventListener("click", () => {
      document.body.classList.add("mobile-panel-closed");
      resizeCanvas();
      render();
    });
  }
  if (mobilePanelOpen) {
    mobilePanelOpen.addEventListener("click", () => {
      document.body.classList.remove("mobile-panel-closed");
      resizeCanvas();
      render();
    });
  }

  // Keep the touch bar anchored to the visible screen while mobile browsers
  // pinch-zoom and pan their visual viewport.
  function syncMobileTransportViewport() {
    if (!mobileTransport || !window.visualViewport) return;
    const viewport = window.visualViewport;
    if (viewport.scale <= 1.01) {
      mobileTransport.style.transform = "";
      return;
    }
    const x = viewport.offsetLeft;
    const y = viewport.offsetTop + viewport.height - document.documentElement.clientHeight;
    mobileTransport.style.transform =
      `translate(${x}px, ${y}px) scale(${1 / viewport.scale})`;
  }
  if (window.visualViewport) {
    window.visualViewport.addEventListener("resize", syncMobileTransportViewport);
    window.visualViewport.addEventListener("scroll", syncMobileTransportViewport);
  }
  window.addEventListener("resize", syncMobileTransportViewport);
  syncMobileTransportViewport();

  if (mobileStepBackBtn) {
    mobileStepBackBtn.addEventListener("click", () => el("stepBackBtn").click());
  }
  if (mobileCandidateBackBtn) {
    mobileCandidateBackBtn.addEventListener("click", () => el("candidateBackBtn").click());
  }
  if (mobilePlayBtn) {
    mobilePlayBtn.addEventListener("click", () => playBtn.click());
  }
  if (mobileCandidateForwardBtn) {
    mobileCandidateForwardBtn.addEventListener("click", () => el("candidateForwardBtn").click());
  }
  if (mobileStepForwardBtn) {
    mobileStepForwardBtn.addEventListener("click", () => el("stepForwardBtn").click());
  }

  // Populate dropdown on page load
  (async () => {
    try {
      console.log('[Traces] Fetching trace list from /api/traces...');
      const resp = await fetch("/api/traces");
      console.log('[Traces] Response status:', resp.status, resp.statusText);
      if (resp.ok) {
        const data = await resp.json();
        console.log('[Traces] Received data:', data);
        if (data.traces && data.traces.length > 0) {
          tracesList = data.traces;
          let addedDivider = false;
          data.traces.forEach(t => {
            const id = t.id !== undefined ? t.id : t;
            const n  = t.n_points;
            traceMetaById[id] = t;

            // Insert a visual separator between sample traces and regular ones.
            const isSample = (id === 51 || id === 52 || id === 53);
            if (!isSample && !addedDivider) {
              const sep = document.createElement("option");
              sep.disabled = true;
              sep.textContent = "──────────────────";
              traceSelect.appendChild(sep);
              addedDivider = true;
            }

            const opt = document.createElement("option");
            opt.value = id;
            if (t.label) {
              opt.textContent = t.label;
            } else {
              opt.textContent = n != null
                ? `Trace ${id}  (${n.toLocaleString()} pts)`
                : `Trace ${id}`;
            }
            traceSelect.appendChild(opt);
          });
          renderTracePicker();
          console.log(`[Traces] Added ${data.traces.length} traces to dropdown`);
        } else {
          console.warn('[Traces] No traces in response');
        }
      } else {
        console.error('[Traces] Failed to fetch:', resp.status, resp.statusText);
      }
    } catch (err) {
      console.error("Failed to load trace list:", err);
    }
  })();

  // Load trace automatically when dropdown selection changes
  traceSelect.addEventListener("change", async () => {
    const traceId = traceSelect.value;
    if (!traceId) {
      // Clear trace if deselected
      state.trace = null;
      currentFile = null;
      currentTraceId = null;
      render();
      uploadStatus.textContent = "";
      dropHint.style.display = "flex";
      document.body.classList.remove("trace-loaded-mobile");
      if (playbackBar) playbackBar.classList.remove("visible");
      if (playbackToggle) playbackToggle.classList.remove("visible");
      syncPreloadedTrigger();
      return;
    }

    currentTraceId = traceId;
    currentFile = null;
    trajectoryFileName.textContent = "no file chosen";
    syncPreloadedTrigger();
    if (headerBody && headerBody.classList.contains("open")) closeHeader();

    // Auto-fill ε/δ from meta if available (sample traces have recommended values).
    const meta = traceMetaById[traceId];
    if (meta && meta.epsilon != null) epsilonInput.value = meta.epsilon;
    if (meta && meta.delta   != null) deltaInput.value   = meta.delta;

    const name = (meta && meta.label) ? meta.label : `Trace ${traceId}`;
    uploadStatus.textContent = `${name} selected. Click 'Load Trace' to generate.`;
    uploadStatus.style.color = "#8b93a3";
  });

  function computeBBox(points) {
    if (!points || points.length === 0) return [0, 0, 100, 100];
    let xmin = Infinity, ymin = Infinity, xmax = -Infinity, ymax = -Infinity;
    for (const [x, y] of points) {
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
    }
    return [xmin, ymin, xmax, ymax];
  }


  // -------------------------------------------------------------------------
  //  Params bar / status grid
  // -------------------------------------------------------------------------

  function renderParamsBar() {
    const t = state.trace;
    if (!t) { paramsBar.innerHTML = ""; return; }
    const fmt = (v) => (typeof v === "number" ? (Math.abs(v) >= 1000 ? v.toFixed(1) : v.toFixed(4)) : v);
    let computedFrechetDisplay = '';
    if (state.computedFrechet != null) {
      console.log('Rendering computed Fréchet:', state.computedFrechet);
      computedFrechetDisplay = `<span style="color: #38BDF8; font-weight: 600;">Computed Fréchet distance <b style="color: #22D3EE;">${fmt(state.computedFrechet)}</b></span>`;
    } else if (state.computingFrechet) {
      computedFrechetDisplay = `<span style="color:var(--text-dim)">Computed Fréchet distance <b>computing…</b></span>`;
    }
    const frechetDisplay = t.frechet_distance != null 
      ? `<span style="color: #C4612F; font-weight: 600;">Actual Fr&eacute;chet distance <b style="color: #A94E22;">${fmt(t.frechet_distance)}</b></span>` 
      : '';
    paramsBar.innerHTML = `
      <span>\\(\\varepsilon\\) <b>${fmt(t.eps)}</b></span>
      <span>\\(\\delta\\) <b>${fmt(t.delta)}</b></span>
      <span>\\(\\text{len}_\\text{grid}\\) <b>${fmt(t.grid_val)}</b></span>
      <span>\\(R\\) (disk radius) <b>${fmt(t.r_val)}</b></span>
      <span>a-priori Fréchet bound <b>${fmt(t.expected_frechet)}</b></span>
      ${computedFrechetDisplay}
      ${frechetDisplay}
      <span style="color: #38BDF8; font-weight: 600;">Simplification time <b style="color: #22D3EE;">${fmt(t.time_ms)} ms</b></span>
      <span>|stream| <b>${t.stream.length}</b></span>
      <span>|simplified| <b>${t.simplified.length}</b></span>
      <span>ratio <b>${(100 * t.simplified.length / t.stream.length).toFixed(1)}%</b></span>
    `;
    if (window.MathJax) {
      MathJax.typesetPromise([paramsBar]).catch(() => {});
    }
  }

  function currentPrefix() {
    if (!state.trace) return null;
    return state.trace.prefixes[state.prefixIdx] || null;
  }

  function currentStep() {
    const pfx = currentPrefix();
    if (!pfx || state.stepIdx < 0) return null;
    return pfx.steps[state.stepIdx] || null;
  }

  function renderStatus() {
    const pfx = currentPrefix();
    const step = currentStep();
    if (!pfx) {
      statusIndices.innerHTML = "";
      statusGrid.innerHTML = "";
      return;
    }

    const curIdx   = step ? step.stream_idx : "—";

    // p always shows the stream index of the prefix's anchor vertex, NOT the
    // candidate grid index.  p0_idx is set once when the prefix is created
    // (it is the index into the full stream where this prefix begins).
    let startPoint = pfx.p0;
    let startIdx = pfx.p0_idx !== undefined ? pfx.p0_idx : "—";
    state.currentStartPoint = startPoint;

    // Big index numbers
    statusIndices.innerHTML = `
      <div class="status-idx-block">
        <span class="idx-label" style="color:#ff9f43">\\(p\\) (start)</span>
        <span class="idx-num" style="color:#ff9f43">#${startIdx}</span>
        <span class="idx-coord">${ptStr(startPoint)}</span>
      </div>
      <div class="status-idx-block">
        <span class="idx-label" style="color:#ff7ae8">\\(v_i\\) (current)</span>
        <span class="idx-num" style="color:#ff7ae8">${step ? "#" + curIdx : "—"}</span>
        <span class="idx-coord">${step ? ptStr(step.pi) : ""}</span>
      </div>`;
    if (window.MathJax) {
      MathJax.typesetPromise([statusIndices]).catch(() => {});
    }

    // Detail rows — present-state only, no future end vertex
    const rows = [];
    if (step) {
      const alive = step.candidates.filter((c) => c.alive).length;
      rows.push(["alive", `<b style="color:#3ddc97">${alive}</b> / ${pfx.P.length}`]);
    } else {
      rows.push(["step", "before first"]);
      rows.push(["alive", `<b style="color:#3ddc97">${pfx.P.length}</b> / ${pfx.P.length}`]);
    }

    statusGrid.innerHTML = rows
      .map(([k, v]) => `<span>${k}</span><span class="mono"><b>${v}</b></span>`)
      .join("");
    if (window.MathJax) {
      MathJax.typesetPromise([statusGrid]).catch(() => {});
    }
  }

  function ptStr(p) {
    if (!p) return "-";
    return `(${p[0].toFixed(1)}, ${p[1].toFixed(1)})`;
  }

  // -------------------------------------------------------------------------
  //  Camera / coordinate transform
  //
  //  World Y grows "up" in the trajectory data (plain Cartesian); canvas Y
  //  grows down, so we flip Y when projecting.
  // -------------------------------------------------------------------------

  function resizeCanvas() {
    const rect = canvasContainer.getBoundingClientRect();
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.round(rect.width * dpr));
    canvas.height = Math.max(1, Math.round(rect.height * dpr));
    canvas.style.width = rect.width + "px";
    canvas.style.height = rect.height + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  }

  function worldToScreen(x, y) {
    const rect = canvasContainer.getBoundingClientRect();
    const sx = (x - state.cam.x) * state.cam.scale + rect.width / 2;
    const sy = -(y - state.cam.y) * state.cam.scale + rect.height / 2;
    return [sx, sy];
  }

  function screenToWorld(sx, sy) {
    const rect = canvasContainer.getBoundingClientRect();
    const x = (sx - rect.width / 2) / state.cam.scale + state.cam.x;
    const y = -(sy - rect.height / 2) / state.cam.scale + state.cam.y;
    return [x, y];
  }

  function fitToBBox(bbox) {
    if (!bbox) return;
    let xmin = bbox[0], ymin = bbox[1], xmax = bbox[2], ymax = bbox[3];
    if (![xmin, ymin, xmax, ymax].every(Number.isFinite)) {
      xmin = -1; ymin = -1; xmax = 1; ymax = 1;
    }
    const rect = canvasContainer.getBoundingClientRect();
    const w = Math.max(1e-6, xmax - xmin);
    const h = Math.max(1e-6, ymax - ymin);
    const pad = 0.08;
    const scale = Math.min(rect.width / (w * (1 + pad)), rect.height / (h * (1 + pad)));
    state.cam.scale = scale > 0 && isFinite(scale) ? scale : 1;
    state.cam.x = (xmin + xmax) / 2;
    state.cam.y = (ymin + ymax) / 2;
  }

  function fitToData() {
    const data = state.trace || state.preview;
    if (!data) return;
    fitToBBox(data.bbox);
  }

  el("fitBtn").addEventListener("click", () => { fitToData(); render(); });

  // Pan (drag).
  canvas.addEventListener("mousedown", (e) => {
    if (!state.trace && !state.preview) return;
    state.dragging = true;
    canvas.classList.add("dragging");
    state.dragStart = [e.clientX, e.clientY];
    state.camStart = { x: state.cam.x, y: state.cam.y };
  });
  window.addEventListener("mousemove", (e) => {
    if (!state.dragging) return;
    const dx = e.clientX - state.dragStart[0];
    const dy = e.clientY - state.dragStart[1];
    state.cam.x = state.camStart.x - dx / state.cam.scale;
    state.cam.y = state.camStart.y + dy / state.cam.scale;
    render();
  });
  window.addEventListener("mouseup", () => {
    state.dragging = false;
    canvas.classList.remove("dragging");
  });

  // Zoom (wheel), anchored at the cursor.
  canvasWrap.addEventListener(
    "wheel",
    (e) => {
      if (!state.trace && !state.preview) return;
      e.preventDefault();
      const rect = canvasWrap.getBoundingClientRect();
      const sx = e.clientX - rect.left;
      const sy = e.clientY - rect.top;
      const before = screenToWorld(sx, sy);
      const factor = Math.exp(-e.deltaY * 0.0015);
      state.cam.scale = Math.min(1e6, Math.max(1e-6, state.cam.scale * factor));
      const after = screenToWorld(sx, sy);
      state.cam.x += before[0] - after[0];
      state.cam.y += before[1] - after[1];
      render();
    },
    { passive: false }
  );

  window.addEventListener("resize", () => { resizeCanvas(); render(); });
  window.addEventListener("orientationchange", () => { setTimeout(() => { resizeCanvas(); render(); }, 100); });

  // -------------------------------------------------------------------------
  //  Prefix / step navigation
  // -------------------------------------------------------------------------

  function setupSliders() {
    // Sliders removed; labels only
  }

  function updateStepSliderBounds() {
    // Sliders removed; labels only
  }

  function clampPrefix(i) {
    const t = state.trace;
    if (!t) return 0;
    return Math.max(0, Math.min(t.prefixes.length - 1, i));
  }

  function goToPrefix(i, keepStep) {
    if (!state.trace) return;
    state.prefixIdx = clampPrefix(i);
    const pfx = currentPrefix();
    if (!keepStep) state.stepIdx = pfx && pfx.steps.length ? 0 : -1;
    else state.stepIdx = pfx ? Math.max(-1, Math.min(pfx.steps.length - 1, state.stepIdx)) : -1;
    state.candidateIdx = 0; // Reset candidate cycling when changing prefix
    render();
    // Update button states - disable step back if at prefix 0 and step 0 or before
    const stepBackBtn = el("stepBackBtn");
    stepBackBtn.disabled = (state.prefixIdx === 0 && state.stepIdx <= 0);
  }

  function goToStep(i) {
    const pfx = currentPrefix();
    if (!pfx) return;
    const n = pfx.steps.length;

    // Going forward past the last step — advance to next prefix.
    if (i > n - 1) {
      if (state.prefixIdx < state.trace.prefixes.length - 1) {
        goToPrefix(state.prefixIdx + 1, false);
      }
      return;
    }

    // Going backward past the first step.
    if (i < -1) {
      if (state.prefixIdx > 0) {
        goToPrefix(state.prefixIdx - 1, true);
        const prevPfx = currentPrefix();
        state.stepIdx = prevPfx.steps.length - 1;
        state.candidateIdx = 0;
        render();
        // Update button states
        const stepBackBtn = el("stepBackBtn");
        stepBackBtn.disabled = (state.prefixIdx === 0 && state.stepIdx <= 0);
      }
      return;
    }

    state.stepIdx = i;
    state.candidateIdx = 0;
    render();
    
    // Update button states - disable step back if at prefix 0 and step 0 or before
    const stepBackBtn = el("stepBackBtn");
    stepBackBtn.disabled = (state.prefixIdx === 0 && state.stepIdx <= 0);
  }

  // Playback-bar prefix buttons
  el("prefixBackBtn").addEventListener("click", () => { stopPlaying(); goToPrefix(state.prefixIdx - 1); });
  el("prefixForwardBtn").addEventListener("click", () => { stopPlaying(); goToPrefix(state.prefixIdx + 1); });

  // Candidate navigation
  el("candidateBackBtn").addEventListener("click", () => {
    const step = currentStep();
    if (step && step.candidates.length > 0) {
      state.candidateIdx = (state.candidateIdx - 1 + step.candidates.length) % step.candidates.length;
      render();
    }
  });

  el("candidateForwardBtn").addEventListener("click", () => {
    const step = currentStep();
    if (step && step.candidates.length > 0) {
      state.candidateIdx = (state.candidateIdx + 1) % step.candidates.length;
      render();
    }
  });

  el("candidateInput").addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      const step = currentStep();
      if (!step || step.candidates.length === 0) return;
      const val = parseInt(el("candidateInput").value);
      if (!isNaN(val) && val >= 0 && val < step.candidates.length) {
        state.candidateIdx = val;
        render();
      }
    }
  });

  function advance() {
    const pfx = currentPrefix();
    if (!pfx) return;
    const step = currentStep();
    
    // First, try to cycle through candidates at current step
    if (step && step.candidates.length > 1) {
      const nextCandidateIdx = (state.candidateIdx + 1) % step.candidates.length;
      if (nextCandidateIdx !== 0) {
        // Still more candidates to show at this step
        state.candidateIdx = nextCandidateIdx;
        render();
        return;
      }
    }
    
    // All candidates shown, advance to next step
    if (state.stepIdx < pfx.steps.length - 1) {
      goToStep(state.stepIdx + 1);
    } else {
      // Reuse goToStep's all-dead-pause logic.
      goToStep(pfx.steps.length); // one past the end
    }
  }

  function startPlaying() {
    if (!state.trace || state.playing) return;
    state.playing = true;
    updatePlayButton();
    const baseDelay = 280; // Base delay in ms
    const tick = () => {
      advance();
      const delay = baseDelay / currentSpeedMultiplier;
      state.playTimer = setTimeout(tick, delay);
    };
    tick();
  }
  function stopPlaying() {
    state.playing = false;
    updatePlayButton();
    if (state.playTimer) clearTimeout(state.playTimer);
    state.playTimer = null;
  }
  playBtn.addEventListener("click", () => {
    if (state.playing) {
      stopPlaying();
    } else {
      // If at the end, restart from beginning
      const pfx = currentPrefix();
      if (pfx && state.prefixIdx === state.trace.prefixes.length - 1 && 
          state.stepIdx === pfx.steps.length - 1) {
        goToPrefix(0);
        goToStep(0);
      }
      startPlaying();
    }
  });

  el("stepBackBtn").addEventListener("click", () => {
    stopPlaying();
    goToStep(state.stepIdx - 1);
  });

  el("stepForwardBtn").addEventListener("click", () => {
    stopPlaying();
    goToStep(state.stepIdx + 1);
  });

  // Step input field - allow direct navigation
  stepInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      const val = stepInput.value.trim();
      const match = val.match(/^(\d+)/);
      if (match) {
        const targetStep = parseInt(match[1], 10) - 2; // Convert 1-based display to 0-based stepIdx (-1 is "before first")
        stopPlaying();
        goToStep(targetStep);
      }
    }
  });

  // Segment input field - allow direct navigation
  segmentInput.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      const val = segmentInput.value.trim();
      const match = val.match(/^(\d+)/);
      if (match) {
        const targetPrefix = parseInt(match[1], 10) - 1; // Convert 1-based display to 0-based prefixIdx
        stopPlaying();
        goToPrefix(targetPrefix);
      }
    }
  });

  // Speed preset buttons
  speedPresets.forEach((btn) => {
    btn.addEventListener("click", () => {
      const speed = parseFloat(btn.dataset.speed);
      currentSpeedMultiplier = speed;
      speedPresets.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      // If already playing, restart with new speed
      if (state.playing) {
        stopPlaying();
        startPlaying();
      }
    });
  });

  // -------------------------------------------------------------------------
  //  Fréchet distance
  // -------------------------------------------------------------------------

  function setFrechetResult(html) {
    frechetResult.innerHTML = html;
  }

  async function startFrechetComputation() {
    if (!state.trace) return;
    console.log('startFrechetComputation called - currentTraceId:', currentTraceId, 'currentFile:', currentFile);
    state.computingFrechet = true;
    renderParamsBar();

    try {
      let body;
      if (currentTraceId) {
        body = JSON.stringify({ trace_id: currentTraceId, eps: state.trace.eps, delta: state.trace.delta });
        console.log('Using trace_id for Fréchet:', currentTraceId);
      } else if (currentFile) {
        const text = await currentFile.text();
        body = JSON.stringify({ file_content: text, eps: state.trace.eps, delta: state.trace.delta });
        console.log('Using file_content for Fréchet');
      } else {
        console.log('No currentTraceId or currentFile - cannot start Fréchet computation');
        state.computingFrechet = false;
        renderParamsBar();
        return;
      }
      console.log('Posting to /api/frechet...');
      const resp = await fetch("/api/frechet", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body,
      });
      if (!resp.ok) {
        throw new Error(await apiErrorMessage(resp));
      }
      const data = await resp.json();
      state.computedFrechet = data.distance;
      state.computingFrechet = false;
      console.log('Fréchet computation finished:', data.distance);
      renderParamsBar();
    } catch (e) {
      console.error('Fréchet computation error:', e);
      state.computingFrechet = false;
      renderParamsBar();
    }
  }

  window.addEventListener("keydown", (e) => {
    console.log("Key pressed:", e.key, "Target:", e.target.tagName, "Has trace:", !!state.trace);
    if (!state.trace) return;
    // Ignore if an input, textarea, or select is focused
    const tag = e.target.tagName;
    if (tag === "INPUT" || tag === "TEXTAREA" || tag === "SELECT") {
      console.log("Ignoring - form element focused");
      return;
    }
    if (e.key === " ") { 
      console.log("Space pressed - toggling playback");
      e.preventDefault(); 
      if (state.playing) {
        stopPlaying();
      } else {
        // If at the end, restart from beginning
        const pfx = currentPrefix();
        if (pfx && state.prefixIdx === state.trace.prefixes.length - 1 && 
            state.stepIdx === pfx.steps.length - 1) {
          goToPrefix(0);
          goToStep(0);
        }
        startPlaying();
      }
      return;
    }
    if (e.key === "ArrowRight") {
      console.log("Arrow right pressed, shift:", e.shiftKey);
      e.preventDefault();
      if (e.shiftKey) goToPrefix(state.prefixIdx + 1);
      else goToStep(state.stepIdx + 1);
    } else if (e.key === "ArrowLeft") {
      console.log("Arrow left pressed, shift:", e.shiftKey);
      e.preventDefault();
      if (e.shiftKey) goToPrefix(state.prefixIdx - 1);
      else goToStep(state.stepIdx - 1);
    } else if (e.key === "c" || e.key === "C") {
      console.log("C pressed - cycling candidate");
      e.preventDefault();
      el("candidateForwardBtn").click();
    } else if (e.key === "x" || e.key === "X") {
      console.log("X pressed - cycling candidate in reverse");
      e.preventDefault();
      el("candidateBackBtn").click();
    }
  });

  for (const t of Object.values(toggles)) t.addEventListener("change", render);

  // -------------------------------------------------------------------------
  //  Drawing helpers
  // -------------------------------------------------------------------------

  function drawPolyline(pts, close) {
    if (!pts || pts.length < 2) return;
    ctx.beginPath();
    const [x0, y0] = worldToScreen(pts[0][0], pts[0][1]);
    ctx.moveTo(x0, y0);
    for (let i = 1; i < pts.length; ++i) {
      const [x, y] = worldToScreen(pts[i][0], pts[i][1]);
      ctx.lineTo(x, y);
    }
    if (close) ctx.closePath();
  }

  function fillPolygon(pts, fillStyle, strokeStyle, lineWidth) {
    if (!pts || pts.length < 3) return;
    drawPolyline(pts, true);
    if (fillStyle) { ctx.fillStyle = fillStyle; ctx.fill(); }
    if (strokeStyle) {
      ctx.strokeStyle = strokeStyle;
      ctx.lineWidth = lineWidth || 1.5;
      ctx.stroke();
    }
  }

  function strokePath(pts, strokeStyle, lineWidth) {
    drawPolyline(pts, false);
    ctx.strokeStyle = strokeStyle;
    ctx.lineWidth = lineWidth || 1.5;
    ctx.stroke();
  }

  // Ray-casting point-in-polygon (world coords). Returns true if pt is inside poly.
  function pointInPolygon(pt, poly) {
    if (!poly || poly.length < 3) return false;
    const [x, y] = pt;
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const [xi, yi] = poly[i], [xj, yj] = poly[j];
      if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi))
        inside = !inside;
    }
    return inside;
  }

  function dot(p, radius, fillStyle, strokeStyle) {
    const [x, y] = worldToScreen(p[0], p[1]);
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, Math.PI * 2);
    if (fillStyle) { ctx.fillStyle = fillStyle; ctx.fill(); }
    if (strokeStyle) { ctx.strokeStyle = strokeStyle; ctx.lineWidth = 1.5; ctx.stroke(); }
  }

  function cross(p, radius, strokeStyle) {
    const [x, y] = worldToScreen(p[0], p[1]);
    ctx.strokeStyle = strokeStyle;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(x - radius, y - radius); ctx.lineTo(x + radius, y + radius);
    ctx.moveTo(x - radius, y + radius); ctx.lineTo(x + radius, y - radius);
    ctx.stroke();
  }

  // -------------------------------------------------------------------------
  //  Axis ticks
  // -------------------------------------------------------------------------

  function drawAxisTicks() {
    const t = state.trace;
    if (!t || !t.bbox) return;

    const [xmin, ymin, xmax, ymax] = t.bbox;
    const xRange = xmax - xmin;
    const yRange = ymax - ymin;

    // Determine tick spacing based on data size
    const getTickSpacing = (range) => {
      const roughTicks = 8;
      const roughSpacing = range / roughTicks;
      const magnitude = Math.pow(10, Math.floor(Math.log10(roughSpacing)));
      const normalized = roughSpacing / magnitude;
      let nice;
      if (normalized < 1.5) nice = 1;
      else if (normalized < 3) nice = 2;
      else if (normalized < 7) nice = 5;
      else nice = 10;
      return nice * magnitude;
    };

    const xSpacing = getTickSpacing(xRange);
    const ySpacing = getTickSpacing(yRange);

    // Position axes at the data bounds instead of world origin (0,0)
    const xAxisY = ymin;
    const yAxisX = xmin;

    ctx.strokeStyle = "#444";
    ctx.fillStyle = "#999";
    ctx.lineWidth = 1;
    ctx.font = "11px monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";

    // X-axis ticks
    const xStart = Math.ceil(xmin / xSpacing) * xSpacing;
    for (let x = xStart; x <= xmax; x += xSpacing) {
      const [sx, sy] = worldToScreen(x, xAxisY);
      ctx.beginPath();
      ctx.moveTo(sx, sy - 5);
      ctx.lineTo(sx, sy + 5);
      ctx.stroke();
      ctx.fillText(x.toFixed(0), sx, sy + 8);
    }

    // Y-axis ticks
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    const yStart = Math.ceil(ymin / ySpacing) * ySpacing;
    for (let y = yStart; y <= ymax; y += ySpacing) {
      const [sx, sy] = worldToScreen(yAxisX, y);
      ctx.beginPath();
      ctx.moveTo(sx - 5, sy);
      ctx.lineTo(sx + 5, sy);
      ctx.stroke();
      ctx.fillText(y.toFixed(0), sx - 8, sy);
    }

    // Draw axes at the data bounds
    ctx.strokeStyle = "#666";
    ctx.lineWidth = 1.5;
    // X-axis at bottom of data
    const [x0Left, y0] = worldToScreen(xmin, xAxisY);
    const [x0Right, _] = worldToScreen(xmax, xAxisY);
    ctx.beginPath();
    ctx.moveTo(x0Left, y0);
    ctx.lineTo(x0Right, y0);
    ctx.stroke();
    // Y-axis at left of data
    const [x0, y0Bottom] = worldToScreen(yAxisX, ymin);
    const [__, y0Top] = worldToScreen(yAxisX, ymax);
    ctx.beginPath();
    ctx.moveTo(x0, y0Bottom);
    ctx.lineTo(x0, y0Top);
    ctx.stroke();

    // Draw bounding box (always visible, solid lines)
    ctx.strokeStyle = "#888";
    ctx.lineWidth = 1.5;
    const [bx1, by1] = worldToScreen(xmin, ymin);
    const [bx2, by2] = worldToScreen(xmax, ymax);
    ctx.strokeRect(bx1, by2, bx2 - bx1, by1 - by2);
  }

  // -------------------------------------------------------------------------
  //  Main render
  // -------------------------------------------------------------------------

  function render() {
    resizeCanvas();
    const rect = canvasWrap.getBoundingClientRect();
    ctx.clearRect(0, 0, rect.width, rect.height);
    renderStatus();

    const t = state.trace;
    if (!t) {
      if (state.preview && state.preview.points.length) {
        strokePath(state.preview.points, "#4f9dff", 1.6);
        const stride = Math.max(1, Math.ceil(state.preview.points.length / 1200));
        for (let i = 0; i < state.preview.points.length; i += stride) {
          dot(state.preview.points[i], 1.25, "#7bb6ff", null);
        }
      }
      return;
    }

    const pfx = currentPrefix();
    const step = currentStep();

    // Draw axis ticks
    drawAxisTicks();

    // Update navigation input fields
    segmentInput.value = `${state.prefixIdx + 1} / ${t.prefixes.length}`;
    const n = pfx ? pfx.steps.length : 0;
    stepInput.value = `${state.stepIdx + 2} / ${n + 1}`;
    
    // For candidate display, show position in the cycle pool (alive + justDied candidates)
    if (step && step.candidates.length > 0) {
      const statusOf = (c) =>
        c.alive ? 'alive' : (c.F && c.F.length >= 3 ? 'justDied' : 'dead');
      const allCandidates = step.candidates.map((c, i) => ({ ...c, originalIdx: i }));
      const cyclePool = allCandidates.filter(c => statusOf(c) !== 'dead');
      const displayIdx = cyclePool.length > 0 ? (state.candidateIdx % cyclePool.length) + 1 : 0;
      candidateInput.value = `${displayIdx} / ${cyclePool.length}`;
    } else {
      candidateInput.value = `0 / 0`;
    }
    if (mobileStepBackBtn) {
      mobileStepBackBtn.disabled = state.prefixIdx === 0 && state.stepIdx <= 0;
    }
    if (mobileStepForwardBtn) {
      const atLastPrefix = state.prefixIdx === t.prefixes.length - 1;
      mobileStepForwardBtn.disabled = atLastPrefix && state.stepIdx >= n - 1;
    }
    const hasCandidates = Boolean(step && step.candidates.length);
    if (mobileCandidateBackBtn) mobileCandidateBackBtn.disabled = !hasCandidates;
    if (mobileCandidateForwardBtn) mobileCandidateForwardBtn.disabled = !hasCandidates;

    // 1. Full input stream (faint polyline + small filled dots).
    if (toggles.stream.checked) {
      strokePath(t.stream, "#3d4456", 1.25);
      
      // Collect all simplified output points to avoid drawing them twice
      const simplifiedPoints = new Set();
      for (let i = 0; i < state.prefixIdx; ++i) {
        const pfxOutput = t.prefixes[i].output;
        simplifiedPoints.add(`${pfxOutput[0][0]},${pfxOutput[0][1]}`);
        simplifiedPoints.add(`${pfxOutput[1][0]},${pfxOutput[1][1]}`);
      }
      if (pfx && step && state.stepIdx === pfx.steps.length - 1) {
        simplifiedPoints.add(`${pfx.output[0][0]},${pfx.output[0][1]}`);
        // Do NOT add output[1] if this is the very last prefix - that's the trajectory endpoint
        const isLastPrefix = state.prefixIdx === t.prefixes.length - 1;
        if (!isLastPrefix) {
          simplifiedPoints.add(`${pfx.output[1][0]},${pfx.output[1][1]}`);
        }
      }
      
      // Draw stream dots, but skip those that are part of simplified output
      for (const p of t.stream) {
        const key = `${p[0]},${p[1]}`;
        if (!simplifiedPoints.has(key)) {
          dot(p, 0.8, "#4a5266", null);
        }
      }
    }

    // 3. Simplified output committed so far (all prior prefixes' outputs).
    if (toggles.simplified.checked) {
      const committed = [];
      for (let i = 0; i < state.prefixIdx; ++i) {
        committed.push(t.prefixes[i].output[0], t.prefixes[i].output[1]);
      }
      // Do NOT push pfx.output[0] during mid-prefix steps: it coincides with a
      // P-grid corner, so a green dot would foreshadow which candidate "wins".
      // Exception: once we are on the last step of this prefix the computation is
      // complete, so we can safely reveal the output segment (this also ensures the
      // final prefix's segment is shown — it would otherwise never be promoted into
      // a "prior" prefix by advancing to a non-existent next prefix).
      // BUT: do NOT show output if we're at stepIdx === -1 (before first step of prefix)
      if (pfx && step && state.stepIdx >= 0 && state.stepIdx === pfx.steps.length - 1) {
        committed.push(pfx.output[0], pfx.output[1]);
      }
      // Only draw if we have at least 2 points (a complete segment)
      if (committed.length >= 2) {
        strokePath(committed, "#3ddc97", 2.25);
        // Show dots for all committed points except the trajectory endpoint
        // If at the final step of the final prefix, exclude the last point (trajectory endpoint)
        const isAtFinalStepOfFinalPrefix = pfx && step && 
                                            state.stepIdx === pfx.steps.length - 1 && 
                                            state.prefixIdx === t.prefixes.length - 1;
        const numToExclude = isAtFinalStepOfFinalPrefix ? 1 : 0;
        const dotsToShow = committed.length > numToExclude ? committed.slice(0, committed.length - numToExclude) : [];
        for (const p of dotsToShow) dot(p, 1.8, "#3ddc97", null);
      }
    }

    if (pfx) {
      // Determine "best" candidate index for this step (the one the buffer
      // was pulled from — mirrors the reverse scan in get_longest_stab).
      let bestIdx = -1;
      if (step) {
        for (let i = step.candidates.length - 1; i >= 0; --i) {
          const c = step.candidates[i];
          if (!c.alive || !c.S || c.S.length === 0) continue;
          bestIdx = i;
          break;
        }
      }

      // 4. δ-balls — each independently toggleable.
      {
        const r = t.r_val;
        const drawBall = (cx, cy, fillColor, strokeColor, lw) => {
          const [scx, scy] = worldToScreen(cx, cy);
          const rScreen = r * state.cam.scale;
          ctx.beginPath();
          ctx.arc(scx, scy, rScreen, 0, 2 * Math.PI);
          if (fillColor) { ctx.fillStyle = fillColor; ctx.fill(); }
          ctx.strokeStyle = strokeColor;
          ctx.lineWidth = lw;
          ctx.stroke();
        };
        // Ball at p0: the ball whose grid corners define the P anchor set.
        if (toggles["ball-p0"].checked) {
          drawBall(pfx.p0[0], pfx.p0[1], "rgba(196,122,255,0.07)", "rgba(196,122,255,0.6)", 1.5);
        }
        // Ball at v_i: the ball that conv(G_i) must cover at this step.
        if (step && step.pi && toggles["ball-pi"].checked) {
          drawBall(step.pi[0], step.pi[1], "rgba(255,122,232,0.07)", "rgba(255,122,232,0.85)", 1.75);
        }
      }

      // 4b. G_i delta-disk convex hull for the point being consumed this step.
      if (step && toggles.Gi.checked) {
        fillPolygon(step.Gi, null, "#e8c547", 1.75);
      }

      // 5/6/7/8/10. Candidate display — single shared cycle pool.
      // Cycle pool = alive + just-died (dead but have F data from this step).
      // Previously-dead (no F this step) are dimmed but not cycled through.
      if (step) {
        const allCandidates = step.candidates
          .map((c, i) => ({ ...c, originalIdx: i }));

        // 'alive' | 'justDied' (dead but has F this step) | 'dead' (previously dead)
        const statusOf = c =>
          c.alive ? 'alive' : (c.F && c.F.length >= 3 ? 'justDied' : 'dead');

        const COLORS = {
          alive:    { fFill: "rgba(79,157,255,0.18)",  fStroke: "rgba(79,157,255,0.85)",  ray: "rgba(79,157,255,1)",   lw: 2,    dotR: 3,   dotFill: "#ff9f43", dotStroke: "#ff6b1a" },
          justDied: { fFill: "rgba(255,60,60,0.14)",    fStroke: "rgba(255,60,60,0.9)",    ray: "rgba(255,60,60,1)",    lw: 1.75, dotR: 3,   dotFill: "#ff3c3c", dotStroke: "#cc1a1a" },
          dead:     { fFill: "rgba(255,60,60,0.08)",  fStroke: "rgba(255,60,60,0.7)",   ray: "rgba(255,60,60,0.7)", lw: 1.5,  dotR: 2.5, dotFill: "#ff3c3c", dotStroke: "#cc1a1a" },
        };

        const cyclePool = allCandidates.filter(c => statusOf(c) !== 'dead');
        const viewedCandidate = cyclePool.length > 0
          ? cyclePool[state.candidateIdx % cyclePool.length]
          : null;

        // 5/6. F wedge, tangent rays, S region, anchor highlight.
        if (viewedCandidate) {
          const i = viewedCandidate.originalIdx;
          const c = viewedCandidate;
          const status = statusOf(c);
          const col = COLORS[status];
          const anchor = pfx.P[i];

          if (toggles.F.checked && c.F && c.F.length >= 3) {
            const isBbox = (() => {
              if (c.F.length !== 4 || !state.trace) return false;
              const [xmin, ymin, xmax, ymax] = state.trace.bbox;
              const corners = [[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]];
              return c.F.every(fp =>
                corners.some(co =>
                  Math.abs(fp[0]-co[0]) < 1 && Math.abs(fp[1]-co[1]) < 1));
            })();

            if (status === 'alive' && isBbox) {
              ctx.save();
              ctx.strokeStyle = "rgba(79,157,255,0.7)";
              ctx.lineWidth = 1.5;
              ctx.setLineDash([6, 4]);
              fillPolygon(c.F, null, "rgba(79,157,255,0.7)", 1.5);
              ctx.setLineDash([]);
              ctx.restore();
            } else {
              fillPolygon(c.F,
                isBbox ? null : col.fFill,
                isBbox ? col.fStroke.replace(/[\d.]+\)$/, "0.45)") : col.fStroke,
                col.lw);
            }

            // Tangent rays from anchor → the two extremal vertices of F.
            if (anchor && c.F.length >= 1 && !pointInPolygon(anchor, c.F)) {
              const [ax, ay] = [anchor[0], anchor[1]];
              const angles = c.F.map(sp => Math.atan2(sp[1] - ay, sp[0] - ax));
              const indexed = angles.map((a, ii) => ({ a, ii }))
                .sort((u, v) => u.a - v.a);
              let maxGap = -1, gapAfter = 0;
              for (let k = 0; k < indexed.length; k++) {
                const next = (k + 1) % indexed.length;
                let gap = indexed[next].a - indexed[k].a;
                if (gap < 0) gap += 2 * Math.PI;
                if (gap > maxGap) { maxGap = gap; gapAfter = k; }
              }
              const t1 = c.F[indexed[gapAfter].ii];
              const t2 = c.F[indexed[(gapAfter + 1) % indexed.length].ii];
              ctx.save();
              ctx.strokeStyle = col.ray;
              ctx.lineWidth = col.lw;
              ctx.setLineDash([5, 5]);
              const [asx, asy] = worldToScreen(ax, ay);
              for (const tp of [t1, t2]) {
                const [tx, ty] = worldToScreen(tp[0], tp[1]);
                ctx.beginPath(); ctx.moveTo(asx, asy); ctx.lineTo(tx, ty); ctx.stroke();
              }
              ctx.setLineDash([]);
              ctx.restore();
            }
          }

          // F_Si wedge for the viewed candidate
          if (toggles["F-Si"].checked && c.F_Si && c.F_Si.length >= 3) {
            const isBbox = (() => {
              if (c.F_Si.length !== 4 || !state.trace) return false;
              const [xmin, ymin, xmax, ymax] = state.trace.bbox;
              const corners = [[xmin,ymin],[xmax,ymin],[xmax,ymax],[xmin,ymax]];
              return c.F_Si.every(fp =>
                corners.some(co =>
                  Math.abs(fp[0]-co[0]) < 1 && Math.abs(fp[1]-co[1]) < 1));
            })();

            // Use cyan color for F_Si
            const fSiCol = {
              fFill: "rgba(34,211,238,0.12)",
              fStroke: "rgba(34,211,238,0.85)",
              ray: "rgba(34,211,238,0.6)",
              lw: 1.5
            };

            if (status === 'alive' && isBbox) {
              ctx.save();
              ctx.strokeStyle = "rgba(34,211,238,0.7)";
              ctx.lineWidth = 1.5;
              ctx.setLineDash([6, 4]);
              fillPolygon(c.F_Si, null, "rgba(34,211,238,0.7)", 1.5);
              ctx.setLineDash([]);
              ctx.restore();
            } else {
              fillPolygon(c.F_Si,
                isBbox ? null : fSiCol.fFill,
                isBbox ? fSiCol.fStroke.replace(/[\d.]+\)$/, "0.45)") : fSiCol.fStroke,
                fSiCol.lw);
            }

            // Tangent rays from anchor → the two extremal vertices of F_Si.
            if (anchor && c.F_Si.length >= 1 && !pointInPolygon(anchor, c.F_Si)) {
              const [ax, ay] = [anchor[0], anchor[1]];
              const angles = c.F_Si.map(sp => Math.atan2(sp[1] - ay, sp[0] - ax));
              const indexed = angles.map((a, ii) => ({ a, ii }))
                .sort((u, v) => u.a - v.a);
              let maxGap = -1, gapAfter = 0;
              for (let k = 0; k < indexed.length; k++) {
                const next = (k + 1) % indexed.length;
                let gap = indexed[next].a - indexed[k].a;
                if (gap < 0) gap += 2 * Math.PI;
                if (gap > maxGap) { maxGap = gap; gapAfter = k; }
              }
              const t1 = c.F_Si[indexed[gapAfter].ii];
              const t2 = c.F_Si[indexed[(gapAfter + 1) % indexed.length].ii];
              ctx.save();
              ctx.strokeStyle = fSiCol.ray;
              ctx.lineWidth = fSiCol.lw;
              ctx.setLineDash([5, 5]);
              const [asx, asy] = worldToScreen(ax, ay);
              for (const tp of [t1, t2]) {
                const [tx, ty] = worldToScreen(tp[0], tp[1]);
                ctx.beginPath(); ctx.moveTo(asx, asy); ctx.lineTo(tx, ty); ctx.stroke();
              }
              ctx.setLineDash([]);
              ctx.restore();
            }
          }

          // S region — only when alive.
          if (status === 'alive' && toggles.S.checked && c.S && c.S.length >= 3) {
            fillPolygon(c.S, "rgba(167,139,250,0.18)", "#a78bfa", 2);
          }

          if (toggles.P.checked && anchor) {
            dot(anchor, col.dotR, col.dotFill, col.dotStroke);
          }
        }

        // 7. Dead candidate anchors — red dots (only when toggle is checked).
        if (toggles["dead-candidates"].checked) {
          for (const c of allCandidates) {
            if (c.alive) continue;
            const anchor = pfx.P[c.originalIdx];
            if (anchor) {
              const [ax, ay] = worldToScreen(anchor[0], anchor[1]);
              ctx.fillStyle = "#888888";
              ctx.beginPath();
              ctx.arc(ax, ay, 1.2, 0, 2 * Math.PI);
              ctx.fill();
            }
          }
        }

        // 8. P candidate anchors (non-viewed).
        if (toggles.P.checked) {
          const viewedOrigIdx = viewedCandidate ? viewedCandidate.originalIdx : -1;
          const showDead = toggles["dead-candidates"].checked;
          for (let ii = 0; ii < pfx.P.length; ++ii) {
            if (ii === viewedOrigIdx) continue;
            const sc = step.candidates[ii];
            const st = sc ? statusOf(sc) : 'alive';
            
            // Skip dead/justDied candidates if toggle is off
            if (!showDead && (st === 'justDied' || st === 'dead')) continue;
            
            dot(pfx.P[ii], 1.2,
              st === 'alive' ? "#ff9f43" : st === 'justDied' ? "#ff3c3c" : "#7a4a1f",
              null);
          }
        }

        // 10. "N/M" label anchored to top-left of the entire P grid (doesn't move when cycling).
        if (viewedCandidate) {
          const currentCycleIdx = state.candidateIdx % cyclePool.length;
          const status = statusOf(viewedCandidate);
          
          // Find the top-left corner of all P anchors
          let minX = Infinity, minY = Infinity;
          for (const pt of pfx.P) {
            const [sx, sy] = worldToScreen(pt[0], pt[1]);
            if (sx < minX) minX = sx;
            if (sy < minY) minY = sy;
          }
          
          const headerText = `${currentCycleIdx + 1}/${cyclePool.length}`;
          ctx.save();
          ctx.font = "bold 12px -apple-system, BlinkMacSystemFont, sans-serif";
          const tw = ctx.measureText(headerText).width;
          // Position very close to the top-left corner of the orange grid
          const hx = minX - tw + 30;  // 8px left of leftmost dot
          const hy = minY;           // Aligned with topmost dot
          ctx.fillStyle = status === 'dead' ? "rgba(255,60,60,0.85)"
            : status === 'justDied' ? "rgba(255,60,60,0.95)"
            : "#ff9f43";
          ctx.fillText(headerText, hx, hy);
          ctx.restore();
        }
      }

      // 9. Small canvas labels — v_i with actual stream index in pink, and p in orange.
      if (state.currentStartPoint) {
        // Add orange "p" label for the active anchor point
        const [pax, pay] = worldToScreen(state.currentStartPoint[0], state.currentStartPoint[1]);
        ctx.save();
        ctx.font = "bold 11px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.fillStyle = "#ff9f43";
        ctx.fillText(`p`, pax + 7, pay - 6);
        ctx.restore();
      }
      if (step) {
        const [pix, piy] = worldToScreen(step.pi[0], step.pi[1]);
        ctx.save();
        ctx.font = "bold 11px -apple-system, BlinkMacSystemFont, sans-serif";
        ctx.fillStyle = "#ff7ae8";
        ctx.fillText(`v${step.stream_idx}`, pix + 7, piy - 6);
        ctx.restore();
      }
    }
  }

  // -------------------------------------------------------------------------
  //  Resizer for adjustable sidebar
  // -------------------------------------------------------------------------

  const resizer = el("resizer");
  const sidebar = el("sidebar");
  let isResizing = false;

  resizer.addEventListener("mousedown", (e) => {
    isResizing = true;
    e.preventDefault();
  });

  window.addEventListener("mousemove", (e) => {
    if (!isResizing) return;
    const containerWidth = document.getElementById("mainContent").offsetWidth;
    const newWidth = containerWidth - e.clientX;
    if (newWidth >= 250 && newWidth <= 600) {
      sidebar.style.width = newWidth + "px";
      resizer.style.right = newWidth + "px";
      resizeCanvas();
      render();
    }
  });

  window.addEventListener("mouseup", () => {
    isResizing = false;
  });

  // -------------------------------------------------------------------------
  //  Boot
  // -------------------------------------------------------------------------

  resizeCanvas();
  render();
})();
