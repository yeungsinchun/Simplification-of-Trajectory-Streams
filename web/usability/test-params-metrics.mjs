#!/usr/bin/env node
/**
 * Focused browser test: params-bar spinner placement for
 * Computed Fréchet distance and Simplification time.
 *
 * Drives the real viewer (index.html + viewer.js) against a mock NDJSON
 * stream so loading and loaded states are observable. Writes screenshots
 * to EVIDENCE_DIR (required for evidence collection).
 *
 *   EVIDENCE_DIR=/path/to/evidence node web/usability/test-params-metrics.mjs
 */
import http from "node:http";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";
import { spawnSync } from "node:child_process";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const WEB_DIR = path.resolve(__dirname, "..");
const EVIDENCE_DIR = process.env.EVIDENCE_DIR;
const CHROME =
  process.env.CHROME_PATH ||
  "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";

if (!EVIDENCE_DIR) {
  console.error("EVIDENCE_DIR is required");
  process.exit(2);
}
fs.mkdirSync(EVIDENCE_DIR, { recursive: true });

function loadPuppeteer() {
  const tmp = path.join(os.tmpdir(), "params-metrics-puppeteer");
  fs.mkdirSync(tmp, { recursive: true });
  const pkg = path.join(tmp, "package.json");
  if (!fs.existsSync(pkg)) {
    fs.writeFileSync(pkg, '{"private":true}\n');
  }
  const requireFromTmp = createRequire(path.join(tmp, "package.json"));
  try {
    return requireFromTmp("puppeteer-core");
  } catch {
    const install = spawnSync(
      "npm",
      ["install", "--no-save", "--no-fund", "--no-audit", "puppeteer-core@24"],
      { cwd: tmp, stdio: "inherit", env: process.env },
    );
    if (install.status !== 0) {
      throw new Error("npm install puppeteer-core failed");
    }
    return requireFromTmp("puppeteer-core");
  }
}

const MIME = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".css": "text/css; charset=utf-8",
  ".png": "image/png",
  ".md": "text/markdown; charset=utf-8",
};

function json(res, status, body) {
  res.writeHead(status, {
    "Content-Type": "application/json",
    "Cache-Control": "no-store",
  });
  res.end(JSON.stringify(body));
}

function serveStatic(req, res) {
  const url = new URL(req.url, "http://127.0.0.1");
  let rel = decodeURIComponent(url.pathname);
  if (rel === "/") rel = "/index.html";
  const filePath = path.normalize(path.join(WEB_DIR, rel));
  if (!filePath.startsWith(WEB_DIR)) {
    res.writeHead(403);
    res.end("forbidden");
    return;
  }
  fs.readFile(filePath, (err, data) => {
    if (err) {
      res.writeHead(404);
      res.end("not found");
      return;
    }
    res.writeHead(200, { "Content-Type": MIME[path.extname(filePath)] || "application/octet-stream" });
    res.end(data);
  });
}

function createMockServer() {
  const streamHolders = [];
  const frechetHolders = [];

  const header = {
    type: "header",
    eps: 0.9,
    delta: 500,
    grid_val: 1.25,
    r_val: 0.45,
    expected_frechet: 1.8,
    bbox: [0, 0, 10, 10],
    stream: [[0, 0], [4, 3], [10, 0]],
  };
  const prefix = {
    type: "prefix",
    data: {
      end_idx: 2,
      p0: [0, 0],
      p0_idx: 0,
      output: [[0, 0], [10, 0]],
      P: [[0, 0]],
      steps: [{ stream_idx: 1, pi: [4, 3], candidates: [] }],
    },
  };
  const done = {
    type: "done",
    time_ms: 42.18,
    simplified: [[0, 0], [10, 0]],
    frechet_distance: 0.5,
  };

  const server = http.createServer((req, res) => {
    const url = new URL(req.url, "http://127.0.0.1");

    if (url.pathname === "/api/traces") {
      json(res, 200, {
        traces: [{ id: 1, n_points: 3, label: "Tiny mock trace", epsilon: 0.9, delta: 500 }],
      });
      return;
    }

    if (url.pathname === "/api/trace/1" && req.method === "GET") {
      res.writeHead(200, {
        "Content-Type": "application/x-ndjson",
        "Cache-Control": "no-store",
        "X-Accel-Buffering": "no",
      });
      res.write(JSON.stringify(header) + "\n");
      const release = () => {
        if (res.writableEnded) return;
        res.write(JSON.stringify(prefix) + "\n");
        res.write(JSON.stringify(done) + "\n");
        res.end();
      };
      streamHolders.push({ release });
      return;
    }

    if (url.pathname === "/api/frechet" && req.method === "POST") {
      let body = "";
      req.on("data", (chunk) => { body += chunk; });
      req.on("end", () => {
        const answer = () => {
          if (res.writableEnded) return;
          json(res, 200, { distance: 0.1842 });
        };
        frechetHolders.push({ answer });
      });
      return;
    }

    serveStatic(req, res);
  });

  return {
    server,
    releaseStream() {
      for (const h of streamHolders.splice(0)) h.release();
    },
    answerFrechet() {
      for (const h of frechetHolders.splice(0)) h.answer();
    },
    listen() {
      return new Promise((resolve) => {
        server.listen(0, "127.0.0.1", () => resolve(server.address().port));
      });
    },
    close() {
      return new Promise((resolve) => server.close(resolve));
    },
  };
}

function assert(cond, message) {
  if (!cond) throw new Error(message);
}

async function metricSnapshot(page, kind) {
  return page.$eval(`.params-metric--${kind}`, (el) => {
    const body = el.querySelector(".params-metric-body");
    const label = el.querySelector(".params-metric-label");
    const value = el.querySelector(".params-metric-value");
    const spinner = el.querySelector(".params-spinner");
    const children = [...body.children].map((c) =>
      [...c.classList].filter((n) => n.startsWith("params-metric-")).join(" ")
    );
    const labelRect = label.getBoundingClientRect();
    const valueRect = value.getBoundingClientRect();
    const spinnerRect = spinner ? spinner.getBoundingClientRect() : null;
    return {
      label: label.textContent.trim(),
      valueText: (value.textContent || "").replace(/\s+/g, " ").trim(),
      loading: value.classList.contains("params-metric-value--loading"),
      hasSpinner: Boolean(spinner),
      childOrder: children,
      labelRight: labelRect.right,
      valueLeft: valueRect.left,
      spinnerLeft: spinnerRect ? spinnerRect.left : null,
      visible: el.getClientRects().length > 0 && getComputedStyle(el).display !== "none",
    };
  });
}

function checkSpinnerInValueSlot(snap, kind) {
  assert(snap.childOrder[0] === "params-metric-label", `${kind}: label is not the first child`);
  assert(
    snap.childOrder.some((c) => c.includes("params-metric-value")),
    `${kind}: missing value slot`,
  );
  assert(snap.hasSpinner, `${kind}: expected spinner while loading`);
  assert(snap.loading, `${kind}: value slot should be in loading state`);
  assert(
    snap.spinnerLeft != null && snap.spinnerLeft >= snap.labelRight - 1,
    `${kind}: spinner is not after the label (spinnerLeft=${snap.spinnerLeft}, labelRight=${snap.labelRight})`,
  );
  assert(
    snap.valueLeft >= snap.labelRight - 1,
    `${kind}: value slot is not after the label`,
  );
}

function checkLoadedValue(snap, kind, expectedSubstring) {
  assert(!snap.hasSpinner, `${kind}: spinner still present after load`);
  assert(!snap.loading, `${kind}: still marked loading`);
  assert(
    snap.valueText.includes(expectedSubstring),
    `${kind}: expected value containing ${expectedSubstring}, got ${JSON.stringify(snap.valueText)}`,
  );
  assert(
    snap.valueLeft >= snap.labelRight - 1,
    `${kind}: loaded value is not after the label`,
  );
}

async function screenshot(page, name, clipSelector) {
  const dest = path.join(EVIDENCE_DIR, name);
  if (clipSelector) {
    const handle = await page.$(clipSelector);
    if (handle) {
      const box = await handle.boundingBox();
      if (box && box.width > 0 && box.height > 0) {
        const pad = 8;
        await page.screenshot({
          path: dest,
          clip: {
            x: Math.max(0, box.x - pad),
            y: Math.max(0, box.y - pad),
            width: box.width + pad * 2,
            height: box.height + pad * 2,
          },
        });
        return dest;
      }
    }
  }
  await page.screenshot({ path: dest });
  return dest;
}

async function loadTraceAndWaitForHeader(page) {
  await page.waitForFunction(
    () => document.querySelector("#traceSelect")?.options.length > 1,
    { timeout: 15000 },
  );
  await page.select("#traceSelect", "1");
  await page.click("#loadBtn");
  await page.waitForSelector(".params-metric--frechet", { timeout: 15000 });
  await page.waitForSelector(".params-metric--time", { timeout: 15000 });
}

async function run() {
  const puppeteer = loadPuppeteer();
  const mock = createMockServer();
  const port = await mock.listen();
  const origin = `http://127.0.0.1:${port}`;
  const artifacts = [];
  const notes = [];

  const browser = await puppeteer.launch({
    executablePath: CHROME,
    headless: "new",
    args: ["--hide-scrollbars", "--disable-gpu", "--no-sandbox"],
  });

  try {
    // --- Fixture: documented spinner placement at desktop and mobile ---
    const fixture = await browser.newPage();
    await fixture.setViewport({ width: 1280, height: 900, deviceScaleFactor: 2 });
    await fixture.goto(`${origin}/usability/params-metric-states.html`, {
      waitUntil: "domcontentloaded",
    });
    artifacts.push(await screenshot(fixture, "fixture-desktop-params-metrics.png"));
    await fixture.setViewport({ width: 390, height: 1200, deviceScaleFactor: 2 });
    await fixture.reload({ waitUntil: "domcontentloaded" });
    artifacts.push(await screenshot(fixture, "fixture-mobile-params-metrics.png"));
    await fixture.close();

    // --- Live desktop: header-only stream holds both chips with in-slot spinners ---
    const desktop = await browser.newPage();
    desktop.on("pageerror", (err) => console.error("desktop pageerror:", err));
    desktop.on("console", (msg) => {
      if (msg.type() === "error") console.error("desktop console:", msg.text());
    });
    await desktop.setViewport({ width: 1280, height: 800, deviceScaleFactor: 2 });
    await desktop.goto(origin, { waitUntil: "domcontentloaded" });

    const typo = await desktop.$eval(".desktop-instructions h2", (el) => el.textContent.trim());
    notes.push({ stillOpen: "desktop-intruction-typo", text: typo });
    artifacts.push(await screenshot(desktop, "desktop-drop-hint-intruction.png", "#dropHint"));
    assert(typo === "Intruction", `expected documented typo Intruction, got ${JSON.stringify(typo)}`);

    await loadTraceAndWaitForHeader(desktop);
    const loadingFrechet = await metricSnapshot(desktop, "frechet");
    const loadingTime = await metricSnapshot(desktop, "time");
    checkSpinnerInValueSlot(loadingFrechet, "frechet");
    checkSpinnerInValueSlot(loadingTime, "time");
    assert(loadingFrechet.visible && loadingTime.visible, "desktop loading: both metrics must be visible");
    const bothDuringLoad = await desktop.$$eval(".params-metric--blue", (els) => els.length);
    assert(bothDuringLoad === 2, `desktop loading: expected 2 blue metrics, got ${bothDuringLoad}`);
    artifacts.push(await screenshot(desktop, "desktop-live-loading-params.png", "#paramsBar"));
    artifacts.push(await screenshot(desktop, "desktop-live-loading-header.png", "header"));

    const timeLeftWhileFrechetSpinning = loadingTime.valueLeft;

    mock.releaseStream();
    try {
      await desktop.waitForFunction(() => {
        const time = document.querySelector(".params-metric--time .params-metric-value");
        return time && /ms/.test(time.textContent || "") && !time.querySelector(".params-spinner");
      }, { timeout: 10000 });
    } catch (err) {
      const html = await desktop.$eval("#paramsBar", (el) => el.innerHTML).catch(() => "<missing>");
      console.error("paramsBar after releaseStream:", html);
      throw err;
    }

    const midFrechet = await metricSnapshot(desktop, "frechet");
    const midTime = await metricSnapshot(desktop, "time");
    checkSpinnerInValueSlot(midFrechet, "frechet");
    checkLoadedValue(midTime, "time", "ms");
    artifacts.push(await screenshot(desktop, "desktop-live-time-loaded-frechet-spinning.png", "#paramsBar"));

    mock.answerFrechet();
    await desktop.waitForFunction(() => {
      const frechet = document.querySelector(".params-metric--frechet .params-metric-value");
      return frechet && !frechet.querySelector(".params-spinner") && /[0-9]/.test(frechet.textContent || "");
    }, { timeout: 10000 });

    const loadedFrechet = await metricSnapshot(desktop, "frechet");
    const loadedTime = await metricSnapshot(desktop, "time");
    checkLoadedValue(loadedFrechet, "frechet", "0.1842");
    checkLoadedValue(loadedTime, "time", "ms");
    artifacts.push(await screenshot(desktop, "desktop-live-loaded-params.png", "#paramsBar"));
    artifacts.push(await screenshot(desktop, "desktop-live-loaded-header.png", "header"));

    const shift = Math.abs(loadedTime.valueLeft - timeLeftWhileFrechetSpinning);
    notes.push({ timeValueSlotShiftPx: shift });
    await desktop.close();

    // --- Live mobile: spinner markup exists, but headerBody is hidden after load ---
    const mobile = await browser.newPage();
    await mobile.setViewport({ width: 390, height: 844, deviceScaleFactor: 2, isMobile: true, hasTouch: true });
    await mobile.goto(origin, { waitUntil: "domcontentloaded" });
    artifacts.push(await screenshot(mobile, "mobile-live-start.png"));

    await loadTraceAndWaitForHeader(mobile);
    const headerBodyDisplay = await mobile.$eval("#headerBody", (el) => getComputedStyle(el).display);
    const paramsVisible = await mobile.$eval("#paramsBar", (el) => {
      const rects = el.getClientRects();
      return rects.length > 0 && getComputedStyle(el).display !== "none" && getComputedStyle(el).visibility !== "hidden";
    });
    const mobileFrechet = await metricSnapshot(mobile, "frechet");
    const mobileTime = await metricSnapshot(mobile, "time");
    checkSpinnerInValueSlot(mobileFrechet, "frechet");
    checkSpinnerInValueSlot(mobileTime, "time");
    notes.push({
      stillOpen: "mobile-params-bar-hidden-after-load",
      headerBodyDisplay,
      paramsVisible,
    });
    assert(
      headerBodyDisplay === "none" && paramsVisible === false,
      "mobile still-open finding should still hold: headerBody hidden after load",
    );
    artifacts.push(await screenshot(mobile, "mobile-live-loading-header-hidden.png"));

    mock.releaseStream();
    mock.answerFrechet();
    await mobile.waitForFunction(() => document.body.classList.contains("trace-loaded-mobile"), { timeout: 10000 });
    artifacts.push(await screenshot(mobile, "mobile-live-after-load.png"));
    await mobile.close();

    const resultsPath = path.join(EVIDENCE_DIR, "params-metrics-test-results.json");
    fs.writeFileSync(
      resultsPath,
      JSON.stringify(
        {
          ok: true,
          loadingFrechet,
          loadingTime,
          loadedFrechet,
          loadedTime,
          notes,
          artifacts,
        },
        null,
        2,
      ),
    );
    artifacts.push(resultsPath);
    console.log(JSON.stringify({ ok: true, artifacts, notes }, null, 2));
  } finally {
    await browser.close();
    await mock.close();
  }
}

run().catch((err) => {
  console.error(err);
  process.exit(1);
});
