import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const algorithms = [
  ["ALPS", "#0072B2"],
  ["AWE", "#D55E00"],
  ["GAWE", "#009E73"],
  ["MGAWE", "#CC79A7"],
  ["WCAWE", "#E69F00"],
];

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; ++i) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"') {
        if (text[i + 1] === '"') {
          field += '"';
          ++i;
        } else {
          quoted = false;
        }
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(field);
      field = "";
    } else if (ch === "\n") {
      row.push(field);
      rows.push(row);
      row = [];
      field = "";
    } else if (ch !== "\r") {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field);
    rows.push(row);
  }
  const header = rows.shift().map((v) => v.trim().replace(/^"|"$/g, ""));
  return rows
    .filter((r) => r.some((v) => v.trim() !== ""))
    .map((r) => Object.fromEntries(header.map((h, i) => [h, (r[i] ?? "").trim()])));
}

function readCsv(relPath) {
  return parseCsv(fs.readFileSync(path.join(root, relPath), "utf8"));
}

function magToDb(mag) {
  return 20.0 * Math.log10(Math.max(Number(mag), 1e-300));
}

function interp(xs, ys, x) {
  if (x <= xs[0]) return ys[0];
  if (x >= xs[xs.length - 1]) return ys[ys.length - 1];
  let lo = 0;
  let hi = xs.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid;
    else hi = mid;
  }
  const t = (x - xs[lo]) / (xs[hi] - xs[lo]);
  return ys[lo] * (1.0 - t) + ys[hi] * t;
}

function stat(values) {
  const abs = values.map(Math.abs);
  return {
    max: Math.max(...abs),
    mean: abs.reduce((a, b) => a + b, 0.0) / abs.length,
    rms: Math.sqrt(values.reduce((a, b) => a + b * b, 0.0) / values.length),
  };
}

function csvEscape(value) {
  if (value === null || value === undefined) return "";
  const s = typeof value === "number" ? String(Number.isFinite(value) ? value : "") : String(value);
  return /[",\n]/.test(s) ? `"${s.replaceAll('"', '""')}"` : s;
}

function writeTable(fileName, rows) {
  const headers = Object.keys(rows[0]);
  const lines = [headers.join(",")];
  for (const row of rows) {
    lines.push(headers.map((h) => csvEscape(row[h])).join(","));
  }
  fs.writeFileSync(path.join(root, fileName), `${lines.join("\n")}\n`, "utf8");
}

function readTiming(algo) {
  try {
    const timing = JSON.parse(fs.readFileSync(path.join(root, "result", `result_${algo}`, "timing.json"), "utf8"));
    return {
      totalSeconds: timing.total_elapsed_s ?? timing.total_seconds ?? timing.totalSeconds ?? "",
      peakMemoryMb: timing.peak_memory_mb ?? (
        timing.peak_memory_bytes !== undefined ? timing.peak_memory_bytes / (1024 * 1024) : (
          timing.peakMemoryBytes !== undefined ? timing.peakMemoryBytes / (1024 * 1024) : ""
        )
      ),
    };
  } catch {
    return { totalSeconds: "", peakMemoryBytes: "" };
  }
}

const hfssRows = readCsv("wg_bp_filter_S_parameters.csv");
const hfssFreq = hfssRows.map((r) => Number(r["Freq [GHz]"]));
const hfssS11 = hfssRows.map((r) => Number(r["mag(S(1,1)) []"]));
const hfssS21 = hfssRows.map((r) => Number(r["mag(S(2,1)) []"]));
const hfssS11Db = hfssS11.map(magToDb);
const hfssS21Db = hfssS21.map(magToDb);

const runs = algorithms.map(([name, color]) => {
  const rows = readCsv(`result/result_${name}/s_parameters.csv`);
  const freq = rows.map((r) => Number(r.freq_Hz) / 1e9);
  const s11 = rows.map((r) => Math.hypot(Number(r.S11_real), Number(r.S11_imag)));
  const s21 = rows.map((r) => Math.hypot(Number(r.S21_real), Number(r.S21_imag)));
  const s11Db = s11.map(magToDb);
  const s21Db = s21.map(magToDb);
  const refS11 = freq.map((f) => interp(hfssFreq, hfssS11, f));
  const refS21 = freq.map((f) => interp(hfssFreq, hfssS21, f));
  const refS11Db = refS11.map(magToDb);
  const refS21Db = refS21.map(magToDb);
  const errS11Db = s11Db.map((v, i) => v - refS11Db[i]);
  const errS21Db = s21Db.map((v, i) => v - refS21Db[i]);
  const errS11Lin = s11.map((v, i) => v - refS11[i]);
  const errS21Lin = s21.map((v, i) => v - refS21[i]);
  return { name, color, freq, s11Db, s21Db, errS11Db, errS21Db, errS11Lin, errS21Lin, timing: readTiming(name) };
});

const summary = runs.map((run) => {
  const s11Db = stat(run.errS11Db);
  const s21Db = stat(run.errS21Db);
  const s11Lin = stat(run.errS11Lin);
  const s21Lin = stat(run.errS21Lin);
  return {
    algorithm: run.name,
    points: run.freq.length,
    max_abs_db_s11: s11Db.max,
    mean_abs_db_s11: s11Db.mean,
    rms_db_s11: s11Db.rms,
    max_abs_db_s21: s21Db.max,
    mean_abs_db_s21: s21Db.mean,
    rms_db_s21: s21Db.rms,
    max_abs_lin_s11: s11Lin.max,
    max_abs_lin_s21: s21Lin.max,
    total_seconds: run.timing.totalSeconds,
    peak_memory_mb: run.timing.peakMemoryMb,
  };
});
writeTable("hfss_fastsweep_error_summary.csv", summary);

const comparisonRows = [];
for (let i = 0; i < runs[0].freq.length; ++i) {
  const freq = runs[0].freq[i];
  const row = {
    freq_GHz: freq,
    hfss_s11_db: magToDb(interp(hfssFreq, hfssS11, freq)),
    hfss_s21_db: magToDb(interp(hfssFreq, hfssS21, freq)),
  };
  for (const run of runs) {
    row[`${run.name}_s11_db`] = run.s11Db[i];
    row[`${run.name}_s21_db`] = run.s21Db[i];
    row[`${run.name}_err_s11_db`] = run.errS11Db[i];
    row[`${run.name}_err_s21_db`] = run.errS21Db[i];
  }
  comparisonRows.push(row);
}
writeTable("hfss_fastsweep_comparison.csv", comparisonRows);

const width = 1320;
const height = 980;
const margin = { left: 76, right: 28, top: 64, bottom: 58 };
const gapX = 72;
const gapY = 88;
const panelW = (width - margin.left - margin.right - gapX) / 2;
const panelH = (height - margin.top - margin.bottom - gapY) / 2;
const fMin = 40.0;
const fMax = 43.0;

function range(values, padFrac) {
  let min = Math.min(...values);
  let max = Math.max(...values);
  if (Math.abs(max - min) < 1e-12) {
    min -= 1.0;
    max += 1.0;
  }
  const pad = (max - min) * padFrac;
  return [min - pad, max + pad];
}

function xScale(panel, freq) {
  return panel.x + ((freq - fMin) / (fMax - fMin)) * panel.w;
}

function yScale(panel, value, yrange) {
  return panel.y + panel.h - ((value - yrange[0]) / (yrange[1] - yrange[0])) * panel.h;
}

function pathLine(points, panel, yrange) {
  return points
    .map(([x, y], i) => `${i === 0 ? "M" : "L"}${xScale(panel, x).toFixed(2)},${yScale(panel, y, yrange).toFixed(2)}`)
    .join(" ");
}

function ticks(min, max, count) {
  return Array.from({ length: count }, (_, i) => min + ((max - min) * i) / (count - 1));
}

function esc(text) {
  return String(text).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;");
}

const hfssBandIndexes = hfssFreq.map((f, i) => [f, i]).filter(([f]) => f >= fMin && f <= fMax).map(([, i]) => i);
const yRanges = {
  s11: range([...hfssBandIndexes.map((i) => hfssS11Db[i]), ...runs.flatMap((r) => r.s11Db)], 0.04),
  s21: range([...hfssBandIndexes.map((i) => hfssS21Db[i]), ...runs.flatMap((r) => r.s21Db)], 0.04),
  errS11: range(runs.flatMap((r) => r.errS11Db), 0.12),
  errS21: range(runs.flatMap((r) => r.errS21Db), 0.12),
};

const panels = [
  { x: margin.left, y: margin.top, w: panelW, h: panelH, mode: "s11", title: "|S11| dB: HFSS vs fast sweep", ylabel: "dB" },
  { x: margin.left + panelW + gapX, y: margin.top, w: panelW, h: panelH, mode: "s21", title: "|S21| dB: HFSS vs fast sweep", ylabel: "dB" },
  { x: margin.left, y: margin.top + panelH + gapY, w: panelW, h: panelH, mode: "errS11", title: "S11 dB error vs HFSS", ylabel: "solver - HFSS [dB]" },
  { x: margin.left + panelW + gapX, y: margin.top + panelH + gapY, w: panelW, h: panelH, mode: "errS21", title: "S21 dB error vs HFSS", ylabel: "solver - HFSS [dB]" },
];

let svg = `<?xml version="1.0" encoding="UTF-8"?>\n`;
svg += `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">\n`;
svg += `<rect width="100%" height="100%" fill="#fff"/>\n`;
svg += `<text x="${width / 2}" y="30" text-anchor="middle" font-family="Segoe UI, Arial" font-size="22" font-weight="700">BP Filter 101-point Sweep: HFSS comparison</text>\n`;
svg += `<text x="${width / 2}" y="52" text-anchor="middle" font-family="Segoe UI, Arial" font-size="12" fill="#555">ALPS, AWE, GAWE, MGAWE, WCAWE; PARDISO Release; reference: wg_bp_filter_S_parameters.csv</text>\n`;

for (const panel of panels) {
  const yrange = yRanges[panel.mode];
  svg += `<g font-family="Segoe UI, Arial" font-size="11">\n`;
  svg += `<rect x="${panel.x}" y="${panel.y}" width="${panel.w}" height="${panel.h}" fill="#fff" stroke="#bbb"/>\n`;
  svg += `<text x="${panel.x + panel.w / 2}" y="${panel.y - 18}" text-anchor="middle" font-size="14" font-weight="600">${esc(panel.title)}</text>\n`;
  for (const tick of [40, 41, 42, 43]) {
    const x = xScale(panel, tick);
    svg += `<line x1="${x}" y1="${panel.y}" x2="${x}" y2="${panel.y + panel.h}" stroke="#eee"/>\n`;
    svg += `<text x="${x}" y="${panel.y + panel.h + 18}" text-anchor="middle">${tick.toFixed(0)}</text>\n`;
  }
  for (const tick of ticks(yrange[0], yrange[1], 6)) {
    const y = yScale(panel, tick, yrange);
    svg += `<line x1="${panel.x}" y1="${y}" x2="${panel.x + panel.w}" y2="${y}" stroke="#eee"/>\n`;
    svg += `<text x="${panel.x - 8}" y="${y + 4}" text-anchor="end">${tick.toFixed(panel.mode.startsWith("err") ? 2 : 1)}</text>\n`;
  }
  if (yrange[0] < 0.0 && yrange[1] > 0.0) {
    const y0 = yScale(panel, 0.0, yrange);
    svg += `<line x1="${panel.x}" y1="${y0}" x2="${panel.x + panel.w}" y2="${y0}" stroke="#777" stroke-dasharray="4 4"/>\n`;
  }
  svg += `<text x="${panel.x + panel.w / 2}" y="${panel.y + panel.h + 40}" text-anchor="middle">Frequency [GHz]</text>\n`;
  svg += `<text transform="translate(${panel.x - 54},${panel.y + panel.h / 2}) rotate(-90)" text-anchor="middle">${esc(panel.ylabel)}</text>\n`;
  if (panel.mode === "s11" || panel.mode === "s21") {
    const hfssPoints = hfssBandIndexes.map((i) => [hfssFreq[i], panel.mode === "s11" ? hfssS11Db[i] : hfssS21Db[i]]);
    svg += `<path d="${pathLine(hfssPoints, panel, yrange)}" fill="none" stroke="#111" stroke-width="2.4"/>\n`;
    for (const run of runs) {
      const points = run.freq.map((f, i) => [f, panel.mode === "s11" ? run.s11Db[i] : run.s21Db[i]]);
      svg += `<path d="${pathLine(points, panel, yrange)}" fill="none" stroke="${run.color}" stroke-width="1.35" opacity="0.95"/>\n`;
    }
  } else {
    for (const run of runs) {
      const points = run.freq.map((f, i) => [f, panel.mode === "errS11" ? run.errS11Db[i] : run.errS21Db[i]]);
      svg += `<path d="${pathLine(points, panel, yrange)}" fill="none" stroke="${run.color}" stroke-width="1.5" opacity="0.95"/>\n`;
    }
  }
  svg += `</g>\n`;
}

let legendX = margin.left;
const legendY = height - 22;
svg += `<g font-family="Segoe UI, Arial" font-size="12">\n`;
svg += `<line x1="${legendX}" y1="${legendY}" x2="${legendX + 28}" y2="${legendY}" stroke="#111" stroke-width="2.4"/><text x="${legendX + 36}" y="${legendY + 4}">HFSS</text>\n`;
legendX += 92;
for (const run of runs) {
  svg += `<line x1="${legendX}" y1="${legendY}" x2="${legendX + 28}" y2="${legendY}" stroke="${run.color}" stroke-width="2"/><text x="${legendX + 36}" y="${legendY + 4}">${run.name}</text>\n`;
  legendX += 108;
}
svg += `</g>\n</svg>\n`;

fs.writeFileSync(path.join(root, "hfss_fastsweep_comparison.svg"), svg, "utf8");
console.log("Wrote hfss_fastsweep_comparison.svg");
console.log("Wrote hfss_fastsweep_comparison.csv");
console.log("Wrote hfss_fastsweep_error_summary.csv");
