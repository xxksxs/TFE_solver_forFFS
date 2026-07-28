import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const outputDir = path.join(root, "result", "result_COMPARE_IOStructure");
const hfssPath = path.join(root, "IOStructure_S_parameters.csv");
const includeAlps = process.argv.includes("--include-alps");
const runs = [
  ...(includeAlps ? [{
    name: "ALPS",
    color: "#009E73",
    directory: path.join(root, "result", "result_ALPS"),
    configuration: "single 95 GHz expansion, two-sided Lanczos-Pade, order 12",
    candidateCount: 24,
  }] : []),
  {
    name: "AWE",
    color: "#CC79A7",
    directory: path.join(root, "result", "result_AWE"),
    configuration: "single expansion at 95 GHz, Pade [11/12], 24 moments",
    candidateCount: 24,
  },
  {
    name: "GAWE",
    color: "#E69F00",
    directory: path.join(root, "result", "result_GAWE"),
    configuration: "single expansion at 95 GHz, Galerkin order 12",
    candidateCount: 12,
  },
  {
    name: "WCAWE",
    color: "#D55E00",
    directory: path.join(root, "result", "result_WCAWE"),
    configuration: "single expansion at 95 GHz, order 12",
    candidateCount: 12,
  },
  {
    name: "MGAWE",
    color: "#0072B2",
    directory: path.join(root, "result", "result_MGAWE"),
    configuration: "3 expansions at 90/95/100 GHz, local order 4",
    candidateCount: 12,
  },
];
const outputStem = includeAlps
  ? "ios_alps_awe_gawe_wcawe_mgawe"
  : "ios_awe_gawe_wcawe_mgawe";

// 解析含引号字段的 CSV，避免依赖字段中不出现逗号这一偶然条件。
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
  if (field.length > 0 || row.length > 0) {
    row.push(field);
    rows.push(row);
  }
  const header = rows.shift().map((value) => value.trim().replace(/^"|"$/g, ""));
  return rows
    .filter((values) => values.some((value) => value.trim() !== ""))
    .map((values) => Object.fromEntries(header.map((key, index) => [key, (values[index] ?? "").trim()])));
}

// 读取 CSV 并保留列名，后续会对关键字段做显式数值校验。
function readCsv(filePath) {
  return parseCsv(fs.readFileSync(filePath, "utf8"));
}

// 把复数实部和虚部转换为线性幅值。
function complexMagnitude(real, imag) {
  return Math.hypot(Number(real), Number(imag));
}

// 将线性幅值转换为 dB，并限制对数下界以避免负无穷。
function magnitudeToDb(value) {
  return 20.0 * Math.log10(Math.max(Number(value), 1e-300));
}

// 计算带符号误差对应的最大绝对值、平均绝对值和均方根。
function statistics(values) {
  const absolute = values.map(Math.abs);
  return {
    maxAbs: Math.max(...absolute),
    meanAbs: absolute.reduce((sum, value) => sum + value, 0.0) / absolute.length,
    rms: Math.sqrt(values.reduce((sum, value) => sum + value * value, 0.0) / values.length),
  };
}

// 计算逐点有符号相对误差；分母下限只用于防止参考值为数值零。
function relativeErrors(candidate, reference, denominatorFloor = 1e-12) {
  return candidate.map(
    (value, index) => (value - reference[index]) / Math.max(Math.abs(reference[index]), denominatorFloor),
  );
}

// 全局相对 L2 误差比逐点最大值更不易被单个深陷点支配。
function relativeL2Error(candidate, reference) {
  const errorSquared = candidate.reduce(
    (sum, value, index) => sum + (value - reference[index]) ** 2,
    0.0,
  );
  const referenceSquared = reference.reduce((sum, value) => sum + value ** 2, 0.0);
  return Math.sqrt(errorSquared / Math.max(referenceSquared, 1e-24));
}

// 对 CSV 单元格做必要转义。
function csvEscape(value) {
  const text = typeof value === "number"
    ? (Number.isFinite(value) ? String(value) : "")
    : String(value ?? "");
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

// 以对象键作为列名写出统一格式的 CSV。
function writeCsv(filePath, rows) {
  if (rows.length === 0) {
    throw new Error(`Cannot write empty CSV: ${filePath}`);
  }
  const headers = Object.keys(rows[0]);
  const lines = [
    headers.join(","),
    ...rows.map((row) => headers.map((header) => csvEscape(row[header])).join(",")),
  ];
  fs.writeFileSync(filePath, `${lines.join("\n")}\n`, "utf8");
}

// 从 timing.json 中提取总时间、求解阶段时间和峰值内存。
function readTiming(directory) {
  const timing = JSON.parse(fs.readFileSync(path.join(directory, "timing.json"), "utf8"));
  const solvePhase = timing.phases?.find((phase) => phase.name === "Frequency-domain solve");
  return {
    totalSeconds: Number(timing.total_elapsed_s),
    solveSeconds: Number(solvePhase?.elapsed_s ?? Number.NaN),
    peakMemoryMb: Number(timing.peak_memory_mb),
    offlineBuildSeconds: Number(timing.offline_build_s ?? Number.NaN),
    orthogonalizationSeconds: Number(timing.orthogonalization_s ?? Number.NaN),
    romProjectionSeconds: Number(timing.rom_projection_s ?? Number.NaN),
    onlineSweepSeconds: Number(timing.online_sweep_s ?? Number.NaN),
    symbolicAnalysisCount: Number(timing.symbolic_analysis_count ?? 0),
    symbolicAnalysisSeconds: Number(timing.symbolic_analysis_s ?? 0),
    numericFactorizationCount: Number(timing.numeric_factorization_count ?? 0),
    numericFactorizationSeconds: Number(timing.numeric_factorization_s ?? 0),
    factorizedRhsSolveCount: Number(timing.factorized_rhs_solve_count ?? 0),
    factorizedRhsSolveSeconds: Number(timing.factorized_rhs_solve_s ?? 0),
  };
}

// 读取算法诊断信息，确保比较中使用实际 ROM 维数和正交性。
function readDiagnostics(directory) {
  return JSON.parse(fs.readFileSync(path.join(directory, "diagnostics.json"), "utf8"));
}

// 读取 WCAWE 条件数曲线最后一行，衡量传统矩基与正交基的稳定性差异。
function readFinalBasisCondition(directory) {
  const rows = readCsv(path.join(directory, "basis_condition.csv"));
  const row = rows.at(-1);
  return {
    order: Number(row.order),
    aweConditionProxy: Number(row.awe_condition_proxy),
    wcaweConditionProxy: Number(row.wcawe_condition_proxy),
  };
}

// 验证 HFSS 与两个 ROM 使用完全一致的 101 个频点。
function assertAligned(reference, candidate, name) {
  if (reference.length !== candidate.length) {
    throw new Error(`${name}: point count ${candidate.length} differs from HFSS ${reference.length}`);
  }
  for (let i = 0; i < reference.length; ++i) {
    if (Math.abs(reference[i].freqGHz - candidate[i].freqGHz) > 1e-9) {
      throw new Error(`${name}: frequency mismatch at row ${i + 1}`);
    }
  }
}

const hfss = readCsv(hfssPath).map((row) => ({
  freqGHz: Number(row["Freq [GHz]"]),
  s11: Number(row["mag(S(1,1)) []"]),
  s21: Number(row["mag(S(2,1)) []"]),
}));

if (hfss.length !== 101 || hfss.some((row) => !Number.isFinite(row.freqGHz + row.s11 + row.s21))) {
  throw new Error("HFSS reference must contain 101 finite S11/S21 rows");
}

for (const run of runs) {
  run.rows = readCsv(path.join(run.directory, "s_parameters.csv")).map((row) => ({
    freqGHz: Number(row.freq_Hz) / 1e9,
    s11: complexMagnitude(row.S11_real, row.S11_imag),
    s21: complexMagnitude(row.S21_real, row.S21_imag),
  }));
  assertAligned(hfss, run.rows, run.name);
  run.timing = readTiming(run.directory);
  run.diagnostics = readDiagnostics(run.directory);
}

const runByName = Object.fromEntries(runs.map((run) => [run.name, run]));
runByName.WCAWE.basisCondition = readFinalBasisCondition(runByName.WCAWE.directory);
fs.mkdirSync(outputDir, { recursive: true });
for (const staleStem of [
  "ios_wcawe_mgawe",
  "ios_alps_wcawe_mgawe",
  "ios_awe_gawe_wcawe_mgawe",
  "ios_alps_awe_gawe_wcawe_mgawe",
]) {
  for (const suffix of ["comparison.csv", "summary.csv", "report.json", "comparison.svg"]) {
    fs.rmSync(path.join(outputDir, staleStem + "_" + suffix), { force: true });
  }
}

const comparisonRows = hfss.map((reference, index) => {
  const row = {
    freq_GHz: reference.freqGHz,
    HFSS_S11_mag: reference.s11,
    HFSS_S21_mag: reference.s21,
    HFSS_S11_dB: magnitudeToDb(reference.s11),
    HFSS_S21_dB: magnitudeToDb(reference.s21),
  };
  for (const run of runs) {
    const current = run.rows[index];
    row[`${run.name}_S11_mag`] = current.s11;
    row[`${run.name}_S21_mag`] = current.s21;
    row[`${run.name}_S11_dB`] = magnitudeToDb(current.s11);
    row[`${run.name}_S21_dB`] = magnitudeToDb(current.s21);
    row[`${run.name}_S11_error_mag`] = current.s11 - reference.s11;
    row[`${run.name}_S21_error_mag`] = current.s21 - reference.s21;
    row[`${run.name}_S11_relative_error_percent`] =
      100.0 * (current.s11 - reference.s11) / Math.max(Math.abs(reference.s11), 1e-12);
    row[`${run.name}_S21_relative_error_percent`] =
      100.0 * (current.s21 - reference.s21) / Math.max(Math.abs(reference.s21), 1e-12);
    row[`${run.name}_S11_error_dB`] = magnitudeToDb(current.s11) - magnitudeToDb(reference.s11);
    row[`${run.name}_S21_error_dB`] = magnitudeToDb(current.s21) - magnitudeToDb(reference.s21);
  }
  return row;
});

const summaryRows = runs.map((run) => {
  const s11Linear = statistics(run.rows.map((row, index) => row.s11 - hfss[index].s11));
  const s21Linear = statistics(run.rows.map((row, index) => row.s21 - hfss[index].s21));
  const s11Db = statistics(run.rows.map(
    (row, index) => magnitudeToDb(row.s11) - magnitudeToDb(hfss[index].s11),
  ));
  const s21Db = statistics(run.rows.map(
    (row, index) => magnitudeToDb(row.s21) - magnitudeToDb(hfss[index].s21),
  ));
  const s11Relative = statistics(relativeErrors(
    run.rows.map((row) => row.s11),
    hfss.map((row) => row.s11),
  ));
  const s21Relative = statistics(relativeErrors(
    run.rows.map((row) => row.s21),
    hfss.map((row) => row.s21),
  ));
  const s11RelativeL2 = relativeL2Error(
    run.rows.map((row) => row.s11),
    hfss.map((row) => row.s11),
  );
  const s21RelativeL2 = relativeL2Error(
    run.rows.map((row) => row.s21),
    hfss.map((row) => row.s21),
  );
  return {
    algorithm: run.name,
    configuration: run.configuration,
    points: run.rows.length,
    candidate_vectors_or_moments: run.candidateCount,
    rom_dimension: run.diagnostics.rom_dimension,
    retained_columns: run.diagnostics.retained_columns,
    deflated_columns: run.diagnostics.deflated_columns,
    basis_orthogonality_error: run.diagnostics.basis_orthogonality_error,
    lanczos_biorthogonality_error: run.diagnostics.lanczos_biorthogonality_error ?? 0,
    lanczos_tridiagonal_leakage: run.diagnostics.lanczos_tridiagonal_leakage ?? 0,
    lanczos_breakdown_detected: run.diagnostics.lanczos_breakdown_detected ?? false,
    pade_input_pivot_ratio: run.diagnostics.pade_input_pivot_ratio,
    pade_output_pivot_ratio: run.diagnostics.pade_output_pivot_ratio,
    total_seconds: run.timing.totalSeconds,
    solve_seconds: run.timing.solveSeconds,
    offline_build_seconds: run.timing.offlineBuildSeconds,
    orthogonalization_seconds: run.timing.orthogonalizationSeconds,
    rom_projection_seconds: run.timing.romProjectionSeconds,
    online_sweep_seconds: run.timing.onlineSweepSeconds,
    symbolic_analysis_count: run.timing.symbolicAnalysisCount,
    symbolic_analysis_seconds: run.timing.symbolicAnalysisSeconds,
    numeric_factorization_count: run.timing.numericFactorizationCount,
    numeric_factorization_seconds: run.timing.numericFactorizationSeconds,
    factorized_rhs_solve_count: run.timing.factorizedRhsSolveCount,
    factorized_rhs_solve_seconds: run.timing.factorizedRhsSolveSeconds,
    peak_memory_mb: run.timing.peakMemoryMb,
    HFSS_S11_max_abs_mag_error: s11Linear.maxAbs,
    HFSS_S11_mean_abs_mag_error: s11Linear.meanAbs,
    HFSS_S11_rms_mag_error: s11Linear.rms,
    HFSS_S21_max_abs_mag_error: s21Linear.maxAbs,
    HFSS_S21_mean_abs_mag_error: s21Linear.meanAbs,
    HFSS_S21_rms_mag_error: s21Linear.rms,
    HFSS_S11_max_abs_relative_error_percent: 100.0 * s11Relative.maxAbs,
    HFSS_S11_mean_abs_relative_error_percent: 100.0 * s11Relative.meanAbs,
    HFSS_S11_rms_relative_error_percent: 100.0 * s11Relative.rms,
    HFSS_S11_relative_L2_error_percent: 100.0 * s11RelativeL2,
    HFSS_S21_max_abs_relative_error_percent: 100.0 * s21Relative.maxAbs,
    HFSS_S21_mean_abs_relative_error_percent: 100.0 * s21Relative.meanAbs,
    HFSS_S21_rms_relative_error_percent: 100.0 * s21Relative.rms,
    HFSS_S21_relative_L2_error_percent: 100.0 * s21RelativeL2,
    HFSS_S11_max_abs_dB_error: s11Db.maxAbs,
    HFSS_S21_max_abs_dB_error: s21Db.maxAbs,
  };
});

// 计算任意两种 ROM 的线性幅值差，区分模型偏差和算法之间的差异。
function pairwiseDifference(leftName, rightName) {
  const left = runByName[leftName];
  const right = runByName[rightName];
  return {
    S11: statistics(left.rows.map((row, index) => row.s11 - right.rows[index].s11)),
    S21: statistics(left.rows.map((row, index) => row.s21 - right.rows[index].s21)),
  };
}

const pairwiseMagnitudeDifferences = {};
for (let left = 0; left < runs.length; ++left) {
  for (let right = left + 1; right < runs.length; ++right) {
    const leftName = runs[left].name;
    const rightName = runs[right].name;
    pairwiseMagnitudeDifferences[`${leftName}_vs_${rightName}`] =
      pairwiseDifference(leftName, rightName);
  }
}

const report = {
  case: "IOStructure",
  frequencyRangeGHz: [hfss[0].freqGHz, hfss.at(-1).freqGHz],
  points: hfss.length,
  reference: "IOStructure_S_parameters.csv",
  algorithms: summaryRows,
  pairwiseMagnitudeDifferences,
  wcaweConditioning: {
    order: runByName.WCAWE.basisCondition.order,
    aweConditionProxy: runByName.WCAWE.basisCondition.aweConditionProxy,
    wcaweConditionProxy: runByName.WCAWE.basisCondition.wcaweConditionProxy,
    momentReconstructionError: runByName.WCAWE.diagnostics.wcawe_moment_reconstruction_error,
  },
  notes: [
    "Primary accuracy metrics use pointwise relative magnitude error: (ROM - HFSS) / max(abs(HFSS), 1e-12).",
    "Global relative L2 error is also reported because pointwise relative error can be amplified near deep S-parameter nulls.",
    "Absolute magnitude and dB errors are retained as secondary diagnostics.",
    "AWE is a scalar Pade approximation rather than a Galerkin ROM, so rom_dimension is zero by design.",
    "AWE order 12 generates 24 full-system moments; GAWE, WCAWE and MGAWE use 12 candidate vectors and 12 factorized RHS solves.",
    ...(includeAlps ? ["ALPS builds independent S11/S21 two-sided Lanczos-Pade models at scalar order q=12; each reduced model has dimension 12 and the two sequences retain 24 columns in total."] : []),
    "The four-port model requires all outgoing ports for a full passivity balance; S11 and S21 alone are not used as a passivity metric.",
  ],
};

writeCsv(path.join(outputDir, outputStem + "_comparison.csv"), comparisonRows);
writeCsv(path.join(outputDir, outputStem + "_summary.csv"), summaryRows);
fs.writeFileSync(
  path.join(outputDir, outputStem + "_report.json"),
  `${JSON.stringify(report, null, 2)}\n`,
  "utf8",
);

// 将两种算法、HFSS 曲线和相对幅值误差绘制在一张四面板 SVG 中。
const width = 1360;
const height = 940;
const margin = { left: 80, right: 32, top: 78, bottom: 60 };
const gapX = 78;
const gapY = 92;
const panelWidth = (width - margin.left - margin.right - gapX) / 2;
const panelHeight = (height - margin.top - margin.bottom - gapY) / 2;
const fMin = hfss[0].freqGHz;
const fMax = hfss.at(-1).freqGHz;

function paddedRange(values, fraction = 0.08, includeZero = false) {
  let min = Math.min(...values);
  let max = Math.max(...values);
  if (includeZero) {
    min = Math.min(min, 0);
    max = Math.max(max, 0);
  }
  if (Math.abs(max - min) < 1e-15) {
    min -= 1;
    max += 1;
  }
  const padding = (max - min) * fraction;
  return [min - padding, max + padding];
}

function xScale(panel, value) {
  return panel.x + ((value - fMin) / (fMax - fMin)) * panel.width;
}

function yScale(panel, value) {
  return panel.y + panel.height
    - ((value - panel.range[0]) / (panel.range[1] - panel.range[0])) * panel.height;
}

function linePath(panel, values) {
  return values.map((value, index) => {
    const prefix = index === 0 ? "M" : "L";
    return `${prefix}${xScale(panel, hfss[index].freqGHz).toFixed(2)},${yScale(panel, value).toFixed(2)}`;
  }).join(" ");
}

function ticks(min, max, count) {
  return Array.from({ length: count }, (_, index) => min + ((max - min) * index) / (count - 1));
}

const panels = [
  {
    x: margin.left,
    y: margin.top,
    width: panelWidth,
    height: panelHeight,
    title: "|S11| linear magnitude",
    range: paddedRange([
      ...hfss.map((row) => row.s11),
      ...runs.flatMap((run) => run.rows.map((row) => row.s11)),
    ]),
    reference: hfss.map((row) => row.s11),
    values: runs.map((run) => run.rows.map((row) => row.s11)),
    error: false,
  },
  {
    x: margin.left + panelWidth + gapX,
    y: margin.top,
    width: panelWidth,
    height: panelHeight,
    title: "|S21| linear magnitude",
    range: paddedRange([
      ...hfss.map((row) => row.s21),
      ...runs.flatMap((run) => run.rows.map((row) => row.s21)),
    ]),
    reference: hfss.map((row) => row.s21),
    values: runs.map((run) => run.rows.map((row) => row.s21)),
    error: false,
  },
  {
    x: margin.left,
    y: margin.top + panelHeight + gapY,
    width: panelWidth,
    height: panelHeight,
    title: "S11 relative magnitude error versus HFSS [%]",
    range: paddedRange(runs.flatMap(
      (run) => relativeErrors(
        run.rows.map((row) => row.s11),
        hfss.map((row) => row.s11),
      ).map((value) => 100.0 * value),
    ), 0.1, true),
    values: runs.map((run) => relativeErrors(
      run.rows.map((row) => row.s11),
      hfss.map((row) => row.s11),
    ).map((value) => 100.0 * value)),
    error: true,
  },
  {
    x: margin.left + panelWidth + gapX,
    y: margin.top + panelHeight + gapY,
    width: panelWidth,
    height: panelHeight,
    title: "S21 relative magnitude error versus HFSS [%]",
    range: paddedRange(runs.flatMap(
      (run) => relativeErrors(
        run.rows.map((row) => row.s21),
        hfss.map((row) => row.s21),
      ).map((value) => 100.0 * value),
    ), 0.1, true),
    values: runs.map((run) => relativeErrors(
      run.rows.map((row) => row.s21),
      hfss.map((row) => row.s21),
    ).map((value) => 100.0 * value)),
    error: true,
  },
];

const titleAlgorithms = runs.map((run) => run.name).join(", ");
let svg = `<?xml version="1.0" encoding="UTF-8"?>\n`;
svg += `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">\n`;
svg += `<rect width="100%" height="100%" fill="#ffffff"/>\n`;
svg += `<text x="${width / 2}" y="31" text-anchor="middle" font-family="Segoe UI, Arial" font-size="23" font-weight="700">IOStructure: ${titleAlgorithms} versus HFSS</text>\n`;
svg += `<text x="${width / 2}" y="55" text-anchor="middle" font-family="Segoe UI, Arial" font-size="12" fill="#555555">90-100 GHz, 101 points, PARDISO Release, production accuracy configurations</text>\n`;

for (const panel of panels) {
  svg += `<g font-family="Segoe UI, Arial" font-size="11">\n`;
  svg += `<rect x="${panel.x}" y="${panel.y}" width="${panel.width}" height="${panel.height}" fill="#ffffff" stroke="#b8b8b8"/>\n`;
  svg += `<text x="${panel.x + panel.width / 2}" y="${panel.y - 17}" text-anchor="middle" font-size="14" font-weight="600">${panel.title}</text>\n`;
  for (const tick of ticks(fMin, fMax, 6)) {
    const x = xScale(panel, tick);
    svg += `<line x1="${x}" y1="${panel.y}" x2="${x}" y2="${panel.y + panel.height}" stroke="#eeeeee"/>\n`;
    svg += `<text x="${x}" y="${panel.y + panel.height + 19}" text-anchor="middle">${tick.toFixed(0)}</text>\n`;
  }
  for (const tick of ticks(panel.range[0], panel.range[1], 6)) {
    const y = yScale(panel, tick);
    svg += `<line x1="${panel.x}" y1="${y}" x2="${panel.x + panel.width}" y2="${y}" stroke="#eeeeee"/>\n`;
    svg += `<text x="${panel.x - 9}" y="${y + 4}" text-anchor="end">${tick.toExponential(2)}</text>\n`;
  }
  if (panel.error && panel.range[0] <= 0 && panel.range[1] >= 0) {
    const y = yScale(panel, 0);
    svg += `<line x1="${panel.x}" y1="${y}" x2="${panel.x + panel.width}" y2="${y}" stroke="#777777" stroke-dasharray="5 4"/>\n`;
  }
  if (!panel.error) {
    svg += `<path d="${linePath(panel, panel.reference)}" fill="none" stroke="#222222" stroke-width="2.2"/>\n`;
  }
  panel.values.forEach((values, index) => {
    svg += `<path d="${linePath(panel, values)}" fill="none" stroke="${runs[index].color}" stroke-width="1.8"/>\n`;
  });
  svg += `<text x="${panel.x + panel.width / 2}" y="${panel.y + panel.height + 39}" text-anchor="middle">Frequency [GHz]</text>\n`;
  svg += `</g>\n`;
}

const legendY = height - 24;
svg += `<g font-family="Segoe UI, Arial" font-size="12">\n`;
svg += `<line x1="430" y1="${legendY}" x2="466" y2="${legendY}" stroke="#222222" stroke-width="2.2"/><text x="474" y="${legendY + 4}">HFSS</text>\n`;
runs.forEach((run, index) => {
  const x = 545 + index * 145;
  svg += `<line x1="${x}" y1="${legendY}" x2="${x + 36}" y2="${legendY}" stroke="${run.color}" stroke-width="2.2"/><text x="${x + 44}" y="${legendY + 4}">${run.name}</text>\n`;
});
svg += `</g>\n</svg>\n`;

fs.writeFileSync(path.join(outputDir, outputStem + "_comparison.svg"), svg, "utf8");

console.log(JSON.stringify(report, null, 2));
