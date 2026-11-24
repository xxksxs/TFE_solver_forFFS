import fs from "node:fs";
import path from "node:path";
import { randomUUID } from "node:crypto";

const channels = ["S11", "S21"];
const styles = {
  WCAWE: { color: "#b5471c", dash: "", marker: "circle" },
  MGAWE: { color: "#0066a1", dash: ' stroke-dasharray="7 4"', marker: "square" },
};

export function readRunReport(directory, name, allowLegacy) {
  if (path.basename(directory).startsWith(".bpfem-pending-")
      || fs.existsSync(path.join(directory, "run.json.pending"))) {
    throw new Error(`${directory}: solver run has not been published`);
  }
  const report = JSON.parse(fs.readFileSync(path.join(directory, "run.json"), "utf8"));
  if (report.status !== "ok" || report.sweep_strategy !== name.toLowerCase()) {
    throw new Error(`${directory}: expected a successful ${name} run`);
  }
  const validHash = (value) => typeof value === "string" && /^[a-f0-9]{64}$/i.test(value);
  if (report.schema >= 2 || report.inputs !== undefined) {
    if (!validHash(report.inputs?.aedt_sha256) || !validHash(report.inputs?.mesh_sha256)) {
      throw new Error(`${directory}: missing or invalid recorded input SHA-256`);
    }
  } else if (!allowLegacy) {
    throw new Error(`${directory}: legacy report has no input hashes; use --allow-legacy to compare recorded paths`);
  }
  return report;
}

export function validateReportedGrid(report, frequencies, directory) {
  if (report.schema < 2 && !Array.isArray(report.sweep)) return;
  if (!Array.isArray(report.sweep) || report.sweep.length !== frequencies.length
      || frequencies.some((f, i) => !Number.isFinite(f) || report.sweep[i].freq_hz !== f)) {
    throw new Error(`${directory}: frequency samples differ between CSV and run.json`);
  }
}

function readRun(directory, name, readCsv, allowLegacy) {
  const report = readRunReport(directory, name, allowLegacy);
  const file = path.join(directory, "s_parameters.csv");
  const rows = readCsv(file);
  if (!rows.length) throw new Error(`${file}: no frequency samples`);
  const readNumber = (row, column, index) => {
    const value = row[column];
    if (value === undefined || value === "" || !Number.isFinite(Number(value))) {
      throw new Error(`${file}: invalid ${column} at data row ${index + 1}`);
    }
    return Number(value);
  };
  const samples = rows.map((row, index) => ({
    frequencyHz: readNumber(row, "freq_Hz", index),
    ...Object.fromEntries(channels.map((channel) => [channel, [
      readNumber(row, `${channel}_real`, index),
      readNumber(row, `${channel}_imag`, index),
    ]])),
  }));
  if (samples.some(({ frequencyHz }, index) =>
    frequencyHz <= 0 || (index > 0 && frequencyHz <= samples[index - 1].frequencyHz))) {
    throw new Error(`${file}: frequencies must be positive and strictly increasing`);
  }
  const diagnostics = name === "DIRECT" ? null
    : JSON.parse(fs.readFileSync(path.join(directory, "diagnostics.json"), "utf8"));
  if (diagnostics?.reduced_solve_succeeded === false) {
    throw new Error(`${directory}: a reduced solve failed`);
  }
  validateReportedGrid(report, samples.map(({ frequencyHz }) => frequencyHz), directory);
  return { name, directory, report, diagnostics, samples };
}

export function modelSettings({ report }) {
  const { aedt, mesh, basisOrder, portMethod, tfeModesPerPort } = report.options ?? {};
  if (![aedt, mesh, report.cwd].every((value) => typeof value === "string" && value.length > 0)
      || basisOrder == null || basisOrder === "" || !Number.isInteger(Number(basisOrder))
      || !["numerical", "analytic", "tfe"].includes(portMethod)
      || (portMethod === "tfe" && !(Number(tfeModesPerPort) > 0))) {
    throw new Error("Run report is missing the input paths or physical settings");
  }
  const paths = /^[A-Za-z]:[\\/]/.test(report.cwd) ? path.win32 : path;
  return {
    aedt: paths.resolve(report.cwd, aedt),
    mesh: paths.resolve(report.cwd, mesh),
    basisOrder: String(basisOrder),
    portMethod,
    ...(portMethod === "tfe" ? { tfeModesPerPort: String(tfeModesPerPort) } : {}),
  };
}

export function compareModelSettings(reference, model) {
  const settings = modelSettings(reference);
  const candidate = modelSettings(model);
  for (const key of ["basisOrder", "portMethod", "tfeModesPerPort"]) {
    if (settings[key] !== candidate[key]) {
      throw new Error(`${model.name}: physical setting ${key} differs from ${reference.name}`);
    }
  }
  if (reference.report.inputs && model.report.inputs) {
    for (const key of ["aedt_sha256", "mesh_sha256"]) {
      if (reference.report.inputs[key].toLowerCase() !== model.report.inputs[key].toLowerCase()) {
        throw new Error(`${model.name}: ${key} differs from ${reference.name}`);
      }
    }
    return "sha256";
  }
  for (const key of ["aedt", "mesh"]) {
    if (settings[key] !== candidate[key]) {
      throw new Error(`${model.name}: legacy ${key} path differs from ${reference.name}`);
    }
  }
  return "legacy-path-only";
}

function assertComparable(direct, model) {
  const inputVerification = compareModelSettings(direct, model);
  if (model.samples.length !== direct.samples.length) {
    throw new Error(`${model.name}: frequency count differs from Direct`);
  }
  model.samples.forEach(({ frequencyHz }, index) => {
    if (frequencyHz !== direct.samples[index].frequencyHz) {
      throw new Error(`${model.name}: frequency differs from Direct at data row ${index + 1}`);
    }
  });
  return inputVerification;
}

export function publishComparison(root, name, writeOutputs) {
  fs.mkdirSync(root, { recursive: true });
  const staging = fs.mkdtempSync(path.join(root, ".bpfem-compare-pending-"));
  writeOutputs(staging);
  const preferred = path.join(root, name);
  const alternate = `${preferred}-${randomUUID()}`;
  // Published folders are non-empty and a rename cannot replace them.
  if (!fs.existsSync(preferred)) {
    try {
      fs.renameSync(staging, preferred);
      return preferred;
    } catch (error) {
      if (!fs.existsSync(preferred)) throw error;
    }
  }
  fs.renameSync(staging, alternate);
  return alternate;
}

function complexDifference([leftReal, leftImag], [rightReal, rightImag]) {
  return Math.hypot(leftReal - rightReal, leftImag - rightImag);
}

function summarizeErrors(absolute, relative, referenceMagnitudes, frequencies) {
  const valid = relative.filter((value) => value !== null);
  const referenceNorm = Math.hypot(...referenceMagnitudes);
  const maxAbsoluteError = Math.max(...absolute);
  const maxRelativeError = valid.length ? Math.max(...valid) : null;
  return {
    maxAbsoluteError,
    maxAbsoluteErrorFrequencyHz: frequencies[absolute.indexOf(maxAbsoluteError)],
    maxRelativeError,
    maxRelativeErrorFrequencyHz: maxRelativeError === null
      ? null : frequencies[relative.indexOf(maxRelativeError)],
    relativeL2Error: referenceNorm > 0 ? Math.hypot(...absolute) / referenceNorm : null,
    zeroReferenceCount: relative.length - valid.length,
  };
}

function drawComparison(rows, names) {
  const width = 1160;
  const height = 535;
  const top = 111;
  const panelWidth = 474;
  const panelHeight = 302;
  const frequencies = rows.map((row) => row.freq_Hz / 1e9);
  const interval = frequencies.at(-1) - frequencies[0];
  const first = interval > 0 ? frequencies[0] : frequencies[0] - 0.5;
  const last = interval > 0 ? frequencies.at(-1) : frequencies[0] + 0.5;
  const allErrors = channels.flatMap((channel) => names.flatMap((name) =>
    rows.map((row) => row[`${name}_${channel}_relative_error`])));
  const positive = allErrors.filter((value) => value > 0);
  const hasZeros = allErrors.includes(0);
  const minPower = positive.length
    ? Math.floor(Math.log10(Math.min(...positive))) - Number(hasZeros) : -16;
  const maxPower = positive.length
    ? Math.max(minPower + 1, Math.ceil(Math.log10(Math.max(...positive)))) : -15;
  const tickStep = Math.max(1, Math.ceil((maxPower - minPower) / 8));
  const y = (power) => top + (maxPower - power) * panelHeight / (maxPower - minPower);
  const svg = [
    `<?xml version="1.0" encoding="UTF-8"?>`,
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-label="Complex scattering parameter relative errors against Direct FEM">`,
    `<rect width="100%" height="100%" fill="#fff"/>`,
    `<g font-family="Segoe UI, Arial, sans-serif" fill="#16242e">`,
    `<text x="580" y="39" text-anchor="middle" font-size="25" font-weight="600">MGAWE and WCAWE against Direct FEM</text>`,
    `<text x="580" y="67" text-anchor="middle" font-size="14" fill="#465969">Complex relative error · |S_ROM − S_Direct| / |S_Direct|</text>`,
  ];
  channels.forEach((channel, channelIndex) => {
    const left = 90 + channelIndex * 548;
    const x = (frequency) => left + (frequency - first) * panelWidth / (last - first);
    svg.push(`<text x="${left + panelWidth / 2}" y="101" text-anchor="middle" font-size="17" font-weight="600">${channel}</text>`);
    svg.push(`<rect x="${left}" y="${top}" width="${panelWidth}" height="${panelHeight}" fill="#fbfcfd" stroke="#aebbc4"/>`);
    for (let power = minPower; power <= maxPower; power += tickStep) {
      const py = y(power);
      const tick = hasZeros && power === minPower ? "0 (floor)" : `1e${power}`;
      svg.push(`<line x1="${left}" y1="${py}" x2="${left + panelWidth}" y2="${py}" stroke="#dfe6ec"/>`);
      svg.push(`<text x="${left - 9}" y="${py + 4}" text-anchor="end" font-size="12">${tick}</text>`);
    }
    for (let index = 0; index <= 5; ++index) {
      const frequency = first + (last - first) * index / 5;
      const px = x(frequency);
      svg.push(`<line x1="${px}" y1="${top}" x2="${px}" y2="${top + panelHeight}" stroke="#dfe6ec"/>`);
      svg.push(`<text x="${px}" y="${top + panelHeight + 20}" text-anchor="middle" font-size="12">${Number(frequency.toPrecision(7))}</text>`);
    }
    names.forEach((name) => {
      const series = rows.map((row) => row[`${name}_${channel}_relative_error`]);
      const { color, dash, marker } = styles[name];
      const coordinates = series.map((relative, index) => relative === null ? null :
        [x(frequencies[index]).toFixed(2), y(relative === 0 ? minPower : Math.log10(relative)).toFixed(2)]);
      const trace = coordinates.map((point, index) => {
        if (point === null) return "";
        return `${coordinates[index - 1] ? "L" : "M"}${point.join(",")}`;
      }).join(" ");
      svg.push(`<path d="${trace}" fill="none" stroke="${color}" stroke-width="2.2"${dash}/>`);
      coordinates.forEach((point, index) => {
        if (point === null || (rows.length > 21 && series[index] !== 0)) return;
        const [cx, cy] = point;
        svg.push(marker === "circle"
          ? `<circle cx="${cx}" cy="${cy}" r="3.2" fill="white" stroke="${color}" stroke-width="1.7"/>`
          : `<rect x="${Number(cx) - 3}" y="${Number(cy) - 3}" width="6" height="6" fill="${color}"/>`);
      });
    });
    svg.push(`<text x="${left + panelWidth / 2}" y="${top + panelHeight + 44}" text-anchor="middle" font-size="13">Frequency [GHz]</text>`);
  });
  names.forEach((name, index) => {
    const x = 461 + 125 * index;
    const { color, dash } = styles[name];
    svg.push(`<line x1="${x}" y1="480" x2="${x + 35}" y2="480" stroke="${color}" stroke-width="2.2"${dash}/>`);
    svg.push(`<text x="${x + 43}" y="485" font-size="13">${name}</text>`);
  });
  const notes = [
    hasZeros ? "Exact zeros are shown at the plot floor." : "",
    allErrors.includes(null) ? "Relative error is undefined where the Direct magnitude is zero." : "",
  ].filter(Boolean).join(" ");
  if (notes) svg.push(`<text x="580" y="510" text-anchor="middle" font-size="12" fill="#465969">${notes}</text>`);
  svg.push("</g>", "</svg>");
  return `${svg.join("\n")}\n`;
}

export function compareAgainstDirect(resultsRoot, readCsv, writeCsv, runDirectories = {}, allowLegacy = false) {
  const load = (name) => readRun(runDirectories[name] ?? path.join(resultsRoot, `result_${name}`),
    name, readCsv, allowLegacy);
  const direct = load("DIRECT");
  const models = Object.keys(styles).map(load);
  const settings = modelSettings(direct);
  const verification = Object.fromEntries(models.map((model) => [model.name, assertComparable(direct, model)]));
  if (!direct.report.inputs && models.every((model) => model.report.inputs)) {
    compareModelSettings(models[0], models[1]);
  }
  const rows = direct.samples.map(({ frequencyHz }) => ({ freq_Hz: frequencyHz }));
  const frequencies = rows.map((row) => row.freq_Hz);
  const report = {
    modelSettings: settings,
    inputs: direct.report.inputs ?? null,
    inputVerification: verification,
    directDirectory: direct.directory,
    sampleCount: rows.length,
    frequencyRangeHz: [frequencies[0], frequencies.at(-1)],
    relativeErrorDefinition: "abs(S_ROM - S_Direct) / abs(S_Direct); null if abs(S_Direct) = 0",
    directSolver: direct.report.linear_solver_backend,
    models: Object.fromEntries(models.map(({ name, directory, report, diagnostics }) => [name, {
      directory,
      linearSolver: report.linear_solver_backend,
      expansionFrequenciesHz: diagnostics.expansion_frequencies_hz,
      requestedOrder: diagnostics.requested_order,
      localOrders: diagnostics.local_orders,
      romDimension: diagnostics.rom_dimension,
      stoppingReason: diagnostics.stopping_reason ?? null,
      trueResidualAcceptance: diagnostics.true_residual_acceptance ?? null,
      trueResidualTolerance: diagnostics.true_residual_tolerance ?? null,
      unmetFrequenciesHz: diagnostics.unmet_frequencies_hz ?? null,
      errors: {},
    }])),
    pairwiseDifferenceRelativeToDirect: {},
  };
  for (const channel of channels) {
    const referenceMagnitudes = direct.samples.map((row) => Math.hypot(...row[channel]));
    rows.forEach((row, index) => { row[`DIRECT_${channel}_magnitude`] = referenceMagnitudes[index]; });
    const addErrors = (name, absolute) => {
      const relative = absolute.map((value, index) => referenceMagnitudes[index] > 0
        ? value / referenceMagnitudes[index] : null);
      rows.forEach((row, index) => {
        row[`${name}_${channel}_absolute_complex_error`] = absolute[index];
        row[`${name}_${channel}_relative_error`] = relative[index];
      });
      return summarizeErrors(absolute, relative, referenceMagnitudes, frequencies);
    };
    for (const model of models) {
      report.models[model.name].errors[channel] = addErrors(model.name, model.samples.map(
        (sample, index) => complexDifference(sample[channel], direct.samples[index][channel])));
    }
    const [wcawe, mgawe] = models;
    report.pairwiseDifferenceRelativeToDirect[channel] = addErrors("MGAWE_minus_WCAWE", mgawe.samples.map(
      (sample, index) => complexDifference(sample[channel], wcawe.samples[index][channel])));
  }

  return publishComparison(resultsRoot, "result_COMPARE_DIRECT", (directory) => {
    writeCsv(path.join(directory, "direct_complex_errors.csv"), rows);
    fs.writeFileSync(path.join(directory, "direct_complex_error_report.json"), `${JSON.stringify(report, null, 2)}\n`);
    fs.writeFileSync(path.join(directory, "direct_complex_errors.svg"), drawComparison(rows, models.map(({ name }) => name)));
  });
}
