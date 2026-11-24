import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import { createHash } from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { parseArgs } from "node:util";
import { fileURLToPath } from "node:url";

const { values: args } = parseArgs({ options: {
  solver: { type: "string" },
  probe: { type: "string" },
  benchmark: { type: "string" },
  out: { type: "string" },
  "linear-solver": { type: "string", default: "direct" },
  threads: { type: "string", default: "2" },
} });
for (const key of ["solver", "probe", "benchmark", "out"]) {
  assert(args[key], "Required: --solver <exe> --probe <io_parser_tests.exe> --benchmark <root> --out <new-directory>");
}
assert(/^[1-9]\d*$/.test(args.threads), "--threads must be a positive integer");
const solver = path.resolve(args.solver);
const probe = path.resolve(args.probe);
const benchmark = path.resolve(args.benchmark);
const out = path.resolve(args.out);
const fixture = fileURLToPath(new URL("../tests/fixtures/rectangular-waveguide/", import.meta.url));
for (const executable of [solver, probe]) assert(fs.statSync(executable).isFile(), executable);
assert(!fs.existsSync(out), "Choose a new --out directory to preserve previous evidence");
fs.mkdirSync(out, { recursive: true });

const env = { ...process.env, MKL_NUM_THREADS: args.threads, OMP_NUM_THREADS: args.threads };
const settings = {
  basisOrder: 0,
  linearSolver: args["linear-solver"],
  threads: Number(args.threads),
  linearResidualTolerance: 1e-7,
  powerBalanceTolerance: 1e-6,
  legacyDualPowerBalanceTolerance: 5e-4,
  pillboxMinTransmissionMagnitude: 1e-2,
  waveguideCutoffRelativeTolerance: [0.01, 0.06],
  waveguideMaxReflectionMagnitude: 0.05,
  waveguideMaxTransmissionPhaseErrorRad: 0.10,
};
const summary = { schema: 1, status: "running", solver, probe, settings, cases: [] };
const writeJson = (name, object) => fs.writeFileSync(path.join(out, name), JSON.stringify(object, null, 2) + "\n");
const sha256 = file => createHash("sha256").update(fs.readFileSync(file)).digest("hex");
const c0 = 299792458;

// Keep the embedded binary sections intact while making the 5.6–6.0 GHz
// grid independent of the nominal sweep stored in the benchmark project.
function makePillboxInput(original) {
  const copy = path.join(out, "pillbox_5p6_6p0.aedt");
  let text = fs.readFileSync(original).toString("latin1");
  for (const [name, value] of [["RangeStart", "'5.6GHz'"], ["RangeEnd", "'6.0GHz'"], ["RangeCount", "9"]]) {
    const pattern = new RegExp("(^[\\t ]*" + name + "=)(?:'[^'\\r\\n]*'|\\d+)(?=[\\t \\r]*$)", "gm");
    assert.equal([...text.matchAll(pattern)].length, 1, original + ": expected one " + name);
    text = text.replace(pattern, (_, prefix) => prefix + value);
  }
  fs.writeFileSync(copy, Buffer.from(text, "latin1"));
  return copy;
}

function execute(binary, argv, name) {
  const result = spawnSync(binary, argv, { env, encoding: "utf8", maxBuffer: 16 * 1024 * 1024,
    timeout: 180_000, windowsHide: true });
  const log = path.join(out, name);
  fs.writeFileSync(log, (result.stdout ?? "") + "\n" + (result.stderr ?? ""));
  if (result.error) throw result.error;
  assert.equal(result.status, 0, log + ": unexpected exit code");
  return result.stdout;
}

function getRows(file) {
  const [header, ...lines] = fs.readFileSync(file, "utf8").trim().split(/\r?\n/);
  const columns = header.replace(/^\uFEFF/, "").split(",");
  return lines.map(line => {
    const values = line.split(",");
    assert.equal(values.length, columns.length, file + ": malformed row");
    return Object.fromEntries(columns.map((key, index) => {
      assert(values[index].trim().length > 0, file + ": empty " + key);
      const n = Number(values[index]);
      assert(Number.isFinite(n), file + ": nonfinite " + key);
      return [key, n];
    }));
  });
}

function runCase({ name, aedt, mesh, method, frequencyHz, analyticCutoffs,
  original = aedt }) {
  const powerTolerance = method === "tfe" ? settings.legacyDualPowerBalanceTolerance : settings.powerBalanceTolerance;
  const entry = { name, method, input: { aedt, mesh, originalAedt: original,
    aedtSha256: sha256(aedt), originalAedtSha256: sha256(original), meshSha256: sha256(mesh) },
    frequenciesHz: frequencyHz, powerTolerance, status: "running" };
  summary.cases.push(entry);
  try {
    const audit = JSON.parse(execute(probe, ["--port-baseline", aedt, mesh, method], name + "-assembly.log"));
    entry.assembly = audit;
    assert.deepEqual(audit.samples.map(s => s.frequency_hz), frequencyHz);
    const modesPerPort = method === "tfe" ? 2 : 1;
    assert.equal(audit.channels.length, 2 * modesPerPort);
    const physicalPortIds = [...new Set(audit.channels.map(c => c.port_id))];
    assert.equal(physicalPortIds.length, 2);
    const countExcited = audit.channels.filter(c => c.excited).length;
    assert.equal(countExcited, 1);
    const excited = audit.channels[0];
    assert(excited.excited && excited.mode === 1);
    for (const port of physicalPortIds) {
      const modes = audit.channels.filter(c => c.port_id === port);
      assert.deepEqual(modes.map(c => c.mode), Array.from({ length: modesPerPort }, (_, i) => i + 1));
      for (const channel of modes) {
        if (analyticCutoffs) {
          const expected = analyticCutoffs[channel.mode - 1];
          channel.analytic_cutoff_hz = expected;
          channel.analytic_cutoff_relative_error = Math.abs(channel.cutoff_hz / expected - 1);
          assert(channel.analytic_cutoff_relative_error <= settings.waveguideCutoffRelativeTolerance[channel.mode - 1],
            name + ": cutoff differs from analytic rectangular guide");
        }
      }
    }

    const runRoot = path.join(out, name);
    const cli = ["--aedt", aedt, "--mesh", mesh, "--sweep", "direct", "--basis-order", "0",
      "--port-method", method, "--tfe-modes-per-port", String(modesPerPort),
      "--linear-solver", settings.linearSolver, "--tolerance", String(settings.linearResidualTolerance),
      "--max-sweep-points", String(frequencyHz.length), "--no-write-all-fields", "--out", runRoot];
    execute(solver, cli, name + "-direct.log");
    const resultDir = path.join(runRoot, "result_DIRECT");
    const report = JSON.parse(fs.readFileSync(path.join(resultDir, "run.json"), "utf8"));
    const rows = getRows(path.join(resultDir, "s_parameters.csv"));
    assert.equal(report.status, "ok");
    entry.linearSolverBackend = report.linear_solver_backend;
    assert.deepEqual(rows.map(row => row.freq_Hz), frequencyHz);
    assert.deepEqual(report.sweep.map(row => row.freq_hz), frequencyHz);
    entry.direct = rows.map((row, i) => {
      const f = row.freq_Hz;
      const solve = report.sweep[i];
      assert.equal(solve.linear_solve_status, "converged");
      assert(Number.isFinite(solve.residual) && solve.residual <= settings.linearResidualTolerance);
      const observed = audit.channels.filter(c => c.mode === 1 && c.cutoff_hz < f);
      const propagating = audit.channels.filter(c => c.cutoff_hz < f);
      assert.equal(observed.length, 2, name + ": both main modes must propagate");
      assert.equal(propagating.length, observed.length, name + ": unobserved propagating retained modes");
      // Only first modes are observable in the legacy CSV. Here they form
      // the entire propagating set. This is one incident power-normalized column.
      const s11 = Math.hypot(row.S11_real, row.S11_imag);
      const s21 = Math.hypot(row.S21_real, row.S21_imag);
      const observedPower = s11 ** 2 + s21 ** 2;
      assert(Math.abs(observedPower - 1) <= powerTolerance, name + ": power balance");
      const sample = { frequencyHz: f, residual: solve.residual, s11Magnitude: s11, s21Magnitude: s21,
        observedPower, propagatingOutputs: observed.map(c => [c.port_id, c.mode]) };
      if (analyticCutoffs) {
        assert(s11 <= settings.waveguideMaxReflectionMagnitude, name + ": uniform guide reflects too strongly");
        const beta = 2 * Math.PI / c0 * Math.sqrt(f ** 2 - analyticCutoffs[0] ** 2);
        const phaseDifference = Math.atan2(row.S21_imag, row.S21_real) + beta * 0.02;
        // Independently solved real eigenvectors can differ by a sign.
        sample.phaseErrorModPiRad = Math.abs(phaseDifference - Math.round(phaseDifference / Math.PI) * Math.PI);
        assert(sample.phaseErrorModPiRad <= settings.waveguideMaxTransmissionPhaseErrorRad,
          name + ": wrong uniform-guide propagation phase");
      } else {
        assert(s21 >= settings.pillboxMinTransmissionMagnitude, name + ": transmission blocks at dielectric");
      }
      return sample;
    });
    entry.status = "pass";
  } catch (error) {
    entry.status = "fail";
    entry.error = error.message;
    throw error;
  }
}

try {
  execute(probe, [], "io-parser-regression.log");
  summary.ioParserRegression = "pass";
  const waveguideAedt = path.resolve(fixture, "waveguide.aedt");
  const waveguideMesh = path.resolve(fixture, "waveguide.ngmesh");
  const waveguideFrequencies = [9e9, 9.25e9, 9.5e9, 9.75e9, 10e9];
  const cutoffs = [c0 / (2 * 0.020), c0 / (2 * 0.012)]; // TE10, TE01
  const waveguide = { aedt: waveguideAedt, mesh: waveguideMesh,
    frequencyHz: waveguideFrequencies, analyticCutoffs: cutoffs };
  runCase({ ...waveguide, name: "waveguide-single", method: "numerical" });
  runCase({ ...waveguide, name: "waveguide-dual", method: "tfe" });
  const pillbox = path.join(benchmark, "pillbox_window", "project");
  const original = path.join(pillbox, "Project4.aedt");
  const window = { aedt: makePillboxInput(original), mesh: path.join(pillbox, "current.ngmesh"),
    frequencyHz: Array.from({ length: 9 }, (_, i) => 5_600_000_000 + i * 50_000_000), original };
  runCase({ ...window, name: "pillbox-single", method: "numerical" });
  runCase({ ...window, name: "pillbox-dual", method: "tfe" });
  summary.status = "pass";
} catch (error) {
  summary.status = "fail";
  summary.error = error.message;
  console.error(error.message);
  process.exitCode = 1;
} finally {
  writeJson("port-assembly-smoke.json", summary);
  console.log(summary.status + ": " + path.join(out, "port-assembly-smoke.json"));
}
