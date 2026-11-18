#include "bpfem/app/Application.hpp"

#include "bpfem/bc/WavePortBC.hpp"
#include "bpfem/core/Logger.hpp"
#include "bpfem/core/RunEnvironment.hpp"
#include "bpfem/core/RunReport.hpp"
#include "bpfem/core/Utilities.hpp"
#include "bpfem/factory/SparseSolverFactory.hpp"
#include "bpfem/factory/SweepStrategyFactory.hpp"
#include "bpfem/fem/EdgeTopology.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/io/AEDTParser.hpp"
#include "bpfem/io/NGMeshParser.hpp"
#include "bpfem/io/PortFaceResolver.hpp"
#include "bpfem/linalg/IFactorizedSparseSolver.hpp"
#include "bpfem/linalg/ISparseSolver.hpp"
#include "bpfem/linalg/SolveResidual.hpp"
#include "bpfem/post/OutputWriter.hpp"
#include "bpfem/post/ResultExtractor.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"
#include "bpfem/apm/AnalyticPortBuilder.hpp"
#include "bpfem/tfe/TransfiniteElementBuilder.hpp"
#include "bpfem/tfe/TfeSystem.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem {

namespace {

double bytesToMb(std::size_t bytes) {
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

std::vector<PhaseSummary> collectPhaseSummaries(const Logger& log) {
    std::vector<PhaseSummary> out;
    const auto snapshots = log.phaseSnapshots();
    out.reserve(snapshots.size());
    for (const auto& p : snapshots) {
        PhaseSummary s;
        s.name = p.name;
        s.elapsedSec = p.elapsedSec;
        s.endCurrentMb = p.endMem.currentBytes > 0 ? bytesToMb(p.endMem.currentBytes) : 0.0;
        s.endPeakMb = p.endMem.peakBytes > 0 ? bytesToMb(p.endMem.peakBytes) : 0.0;
        out.push_back(std::move(s));
    }
    return out;
}

std::string sweepOutputName(SweepStrategy strategy) {
    switch (strategy) {
        case SweepStrategy::Direct: return "DIRECT";
        case SweepStrategy::Alps:   return "ALPS";
        case SweepStrategy::Awe:    return "AWE";
        case SweepStrategy::Gawe:   return "GAWE";
        case SweepStrategy::Mgawe:  return "MGAWE";
        case SweepStrategy::Wcawe:  return "WCAWE";
    }
    return "UNKNOWN";
}

std::filesystem::path runOutputDirectory(const Options& options) {
    if (options.tfeProxyDirect) return options.outDir / "result_TFE_PROXY_DIRECT";
    return options.outDir / ("result_" + sweepOutputName(options.sweepStrategy));
}

void prepareCleanOutputDirectory(const std::filesystem::path& path) {
    if (path.empty()) {
        throw std::runtime_error("Output directory must not be empty");
    }
    if (std::filesystem::exists(path)) {
        std::filesystem::remove_all(path);
    }
    std::filesystem::create_directories(path);
}

void configurePardisoOutOfCore(const std::filesystem::path& outputDirectory) {
    const auto scratch = std::filesystem::absolute(outputDirectory / "pardiso-ooc");
    std::filesystem::create_directories(scratch);
    const auto prefix = (scratch / "factors").string();
#ifdef _WIN32
    std::size_t existingLength = 0;
    const int query = getenv_s(&existingLength, nullptr, 0, "MKL_PARDISO_OOC_FILE_NAME");
    const int status = query != 0 ? query : existingLength != 0 ? 0
        : _putenv_s("MKL_PARDISO_OOC_FILE_NAME", prefix.c_str());
#else
    const int status = setenv("MKL_PARDISO_OOC_FILE_NAME", prefix.c_str(), 0);
#endif
    if (status != 0) {
        throw std::runtime_error("Could not configure PARDISO out-of-core files");
    }
}

}  // namespace

// Parse the bp_fem_solver command line into an Options struct. Each option
// takes either zero or one value (`--flag` vs `--flag <value>`), with the
// value consumed via the requireValue lambda to give a precise error if the
// user dropped the value. Unknown options throw rather than silently
// continue. The full list of supported flags + valid values is mirrored in
// the --help text emitted at the end of the parser.
Options parseOptions(int argc, char** argv) {
    Options options;
    bool alpsOrderSpecified = false;
    bool alpsLegacyOrderSpecified = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto requireValue = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for " + name);
            }
            return argv[++i];
        };
        if (arg == "--aedt") {
            options.aedt = requireValue(arg);
        } else if (arg == "--mesh") {
            options.mesh = requireValue(arg);
        } else if (arg == "--out") {
            options.outDir = requireValue(arg);
        } else if (arg == "--max-sweep-points") {
            options.maxSweepPoints = std::stoi(requireValue(arg));
        } else if (arg == "--max-iterations") {
            options.maxIterations = std::stoi(requireValue(arg));
        } else if (arg == "--basis-order") {
            options.basisOrder = std::stoi(requireValue(arg));
            if (options.basisOrder < 0 || options.basisOrder > 1) {
                throw std::runtime_error("--basis-order must be 0 or 1");
            }
        } else if (arg == "--field-output-order") {
            options.fieldOutputOrder = std::stoi(requireValue(arg));
            if (options.fieldOutputOrder < 1 || options.fieldOutputOrder > 3) {
                throw std::runtime_error("--field-output-order must be in [1, 3]");
            }
        } else if (arg == "--tolerance") {
            options.tolerance = std::stod(requireValue(arg));
        } else if (arg == "--write-all-fields") {
            options.writeAllFields = true;
        } else if (arg == "--no-write-all-fields") {
            options.writeAllFields = false;
        } else if (arg == "--sweep") {
            const std::string mode = requireValue(arg);
            if (mode == "direct") {
                options.sweepStrategy = SweepStrategy::Direct;
            } else if (mode == "alps") {
                options.sweepStrategy = SweepStrategy::Alps;
            } else if (mode == "awe") {
                options.sweepStrategy = SweepStrategy::Awe;
            } else if (mode == "gawe") {
                options.sweepStrategy = SweepStrategy::Gawe;
            } else if (mode == "mgawe") {
                options.sweepStrategy = SweepStrategy::Mgawe;
            } else if (mode == "wcawe") {
                options.sweepStrategy = SweepStrategy::Wcawe;
            } else {
                throw std::runtime_error("--sweep must be 'direct', 'alps', 'awe', 'gawe', 'mgawe', or 'wcawe'");
            }
        } else if (arg == "--alps-order") {
            if (alpsLegacyOrderSpecified) {
                throw std::runtime_error(
                    "--alps-order and --alps-krylov-order cannot be specified together");
            }
            options.alpsOrder = std::stoi(requireValue(arg));
            alpsOrderSpecified = true;
            if (options.alpsOrder < 1) {
                throw std::runtime_error("--alps-order must be >= 1");
            }
        } else if (arg == "--alps-krylov-order") {
            if (alpsOrderSpecified) {
                throw std::runtime_error(
                    "--alps-order and --alps-krylov-order cannot be specified together");
            }
            options.alpsOrder = std::stoi(requireValue(arg));
            alpsLegacyOrderSpecified = true;
            if (options.alpsOrder < 1) {
                throw std::runtime_error("--alps-krylov-order must be >= 1");
            }
        } else if (arg == "--alps-expansion") {
            options.alpsExpansionFrequencyHz = std::stod(requireValue(arg));
            if (options.alpsExpansionFrequencyHz <= 0.0) {
                throw std::runtime_error("--alps-expansion must be > 0 Hz");
            }
        } else if (arg == "--awe-field-pade") {
            options.aweFieldPade = true;
        } else if (arg == "--awe-order") {
            options.aweOrder = std::stoi(requireValue(arg));
            if (options.aweOrder < 1) {
                throw std::runtime_error("--awe-order must be >= 1");
            }
        } else if (arg == "--awe-expansion") {
            options.aweExpansionFrequencyHz = std::stod(requireValue(arg));
            if (options.aweExpansionFrequencyHz <= 0.0) {
                throw std::runtime_error("--awe-expansion must be > 0 Hz");
            }
        } else if (arg == "--gawe-order") {
            options.gaweOrder = std::stoi(requireValue(arg));
            if (options.gaweOrder < 1) {
                throw std::runtime_error("--gawe-order must be >= 1");
            }
        } else if (arg == "--gawe-expansion") {
            options.gaweExpansionFrequencyHz = std::stod(requireValue(arg));
            if (options.gaweExpansionFrequencyHz <= 0.0) {
                throw std::runtime_error("--gawe-expansion must be > 0 Hz");
            }
        } else if (arg == "--gawe-drop-tolerance") {
            options.gaweDropTolerance = std::stod(requireValue(arg));
            if (options.gaweDropTolerance <= 0.0) {
                throw std::runtime_error("--gawe-drop-tolerance must be > 0");
            }
        } else if (arg == "--mgawe-expansion") {
            const double hz = std::stod(requireValue(arg));
            const int order = std::stoi(requireValue(arg));
            if (!std::isfinite(hz) || hz <= 0.0 || order < 1) {
                throw std::runtime_error("--mgawe-expansion requires a positive frequency in Hz and order");
            }
            options.mgaweExpansions.emplace_back(hz, order);
        } else if (arg == "--mgawe-2001") {
            options.mgawe2001 = true;
        } else if (arg == "--mgawe-adaptive") {
            options.mgaweAdaptive = true;
        } else if (arg == "--mgawe-points") {
            options.mgawePointCount = std::stoi(requireValue(arg));
            if (options.mgawePointCount < 1) {
                throw std::runtime_error("--mgawe-points must be >= 1");
            }
        } else if (arg == "--mgawe-order") {
            options.mgaweOrder = std::stoi(requireValue(arg));
            if (options.mgaweOrder < 1) {
                throw std::runtime_error("--mgawe-order must be >= 1");
            }
        } else if (arg == "--mgawe-drop-tolerance") {
            options.mgaweDropTolerance = std::stod(requireValue(arg));
            if (options.mgaweDropTolerance <= 0.0) {
                throw std::runtime_error("--mgawe-drop-tolerance must be > 0");
            }
        } else if (arg == "--wcawe-check-moments") {
            options.wcaweCheckMoments = true;
        } else if (arg == "--wcawe-order") {
            options.wcaweOrder = std::stoi(requireValue(arg));
            if (options.wcaweOrder < 1) {
                throw std::runtime_error("--wcawe-order must be >= 1");
            }
        } else if (arg == "--wcawe-expansion") {
            options.wcaweExpansionFrequencyHz = std::stod(requireValue(arg));
            if (options.wcaweExpansionFrequencyHz <= 0.0) {
                throw std::runtime_error("--wcawe-expansion must be > 0 Hz");
            }
        } else if (arg == "--wcawe-drop-tolerance") {
            options.wcaweDropTolerance = std::stod(requireValue(arg));
            if (options.wcaweDropTolerance <= 0.0) {
                throw std::runtime_error("--wcawe-drop-tolerance must be > 0");
            }
        } else if (arg == "--port-method") {
            const std::string m = requireValue(arg);
            if (m == "numerical" || m == "npm") {
                options.portMethod = PortMethod::Numerical;
            } else if (m == "analytic") {
                options.portMethod = PortMethod::Analytic;
            } else if (m == "tfe" || m == "transfinite") {
                options.portMethod = PortMethod::Transfinite;
            } else if (m == "tfe-block") {
                options.portMethod = PortMethod::TfeBlock;
            } else {
                throw std::runtime_error("--port-method must be numerical, analytic, tfe, or tfe-block");
            }
        } else if (arg == "--tfe-modes-per-port") {
            options.tfeModesSpecified = true;
            options.tfeModesPerPort = std::stoi(requireValue(arg));
            if (options.tfeModesPerPort < 1) {
                throw std::runtime_error("--tfe-modes-per-port must be >= 1");
            }
        } else if (arg == "--tfe-all-sources") {
            options.tfeAllSources = true;
        } else if (arg == "--tfe-proxy-direct") {
            options.tfeProxyDirect = true;
        } else if (arg == "--tfe-subband-width") {
            options.tfeSubbandWidthHz = std::stod(requireValue(arg));
            if (!std::isfinite(options.tfeSubbandWidthHz) || options.tfeSubbandWidthHz <= 0.0) {
                throw std::runtime_error("--tfe-subband-width requires a positive width in Hz");
            }
        } else if (arg == "--tfe-exact-band") {
            const double left = std::stod(requireValue(arg));
            const double right = std::stod(requireValue(arg));
            if (!std::isfinite(left) || !std::isfinite(right) || left <= 0.0 || right < left) {
                throw std::runtime_error("--tfe-exact-band requires two ordered positive frequencies in Hz");
            }
            options.tfeExactBandsHz.emplace_back(left, right);
        } else if (arg == "--tfe-source") {
            const auto source = requireValue(arg);
            const auto separator = source.find(':');
            if (separator == std::string::npos) throw std::runtime_error("--tfe-source requires PORT:MODE");
            options.tfeSourcePort = std::stoi(source.substr(0, separator));
            options.tfeSourceMode = std::stoi(source.substr(separator + 1));
            if (options.tfeSourcePort < 0 || options.tfeSourceMode < 1) {
                throw std::runtime_error("--tfe-source requires a nonnegative port ID and one-based mode");
            }
        } else if (arg == "--linear-solver") {
            const std::string m = requireValue(arg);
            if (m == "auto") {
                options.linearSolver = LinearSolverKind::Auto;
            } else if (m == "direct") {
                options.linearSolver = LinearSolverKind::Direct;
            } else if (m == "bicgstab") {
                options.linearSolver = LinearSolverKind::BiCGStab;
            } else if (m == "gmres") {
                options.linearSolver = LinearSolverKind::Gmres;
            } else if (m == "mumps") {
                options.linearSolver = LinearSolverKind::Mumps;
            } else {
                throw std::runtime_error("--linear-solver must be auto, direct, bicgstab, gmres or mumps");
            }
        } else if (arg == "--pardiso-mode") {
            const auto mode = requireValue(arg);
            if (mode != "incore" && mode != "ooc") {
                throw std::runtime_error("--pardiso-mode must be 'incore' or 'ooc'");
            }
            options.pardisoOutOfCore = mode == "ooc";
        } else if (arg == "--mumps-mode") {
            options.mumpsMode = requireValue(arg);
            if (options.mumpsMode != "incore" && options.mumpsMode != "ooc"
                && options.mumpsMode != "blr") {
                throw std::runtime_error("--mumps-mode must be incore, ooc or blr");
            }
        } else if (arg == "--mumps-blr-tolerance") {
            options.mumpsBlrTolerance = std::stod(requireValue(arg));
            if (!(options.mumpsBlrTolerance > 0.0 && options.mumpsBlrTolerance < 1.0)) {
                throw std::runtime_error("--mumps-blr-tolerance must be between 0 and 1");
            }
        } else if (arg == "--schwarz-local-solver") {
            const auto value = requireValue(arg);
            if (value != "pardiso" && value != "mumps") {
                throw std::runtime_error("--schwarz-local-solver must be pardiso or mumps");
            }
            options.schwarzLocalMumps = value == "mumps";
        } else if (arg == "--schwarz-subdomains") {
            options.schwarzSubdomains = std::stoi(requireValue(arg));
            if (options.schwarzSubdomains < 1) {
                throw std::runtime_error("--schwarz-subdomains must be >= 1");
            }
        } else if (arg == "--schwarz-overlap") {
            options.schwarzOverlap = std::stoi(requireValue(arg));
            if (options.schwarzOverlap < 0) {
                throw std::runtime_error("--schwarz-overlap must be >= 0");
            }
        } else if (arg == "--precon") {
            const std::string m = requireValue(arg);
            if (m == "none") {
                options.preconditioner = PreconditionerKind::None;
            } else if (m == "jacobi") {
                options.preconditioner = PreconditionerKind::Jacobi;
            } else if (m == "ilu0") {
                options.preconditioner = PreconditionerKind::Ilu0;
            } else if (m == "schwarz") {
                options.preconditioner = PreconditionerKind::Schwarz;
            } else if (m == "pardiso-single") {
                options.preconditioner = PreconditionerKind::PardisoSingle;
            } else {
                throw std::runtime_error("--precon must be none, jacobi, ilu0, schwarz or pardiso-single");
            }
        } else if (arg == "--iterative-spmv") {
            const auto value = requireValue(arg);
            if (value != "native" && value != "mkl") {
                throw std::runtime_error("--iterative-spmv must be native or mkl");
            }
            options.iterativeMklSpmv = value == "mkl";
        } else if (arg == "--gmres-restart") {
            options.gmresRestart = std::stoi(requireValue(arg));
            if (options.gmresRestart < 1) {
                throw std::runtime_error("--gmres-restart must be >= 1");
            }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: bp_fem_solver [--aedt wg_bp_filter.aedt] [--mesh current.ngmesh] [--out result]\n"
                      << "                     [--max-sweep-points 21] [--max-iterations 400] [--tolerance 1e-7]\n"
                      << "                     [--basis-order 0|1] [--field-output-order 1|2|3]\n"
                      << "                     [--write-all-fields|--no-write-all-fields]\n"
                      << "                     [--sweep direct|alps|awe|gawe|mgawe|wcawe]\n"
                      << "                     [--alps-order 12] [--alps-krylov-order <deprecated-alias>] [--alps-expansion <Hz>]\n"
                      << "                     [--awe-order 8] [--awe-expansion <Hz>] [--awe-field-pade]\n"
                      << "                     [--gawe-order 12] [--gawe-expansion <Hz>] [--gawe-drop-tolerance 1e-10]\n"
                      << "                     [--mgawe-points 3] [--mgawe-order 8] [--mgawe-drop-tolerance 1e-10] [--mgawe-adaptive|--mgawe-2001]\n"
                      << "                     [--mgawe-expansion <Hz> <order>] (repeatable)\n"
                      << "                     [--wcawe-order 12] [--wcawe-expansion <Hz>] [--wcawe-drop-tolerance 1e-12] [--wcawe-check-moments]\n"
                      << "                     [--port-method numerical|analytic|tfe|tfe-block] [--tfe-modes-per-port N]\n"
                      << "                     [--tfe-all-sources | --tfe-source PORT:MODE] (tfe-block)\n"
                      << "                     [--tfe-proxy-direct] [--tfe-subband-width HZ] (tfe-block with --sweep mgawe)\n"
                      << "                     [--tfe-exact-band LOW_HZ HIGH_HZ] (repeatable, tfe-block MGAWE)\n"
                      << "                     [--linear-solver auto|direct|bicgstab|gmres|mumps]\n"
                      << "                     [--precon none|jacobi|ilu0|schwarz|pardiso-single]\n"
                      << "                     [--gmres-restart 30]\n"
                      << "                     [--iterative-spmv native|mkl]\n"
                      << "                     [--pardiso-mode incore|ooc] [--mumps-mode incore|ooc|blr]\n"
                      << "                     [--mumps-blr-tolerance 1e-10] [--schwarz-subdomains 8]\n"
                      << "                     [--schwarz-overlap 1] [--schwarz-local-solver pardiso|mumps]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown option: " + arg);
        }
    }
    return options;
}

// Top-level driver:
//   1. Parse Options.
//   2. Load AEDT project + NGMesh.
//   3. Build EdgeTopology + PortModeSolver, populate per-port modes via the
//      requested port-method strategy (analytic APM or numerical multi-mode TFE).
//   4. Construct FEMAssembler and register the boundary-condition list
//      (currently always [WavePortBC]; future BCs are appended here).
//   5. Pull a sparse-solver backend and a sweep strategy from the factories.
//   6. Run the sweep strategy with a SweepContext; the strategy decides
//      whether per-frequency VTU output is meaningful and either invokes
//      the onFieldSolved callback (DirectSweep) or ignores it (AlpsSweep).
//   7. Write s_parameters.csv and field_last.vtu.
//
// All exceptions are caught here; on failure the function returns 1 and the
// process exits with that code via main.cpp.
int runApplication(int argc, char** argv) {
    // The Logger lifetime spans the entire run so its destructor can flush
    // a final summary even if we return early. Both `crashContext` and
    // `report` are populated as the run progresses; on a fatal exception we
    // dump them to run.crash.log + run.json("status":"fatal") before
    // returning non-zero.
    Logger log(std::cout);
    RunReport report;
    std::vector<FrequencyStat> sweepStats;
    std::filesystem::path crashLogPath;
    std::filesystem::path runJsonPath;
    std::filesystem::path timingJsonPath;
    const auto runStart = std::chrono::steady_clock::now();

    try {
        log.banner();

        Options options = parseOptions(argc, argv);
        const std::filesystem::path outputRoot = options.outDir;
        options.outDir = runOutputDirectory(options);
        prepareCleanOutputDirectory(options.outDir);
        if (options.pardisoOutOfCore) {
            configurePardisoOutOfCore(options.outDir);
        }

        // Attach the tee log file before any further info() so every line is
        // captured in run.log.
        log.attachLogFile((options.outDir / "run.log").string());
        crashLogPath = options.outDir / "run.crash.log";
        runJsonPath = options.outDir / "run.json";
        timingJsonPath = options.outDir / "timing.json";
        log.info("Output root: " + outputRoot.string());
        log.info("Output directory: " + options.outDir.string());

        // -------------------- Env header --------------------
        // Captured once and emitted as [env] tagged lines so the log is
        // self-describing: which solver version, which build flags, which
        // host, which CPU, total/avail memory, full command line, cwd.
        const RunEnvironment env = captureRunEnvironment(argc, argv);
        writeEnvHeader(log, env);
        fillReportFromEnvironment(report, env, options);

        // -------------------- Project import --------------------
        log.phase("Project import");
        log.info("Reading AEDT project: " + options.aedt.string());
        ProjectDefinition project = AEDTParser{}.parse(options.aedt);
        log.info("Materials: " + std::to_string(project.materials.size()));
        log.info("Wave ports: " + std::to_string(project.ports.size()));
        log.info("Sweep: " + std::to_string(project.sweep.startHz / 1.0e9) + " GHz to " + std::to_string(project.sweep.endHz / 1.0e9) + " GHz, " + std::to_string(project.sweep.count) + " points");

        log.info("Reading NGMesh: " + options.mesh.string());
        const Mesh mesh = NGMeshParser{}.parse(options.mesh);
        log.info("Mesh points: " + std::to_string(mesh.pointsById.size()));
        log.info("Surface triangles: " + std::to_string(mesh.surfaceTriangles.size()));
        log.info("Tetrahedra: " + std::to_string(mesh.tetrahedra.size()));
        PortFaceResolver::resolve(project, mesh);
        for (const auto& port : project.ports) {
            const std::string source = port.objectId >= 0
                ? "sheet object " + std::to_string(port.objectId)
                : "direct face";
            log.info("Wave port " + std::to_string(port.id)
                     + ": " + source + " -> mesh face " + std::to_string(port.faceId));
        }

        // -------------------- Frequency-domain solve setup --------------------
        log.phase("Frequency-domain solve");
        EdgeTopology topology(mesh, options.basisOrder);
        const int effectiveFieldOrder = options.fieldOutputOrder > 0
            ? options.fieldOutputOrder
            : std::max(1, options.basisOrder + 1);
        log.info("Nedelec basis order: " + std::to_string(options.basisOrder));
        log.info("Field output order: " + std::to_string(effectiveFieldOrder));
        log.info("Nedelec unknowns: " + std::to_string(topology.edgeCount()));
        log.info("Geometric edges: " + std::to_string(topology.geometricEdgeCount()));
        log.info("Local basis functions per tetrahedron: " + std::to_string(topology.localDofCount()));
        PortModeSolver portModeSolver(mesh, topology);
        const bool nativeTfe = options.portMethod == PortMethod::TfeBlock;
        if (options.tfeProxyDirect && (!nativeTfe || options.sweepStrategy != SweepStrategy::Mgawe)) {
            throw std::runtime_error("--tfe-proxy-direct requires tfe-block with --sweep mgawe");
        }
        if (options.tfeSubbandWidthHz > 0.0 && (!nativeTfe || options.sweepStrategy != SweepStrategy::Mgawe)) {
            throw std::runtime_error("--tfe-subband-width requires tfe-block with --sweep mgawe");
        }
        if (!options.tfeExactBandsHz.empty() &&
            (!nativeTfe || options.sweepStrategy != SweepStrategy::Mgawe || options.tfeProxyDirect)) {
            throw std::runtime_error("--tfe-exact-band requires tfe-block with --sweep mgawe");
        }
        if ((options.tfeAllSources || options.tfeSourcePort >= 0) &&
            (!nativeTfe || (options.tfeAllSources && options.tfeSourcePort >= 0))) {
            throw std::runtime_error("Select one source option with --port-method tfe-block");
        }
        if (nativeTfe && options.sweepStrategy != SweepStrategy::Direct && options.sweepStrategy != SweepStrategy::Mgawe) {
            throw std::runtime_error("tfe-block supports --sweep direct or mgawe");
        }
        if (nativeTfe && options.sweepStrategy == SweepStrategy::Mgawe &&
            (options.mgaweAdaptive || options.mgawe2001)) {
            throw std::runtime_error("tfe-block MGAWE uses fixed-point Chebyshev subbands");
        }
        std::string portMethodName;
        switch (options.portMethod) {
            case PortMethod::TfeBlock:
                portMethodName = "TFE block (modal trace replacement)";
                break;
            case PortMethod::Transfinite:
                portMethodName = "transfinite (TFE numerical eigenmodes, "
                                 + std::to_string(options.tfeModesPerPort) + " modes/port)";
                break;
            case PortMethod::Analytic:
                portMethodName = "analytic (TE10 closed-form, rectangular only)";
                break;
            case PortMethod::Numerical:
            default:
                portMethodName = "numerical (2D H(curl) port eigensolve, dominant mode)";
                break;
        }
        log.info("Port-mode method: " + portMethodName);
        if (options.portMethod == PortMethod::Numerical) {
            // Solve K_port v = k_c^2 M_port v on each port face. The first
            // call to PortModeSolver::solve(faceId) lazily triggers
            // computeMode -> computeMultiMode(faceId, 1), so we don't need
            // an explicit builder here. Doing it eagerly keeps timing
            // predictable and surfaces eigensolve failures up-front.
            apm::AnalyticPortBuilder rectDetector(mesh, topology);
            for (const auto& port : project.ports) {
                const auto geo = rectDetector.detectRectangular(port.faceId);
                if (geo.valid) {
                    log.info("Port face " + std::to_string(port.faceId)
                             + " rectangular dimensions (informational): a="
                             + std::to_string(geo.a * 1.0e3) + " mm, b="
                             + std::to_string(geo.b * 1.0e3) + " mm");
                } else {
                    log.info("Port face " + std::to_string(port.faceId)
                             + " non-rectangular cross section; numerical eigensolve still applies");
                }
                const auto& mode = portModeSolver.solve(port.faceId);
                log.info("Port face " + std::to_string(port.faceId) + " numerical mode: k_c^2="
                         + std::to_string(mode.cutoffWavenumberSquared)
                         + " (k_c=" + std::to_string(std::sqrt(mode.cutoffWavenumberSquared))
                         + " rad/m)");
            }
        } else if (options.portMethod == PortMethod::Analytic) {
            apm::AnalyticPortBuilder apmBuilder(mesh, topology);
            for (const auto& port : project.ports) {
                const auto geo = apmBuilder.detectRectangular(port.faceId);
                if (!geo.valid) {
                    throw std::runtime_error("--port-method analytic: face " + std::to_string(port.faceId) + " is not a recognizable rectangular port");
                }
                log.info("Port face " + std::to_string(port.faceId) + " rectangular dimensions: a="
                         + std::to_string(geo.a * 1.0e3) + " mm, b="
                         + std::to_string(geo.b * 1.0e3) + " mm");
                portModeSolver.setPrecomputed(port.faceId, apmBuilder.build(port.faceId));
            }
        } else if (options.portMethod == PortMethod::Transfinite) {
            tfe::TfeOptions tfeOpts;
            tfeOpts.modesPerPort = options.tfeModesPerPort;
            tfe::TransfiniteElementBuilder tfeBuilder(mesh, topology, portModeSolver, tfeOpts);
            for (const auto& port : project.ports) {
                const auto geo = tfeBuilder.detectRectangular(port.faceId);
                if (geo.valid) {
                    log.info("Port face " + std::to_string(port.faceId) + " rectangular dimensions (informational): a="
                             + std::to_string(geo.a * 1.0e3) + " mm, b="
                             + std::to_string(geo.b * 1.0e3) + " mm");
                } else {
                    log.info("Port face " + std::to_string(port.faceId)
                             + " rectangular detection skipped; using numerical eigenmodes only");
                }
                MultiPortMode multi = tfeBuilder.build(port.faceId);
                log.info("Port face " + std::to_string(port.faceId) + " TFE retained modes: "
                         + std::to_string(multi.modes.size()));
                for (std::size_t k = 0; k < multi.modes.size(); ++k) {
                    log.info("  mode " + std::to_string(k) + ": k_c^2="
                             + std::to_string(multi.modes[k].cutoffWavenumberSquared)
                             + " (k_c=" + std::to_string(std::sqrt(multi.modes[k].cutoffWavenumberSquared))
                             + " rad/m)");
                }
                portModeSolver.setMultiMode(port.faceId, std::move(multi));
            }
        }
        FEMAssembler assembler(mesh, project, topology, portModeSolver);
        std::unique_ptr<tfe::TfeSystem> tfeSystem;
        if (nativeTfe) {
            assembler.enableRadiationFaces();
            tfeSystem = std::make_unique<tfe::TfeSystem>(mesh, project, topology,
                assembler, portModeSolver, options.tfeModesSpecified ? options.tfeModesPerPort : 0);
            log.info("TFE interior DOFs " + std::to_string(tfeSystem->interiorCount())
                     + ", replaced trace DOFs " + std::to_string(tfeSystem->traceCount())
                     + ", retained modes " + std::to_string(tfeSystem->channels().size())
                     + ", upper NNZ " + std::to_string(tfeSystem->pattern()->columnIndices.size()));
            for (const auto& channel : tfeSystem->channels()) {
                log.info("TFE port " + std::to_string(channel.portId) + " mode " + std::to_string(channel.mode)
                         + (channel.transverseMagnetic ? " TM" : " TE")
                         + ", cutoff " + std::to_string(channel.cutoffHz) + " Hz");
            }
        } else {
            for (const auto& port : project.ports) {
                const auto& mode = portModeSolver.solve(port.faceId);
                log.info("Port face " + std::to_string(port.faceId) + " mode DOFs: " + std::to_string(mode.edgeDofs.size()));
            }
            std::vector<std::shared_ptr<bc::IBoundaryCondition>> bcs;
            bcs.push_back(std::make_shared<bc::WavePortBC>());
            assembler.setBoundaryConditions(std::move(bcs));
        }
        ResultExtractor extractor(project, portModeSolver);
        std::unique_ptr<linalg::ISparseSolver> sparseSolver = factory::makeSparseSolver(options);
        report.linearSolverBackend = sparseSolver->name();
        log.info("Linear solver backend: " + std::string(sparseSolver->name()));
        if (options.preconditioner != PreconditionerKind::None) {
            std::string preconName = "?";
            switch (options.preconditioner) {
                case PreconditionerKind::Jacobi: preconName = "Jacobi"; break;
                case PreconditionerKind::Ilu0:   preconName = "ILU(0)"; break;
                case PreconditionerKind::Schwarz: preconName = "Schwarz"; break;
                case PreconditionerKind::PardisoSingle: preconName = "PARDISO single"; break;
                case PreconditionerKind::None:   break;
            }
            log.info("Preconditioner: " + preconName);
        }
        auto sweepStrategy = factory::makeSweepStrategy(options, project, assembler, portModeSolver,
                                                         tfeSystem.get());
        report.sweepStrategyName = sweepStrategy->name();
        log.info("Sweep strategy: " + std::string(sweepStrategy->name()));

        const auto frequencies = buildFrequencies(project.sweep, options.maxSweepPoints);

        // ---- Per-frequency stat capture ----
        // The DirectSweep strategy invokes onFrequencyStat with assemble/solve
        // timing; we (a) emit a [stat] log row in the same format as the
        // bulk S-parameter table, and (b) push a FrequencyStat into the JSON
        // sweep array. Field S parameters are filled in below after the
        // sweep returns, when sparams is in scope.
        sweepStats.reserve(frequencies.size());

        sweep::SweepContext sweepCtx{
            project, assembler, portModeSolver, extractor, *sparseSolver, log, options.outDir
        };
        sweepCtx.linearMaxIterations = options.maxIterations;
        sweepCtx.linearTolerance = options.tolerance;
        if (options.writeAllFields) {
            sweepCtx.onFieldSolved = [&](double f, const std::vector<std::complex<double>>& edgeDofs) {
                const std::filesystem::path fieldPath =
                    options.outDir / ("field_" + std::to_string(static_cast<std::uint64_t>(std::llround(f))) + "Hz.vtu");
                OutputWriter::writeVTU(fieldPath, mesh, project, topology, edgeDofs, f, effectiveFieldOrder);
            };
        }
        sweepCtx.onFrequencyStat = [&](double f, double assembleSec, double solveSec,
                                        const SolveResult& result) {
            FrequencyStat fs;
            fs.frequencyHz = f;
            fs.assembleSec = assembleSec;
            fs.solveSec = solveSec;
            fs.iterations = result.iterations;
            fs.residual = result.residual;
            fs.absoluteResidual = result.absoluteResidual;
            fs.linearSolveStatus = linalg::solveStatusName(result.status);
            sweepStats.push_back(fs);
        };

        sweep::SweepResult swept = sweepStrategy->run(frequencies, sweepCtx);
        report.offlineBuildSec = swept.offlineBuildSec;
        report.portLinearizationSec = swept.portLinearizationSec;
        report.lanczosOperatorSec = swept.lanczosOperatorSec;
        report.poleDecompositionSec = swept.poleDecompositionSec;
        report.orthogonalizationSec = swept.orthogonalizationSec;
        report.romProjectionSec = swept.romProjectionSec;
        report.onlineSweepSec = swept.onlineSweepSec;

        if (auto* reusable = dynamic_cast<linalg::IFactorizedSparseSolver*>(sparseSolver.get())) {
            const auto stats = reusable->factorizationStatistics();
            report.symbolicAnalysisCount = stats.symbolicAnalysisCount;
            report.numericFactorizationCount = stats.numericFactorizationCount;
            report.factorizedRhsSolveCount = stats.factorizedRhsSolveCount;
            report.factorizedSolveCallCount = stats.factorizedSolveCallCount;
            report.batchRhsMax = stats.batchRhsMax;
            report.estimatedSolverPeakKiB = stats.estimatedSolverPeakKiB;
            report.symbolicAnalysisSec = stats.symbolicAnalysisSec;
            report.numericFactorizationSec = stats.numericFactorizationSec;
            report.factorizedRhsSolveSec = stats.factorizedRhsSolveSec;

            std::ostringstream solverSummary;
            solverSummary << std::fixed << std::setprecision(6)
                          << sparseSolver->name() << " reuse: symbolic=" << stats.symbolicAnalysisCount
                          << " (" << stats.symbolicAnalysisSec << " s), numeric="
                          << stats.numericFactorizationCount << " ("
                          << stats.numericFactorizationSec << " s), RHS="
                          << stats.factorizedRhsSolveCount << " vectors / "
                          << stats.factorizedSolveCallCount << " calls, batch_max="
                          << stats.batchRhsMax << " ("
                          << stats.factorizedRhsSolveSec << " s)";
            log.info(solverSummary.str());
        }

        const std::vector<SParameterPoint>& sparams = swept.points;
        const std::vector<std::complex<double>>& lastEdgeDofs = swept.lastEdgeDofs;
        const double lastFrequency = swept.lastFrequencyHz;

        // -------------------- S parameter recap table --------------------
        // Single-glance table for users skimming the log. Tagged "[s-params]"
        // so it is grep-able in scripts. The reported power covers only the displayed modes.
        log.info("S parameters (" + std::to_string(sparams.size())
                 + " rows; power sum covers the displayed modes):");
        log.stat("f(GHz)    |S11|dB     |S21|dB     |S11|^2+|S21|^2");
        for (std::size_t i = 0; i < sparams.size(); ++i) {
            const auto& p = sparams[i];
            const double s11mag = std::abs(p.s11);
            const double s21mag = std::abs(p.s21);
            const double s11Db = 20.0 * std::log10(std::max(s11mag, 1.0e-300));
            const double s21Db = 20.0 * std::log10(std::max(s21mag, 1.0e-300));
            const double passivity = s11mag * s11mag + s21mag * s21mag;
            std::ostringstream row;
            row << std::fixed << std::setprecision(3)
                << std::setw(8)  << (p.frequencyHz / 1.0e9)
                << "  "
                << std::setw(10) << s11Db
                << "  "
                << std::setw(10) << s21Db
                << "  "
                << std::setw(7)  << std::setprecision(4) << passivity;
            log.stat(row.str());
            // Mirror into the JSON FrequencyStat array we're already
            // building per-point; if the strategy didn't fill in a stat for
            // this point (e.g. AlpsSweep skips onFrequencyStat), add one
            // with timing zeroed.
            if (i < sweepStats.size()) {
                sweepStats[i].s11Db = s11Db;
                sweepStats[i].s21Db = s21Db;
            } else {
                FrequencyStat fs;
                fs.frequencyHz = p.frequencyHz;
                fs.s11Db = s11Db;
                fs.s21Db = s21Db;
                sweepStats.push_back(fs);
            }
        }
        report.sweep = sweepStats;

        // -------------------- Result export --------------------
        log.phase("Result export");
        const std::filesystem::path csvPath = options.outDir / "s_parameters.csv";
        OutputWriter::writeSParameters(csvPath, sparams);
        if (!swept.modalScattering.empty()) {
            OutputWriter::writeModalScattering(options.outDir / "s_parameters_modal.csv", swept.modalScattering);
        }
        log.info("Wrote " + csvPath.string());
        if (!lastEdgeDofs.empty()) {
            const std::filesystem::path vtuPath = options.outDir / "field_last.vtu";
            OutputWriter::writeVTU(vtuPath, mesh, project, topology, lastEdgeDofs, lastFrequency, effectiveFieldOrder);
            log.info("Wrote " + vtuPath.string());
        }

        // -------------------- Completed --------------------
        log.phase("Completed");
        log.info("Solver finished, writing run summary and JSON sidecar.");
        log.summary();

        report.phases = collectPhaseSummaries(log);
        const auto memFinal = queryProcessMemory();
        if (memFinal.peakBytes > 0) {
            report.peakMemoryMb = bytesToMb(memFinal.peakBytes);
        }
        report.totalElapsedSec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - runStart).count();
        if (!writeTimingReport(timingJsonPath, report)) {
            log.warn("Failed to write timing report to " + timingJsonPath.string());
        } else {
            log.info("Wrote " + timingJsonPath.string());
        }
        if (!writeRunReport(runJsonPath, report)) {
            log.warn("Failed to write JSON sidecar to " + runJsonPath.string());
        } else {
            log.info("Wrote " + runJsonPath.string());
        }

        return 0;
    } catch (const std::exception& ex) {
        // Fatal path: emit a [fatal] line to console + tee log, then dump
        // the last few hundred lines + active phase + active context to a
        // separate run.crash.log; flag the JSON sidecar with status "fatal"
        // and write whatever sweep stats we managed to collect.
        log.error(std::string("[fatal] ") + ex.what());
        std::cerr << "[fatal] " << ex.what() << '\n';
        if (!crashLogPath.empty()) {
            if (log.writeCrashDump(crashLogPath.string(), ex.what())) {
                log.info("Crash details written to " + crashLogPath.string());
            }
        }
        if (!runJsonPath.empty()) {
            report.status = "fatal";
            report.sweep = sweepStats;
            report.phases = collectPhaseSummaries(log);
            const auto memFinal = queryProcessMemory();
            if (memFinal.peakBytes > 0) {
                report.peakMemoryMb = bytesToMb(memFinal.peakBytes);
            }
            report.totalElapsedSec = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - runStart).count();
            if (!timingJsonPath.empty()) {
                writeTimingReport(timingJsonPath, report);
            }
            writeRunReport(runJsonPath, report);
        }
        return 1;
    }
}

}  // namespace fem
