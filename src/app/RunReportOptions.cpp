#include "bpfem/app/Application.hpp"
#include "bpfem/core/RunEnvironment.hpp"
#include "bpfem/core/RunReport.hpp"

#include <charconv>
#include <stdexcept>

namespace fem {

namespace {

std::string formatOption(double value) {
    char buffer[64];
    const auto formatted = std::to_chars(buffer, buffer + sizeof(buffer), value);
    if (formatted.ec != std::errc{}) {
        throw std::runtime_error("Cannot format a floating-point solver option");
    }
    return {buffer, formatted.ptr};
}

}  // namespace

void fillReportFromEnvironment(RunReport& report, const RunEnvironment& env, const Options& options) {
    report.solverVersion = env.solverVersion;
    report.buildString = env.buildString;
    report.runtimeString = env.runtimeString;
    report.hostName = env.hostName;
    report.osDescription = env.osDescription;
    report.cpuDescription = env.cpuDescription;
    report.memoryDescription = env.memoryDescription;
    report.commandLine = env.commandLine;
    report.workingDirectory = env.workingDirectory;

    auto add = [&](const char* k, const std::string& v) {
        report.options.emplace_back(k, v);
    };
    add("aedt", options.aedt.string());
    add("mesh", options.mesh.string());
    add("outDir", options.outDir.string());
    add("maxSweepPoints", std::to_string(options.maxSweepPoints));
    add("maxIterations", std::to_string(options.maxIterations));
    add("basisOrder", std::to_string(options.basisOrder));
    add("fieldOutputOrder", std::to_string(options.fieldOutputOrder));
    add("tolerance", formatOption(options.tolerance));
    add("writeAllFields", options.writeAllFields ? "true" : "false");
    {
        std::string s = "direct";
        switch (options.sweepStrategy) {
            case SweepStrategy::Direct: s = "direct"; break;
            case SweepStrategy::Alps:   s = "alps"; break;
            case SweepStrategy::Awe:    s = "awe"; break;
            case SweepStrategy::Gawe:   s = "gawe"; break;
            case SweepStrategy::Mgawe:  s = "mgawe"; break;
            case SweepStrategy::Wcawe:  s = "wcawe"; break;
        }
        add("sweepStrategy", s);
    }
    add("alpsOrder", std::to_string(options.alpsOrder));
    add("alpsExpansionFrequencyHz", std::to_string(options.alpsExpansionFrequencyHz));
    add("aweOrder", std::to_string(options.aweOrder));
    add("aweExpansionFrequencyHz", std::to_string(options.aweExpansionFrequencyHz));
    add("gaweOrder", std::to_string(options.gaweOrder));
    add("gaweExpansionFrequencyHz", std::to_string(options.gaweExpansionFrequencyHz));
    add("gaweDropTolerance", formatOption(options.gaweDropTolerance));
    add("mgawePointCount", std::to_string(options.mgawePointCount));
    add("mgaweOrder", std::to_string(options.mgaweOrder));
    add("mgaweDropTolerance", formatOption(options.mgaweDropTolerance));
    add("wcaweOrder", std::to_string(options.wcaweOrder));
    add("wcaweExpansionFrequencyHz", std::to_string(options.wcaweExpansionFrequencyHz));
    add("wcaweDropTolerance", formatOption(options.wcaweDropTolerance));
    switch (options.portMethod) {
        case PortMethod::Numerical: report.portFormulation = "numerical"; break;
        case PortMethod::Analytic: report.portFormulation = "analytic"; break;
        case PortMethod::NumericalMultimode: report.portFormulation = "numerical-multimode"; break;
        case PortMethod::TfeBlock: report.portFormulation = "tfe-block"; break;
    }
    add("portMethod", options.portMethod == PortMethod::NumericalMultimode ? "tfe" : report.portFormulation);
    add("tfeModesPerPort", std::to_string(options.tfeModesPerPort));
    add("tfeSubbandWidthHz", std::to_string(options.tfeSubbandWidthHz));
    add("tfeProxyDirect", options.tfeProxyDirect ? "true" : "false");
    {
        std::string ls = "auto";
        switch (options.linearSolver) {
            case LinearSolverKind::Pardiso:  ls = "direct"; break;
            case LinearSolverKind::BiCGStab: ls = "bicgstab"; break;
            case LinearSolverKind::Gmres:    ls = "gmres"; break;
            case LinearSolverKind::Mumps:    ls = "mumps"; break;
            case LinearSolverKind::Auto:     ls = "auto"; break;
        }
        add("linearSolver", ls);
    }
    {
        std::string p = "none";
        switch (options.preconditioner) {
            case PreconditionerKind::Jacobi: p = "jacobi"; break;
            case PreconditionerKind::Ilu0:   p = "ilu0"; break;
            case PreconditionerKind::Schwarz: p = "schwarz"; break;
            case PreconditionerKind::PardisoSingle: p = "pardiso-single"; break;
            case PreconditionerKind::None:   p = "none"; break;
        }
        add("preconditioner", p);
    }
    add("gmresRestart", std::to_string(options.gmresRestart));
    add("iterativeSpmv", options.iterativeMklSpmv ? "mkl" : "native");
    add("pardisoMode", options.pardisoOutOfCore ? "ooc" : "incore");
    add("mumpsMode", options.mumpsMode);
    add("mumpsBlrTolerance", formatOption(options.mumpsBlrTolerance));
    add("schwarzSubdomains", std::to_string(options.schwarzSubdomains));
    add("schwarzOverlap", std::to_string(options.schwarzOverlap));
    add("schwarzLocalSolver", options.schwarzLocalMumps ? "mumps" : "pardiso");
}

}  // namespace fem
