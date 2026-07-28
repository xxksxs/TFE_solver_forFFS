// Implements factory::makeSweepStrategy. See header for the contract.
//
// Routing: read Options::sweepStrategy, instantiate the matching concrete
// ISweepStrategy and copy any strategy-specific configuration (order,
// expansion frequency) into its options object. AlpsOptions::expansionFrequencyHz
// is forwarded as 0.0 when the user did not pass --alps-expansion; AlpsSweep::run
// then uses the stable band-center single-point model, see sweep/AlpsSweep.cpp.
//
// Adding a new sweep strategy requires:
//   1. Implement ISweepStrategy in include/bpfem/sweep/<X>Sweep.hpp + .cpp
//   2. Add a value to SweepStrategy enum in app/Application.hpp
//   3. Add a case below
//   4. Add the CLI string in Application::parseOptions
// No other files need to change.

#include "bpfem/factory/SweepStrategyFactory.hpp"

#include "bpfem/sweep/AlpsSweep.hpp"
#include "bpfem/sweep/AweSweep.hpp"
#include "bpfem/sweep/DirectSweep.hpp"
#include "bpfem/sweep/GaweSweep.hpp"
#include "bpfem/sweep/MgaweSweep.hpp"
#include "bpfem/sweep/WcaweSweep.hpp"

namespace fem::factory {

std::unique_ptr<sweep::ISweepStrategy> makeSweepStrategy(
    const Options& options,
    const ProjectDefinition& project,
    const FEMAssembler& assembler,
    const PortModeSolver& portModeSolver) {
    switch (options.sweepStrategy) {
        case SweepStrategy::Alps: {
            sweep::AlpsOptions alpsOpts;
            alpsOpts.order = options.alpsOrder;
            alpsOpts.expansionFrequencyHz = options.alpsExpansionFrequencyHz;
            // ALPS 保存装配器引用，实际稀疏后端由 SweepContext 注入。
            return std::make_unique<sweep::AlpsSweep>(project, assembler, portModeSolver,
                                                       alpsOpts);
        }
        case SweepStrategy::Awe: {
            sweep::AweOptions aweOpts;
            aweOpts.order = options.aweOrder;
            aweOpts.expansionFrequencyHz = options.aweExpansionFrequencyHz;
            return std::make_unique<sweep::AweSweep>(project, assembler, portModeSolver,
                                                     aweOpts);
        }
        case SweepStrategy::Gawe: {
            sweep::GaweOptions gaweOpts;
            gaweOpts.order = options.gaweOrder;
            gaweOpts.expansionFrequencyHz = options.gaweExpansionFrequencyHz;
            gaweOpts.dropTolerance = options.gaweDropTolerance;
            return std::make_unique<sweep::GaweSweep>(project, assembler, portModeSolver,
                                                      gaweOpts);
        }
        case SweepStrategy::Mgawe: {
            sweep::MgaweOptions mgaweOpts;
            mgaweOpts.expansionPointCount = options.mgawePointCount;
            mgaweOpts.localOrder = options.mgaweOrder;
            mgaweOpts.dropTolerance = options.mgaweDropTolerance;
            return std::make_unique<sweep::MgaweSweep>(project, assembler, portModeSolver,
                                                       mgaweOpts);
        }
        case SweepStrategy::Wcawe: {
            sweep::WcaweOptions wcaweOpts;
            wcaweOpts.order = options.wcaweOrder;
            wcaweOpts.expansionFrequencyHz = options.wcaweExpansionFrequencyHz;
            wcaweOpts.dropTolerance = options.wcaweDropTolerance;
            return std::make_unique<sweep::WcaweSweep>(project, assembler, portModeSolver,
                                                       wcaweOpts);
        }
        case SweepStrategy::Direct:
        default:
            // DirectSweep is parameterless; per-frequency assemble + solve
            // is driven entirely by SweepContext at run() time.
            return std::make_unique<sweep::DirectSweep>();
    }
}

}  // namespace fem::factory
