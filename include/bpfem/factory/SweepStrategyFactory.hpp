#pragma once

#include "bpfem/app/Application.hpp"  // Options, SweepStrategy enum
#include "bpfem/core/Types.hpp"
#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"
#include "bpfem/sweep/ISweepStrategy.hpp"

#include <memory>

namespace fem::factory {

// Builds the sweep strategy that matches Options::sweepStrategy.
//
// DirectSweep is parameterless. AlpsSweep needs the project / assembler /
// port-mode references because its offline phase touches them; we pass them
// in here rather than via SweepContext because they must outlive the run().
std::unique_ptr<sweep::ISweepStrategy> makeSweepStrategy(
    const Options& options,
    const ProjectDefinition& project,
    const FEMAssembler& assembler,
    const PortModeSolver& portModeSolver);

}  // namespace fem::factory

