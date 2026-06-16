#pragma once

#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

#include <vector>

namespace fem::fastsweep {

const PortMode& virtualPortMode(const PortModeSolver& portModeSolver,
                                const FEMAssembler::AffineSystem& affine,
                                int virtualPortIndex);

std::vector<int> dominantVirtualPortsByProject(const PortModeSolver& portModeSolver,
                                               const FEMAssembler::AffineSystem& affine,
                                               int projectPortCount);

}  // namespace fem::fastsweep
