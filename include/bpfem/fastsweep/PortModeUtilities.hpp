#pragma once

#include "bpfem/fem/FEMAssembler.hpp"
#include "bpfem/fem/PortModeSolver.hpp"

#include <vector>

namespace fem::fastsweep {

// 按虚拟端口编号取得对应端口模式；支持 TFE 多模端口。
const PortMode& virtualPortMode(const PortModeSolver& portModeSolver,
                                const FEMAssembler::AffineSystem& affine,
                                int virtualPortIndex);

// 为每个工程端口找出用于 S 参数提取的主导虚拟端口索引。
std::vector<int> dominantVirtualPortsByProject(const PortModeSolver& portModeSolver,
                                               const FEMAssembler::AffineSystem& affine,
                                               int projectPortCount);

}  // namespace fem::fastsweep
