#pragma once

#include "bpfem/sweep/ISweepStrategy.hpp"

namespace fem::sweep {

// Per-frequency direct sweep. For each frequency:
//   1. assembler.assemble(f) -> A, rhs
//   2. solver.solve(A, rhs)  -> field
//   3. extractor.extract(f, field) -> SParameterPoint
//   4. ctx.onFieldSolved(f, field) if set (Application uses it to write VTU)
// The last-evaluated solution is also retained on SweepResult for Application
// to write `field_last.vtu` after the sweep ends.
//
// This is the path that keeps full-space VTU output capability. ALPS does
// not, by design.
class DirectSweep : public ISweepStrategy {
public:
    DirectSweep() = default;

    SweepResult run(const std::vector<double>& frequencies,
                    const SweepContext& ctx) override;

    const char* name() const override { return "direct"; }
};

}  // namespace fem::sweep

