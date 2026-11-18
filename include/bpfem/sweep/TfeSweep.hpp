#pragma once

#include "bpfem/sweep/ISweepStrategy.hpp"
#include "bpfem/sweep/MgaweSweep.hpp"
#include "bpfem/tfe/TfeSystem.hpp"

namespace fem::sweep {

class TfeSweep final : public ISweepStrategy {
public:
    TfeSweep(const tfe::TfeSystem& system, bool allSources, int sourcePort, int sourceMode);
    TfeSweep(const tfe::TfeSystem& system, bool allSources, int sourcePort, int sourceMode,
             MgaweOptions options, bool proxyDirect, double maxSubbandWidthHz,
             std::vector<std::pair<double, double>> exactBandsHz);
    const char* name() const override {
        return proxyDirect_ ? "TFE-Proxy-Direct" : mgawe_ ? "TFE-MGAWE" : "TFE-Direct";
    }
    SweepResult run(const std::vector<double>& frequencies, const SweepContext& ctx) override;

private:
    const tfe::TfeSystem& system_;
    bool allSources_;
    int sourcePort_;
    int sourceMode_;
    bool mgawe_ = false;
    bool proxyDirect_ = false;
    double maxSubbandWidthHz_ = 0.0;
    std::vector<std::pair<double, double>> exactBandsHz_;
    MgaweOptions options_;
};

}  // namespace fem::sweep
