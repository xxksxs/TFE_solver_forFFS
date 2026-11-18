#pragma once

#include "bpfem/fem/FEMAssembler.hpp"

namespace fem::tfe {

struct TfeChannel {
    int portId;
    int mode; // one-based index in the retained spectrum of this port
    bool transverseMagnetic;
    double cutoffHz;
    double relativePermittivity;
    double relativePermeability;
    double cutoffSquared;
};

enum class ModeStatus { Propagating, Evanescent, Cutoff };

class TfeSystem {
public:
    using Vector = std::vector<std::complex<double>>;
    using Trace = std::vector<std::vector<std::pair<int, double>>>;

    TfeSystem(const Mesh& mesh, const ProjectDefinition& project, const EdgeTopology& topology,
              const FEMAssembler& fem, const PortModeSolver& modes, int modeOverride);

    std::size_t interiorCount() const { return interiorCount_; }
    std::size_t traceCount() const { return traceCount_; }
    const std::vector<TfeChannel>& channels() const { return channels_; }
    const std::shared_ptr<const SparsePattern>& pattern() const { return pattern_; }

    ModeStatus status(std::size_t channel, double frequencyHz) const;
    double cutoffGuard(std::size_t channel) const;
    std::complex<double> admittance(std::size_t channel, double frequencyHz) const;
    SparseMatrix assemble(double frequencyHz) const;
    Vector rhs(double frequencyHz, const Vector& incidence) const;
    Vector restoreEdgeField(double frequencyHz, const Vector& state) const;
    Vector outgoing(const Vector& state, const Vector& incidence) const;

private:
    Vector modalScale(double frequencyHz) const;
    std::size_t interiorCount_ = 0;
    std::size_t traceCount_ = 0;
    std::vector<TfeChannel> channels_;
    Trace trace_;
    std::shared_ptr<const SparsePattern> pattern_;
    SparseMatrix stiffness_;
    SparseMatrix mass_;
    SparseMatrix linearFrequency_;
};

}  // namespace fem::tfe
