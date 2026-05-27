#pragma once

// Deprecated location: use <bpfem/sweep/AlpsSweep.hpp> instead.
//
// `bpfem::mor::AlpsSweep` is preserved as a thin alias to the new
// `bpfem::sweep::AlpsSweep` for one release cycle so that any in-tree or
// external consumers compile unchanged. New code should include
// "bpfem/sweep/AlpsSweep.hpp" and use the fem::sweep namespace.

#include "bpfem/sweep/AlpsSweep.hpp"

namespace fem::mor {

using AlpsOptions [[deprecated("use fem::sweep::AlpsOptions from bpfem/sweep/AlpsSweep.hpp")]] =
    fem::sweep::AlpsOptions;

using AlpsSweep [[deprecated("use fem::sweep::AlpsSweep from bpfem/sweep/AlpsSweep.hpp")]] =
    fem::sweep::AlpsSweep;

}  // namespace fem::mor

