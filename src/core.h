/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Linear magnetic cores by the volume integral method.
 *
 * A core is voxelised into rectangular cells. Each cell carries a uniform
 * magnetisation M = chi H with chi = mu_r - 1. The cell magnetisations follow
 * from the linear system
 *
 *   (1/chi) M_v - sum_u G(r_v - r_u) M_u = H0(r_v)
 *
 * where H0 is the field of the winding currents in free space (Biot-Savart, per
 * ampere) and G is the exact H-field tensor of a uniformly magnetised prism
 * (point dipole beyond a few cell sizes). The matrix is symmetric positive
 * definite for mu_r >= 1, so a Jacobi-preconditioned conjugate gradient solve is
 * used. One solve is made per winding at 1 A; everything else is superposition:
 *
 *   B(r)  = B_air(r) + mu0 (H_M(r) + M(r) inside a cell)
 *   dL_ij = sum_v V_v M_j(v) . B_i(v)     (reciprocity, added to the Neumann matrix)
 *
 * Quasi-static and linear: no saturation, no eddy currents, real mu_r only.
 * The vector potential of the cells (needed for E) uses the point-dipole form,
 * which is approximate close to the core.
 */
#pragma once

#include <complex>
#include <iosfwd>
#include <vector>

#include "field.h"
#include "wire.h"

namespace ivsolver
{
    struct core_body
    {
        enum class shape
        {
            box,
            sphere,
            cylinder
        } kind = shape::box;
        double mu_r = 1;
        double voxel = 0;
        vector_3d p0, p1; // box: min / max corner. sphere, cylinder: p0 = centre.
        double radius = 0, length = 0;
        int axis = 2; // cylinder axis 0, 1, 2 = x, y, z
    };

    // Parses the text after the "cores" keyword (the keyword itself included).
    std::vector<core_body> read_cores(std::istream &in);

    struct core_model
    {
        std::vector<vector_3d> centre;
        std::vector<vector_3d> half; // half size of each cell
        std::vector<double> volume;
        std::vector<double> chi;
        std::vector<std::vector<vector_3d>> b0; // [wire][cell] free-space B per ampere
        std::vector<std::vector<vector_3d>> m;  // [wire][cell] magnetisation per ampere
        std::vector<double> residual;           // [wire] relative CG residual
        std::vector<int> iterations;            // [wire]
        bool empty() const { return centre.empty(); }
    };

    constexpr size_t max_core_cells = 4000;

    core_model build_core_model(const std::vector<core_body> &bodies,
                                const std::vector<IPTVisual::wire> &wires,
                                kernel_type kernel, int threads);

    // L (row-major n x n) += core contribution.
    void add_core_inductance(const core_model &core, std::vector<double> &inductance);

    // Core flux density and vector potential at points for the given phasor currents.
    void core_fields(const core_model &core, const std::vector<std::complex<double>> &currents,
                     const std::vector<vector_3d> &points, int threads,
                     std::vector<vector_3d_complex> &B, std::vector<vector_3d_complex> *A);

    // Adds the core to field samples produced by compute_em_field.
    void apply_core_to_samples(const core_model &core, const std::vector<std::complex<double>> &currents,
                               const std::vector<vector_3d> &points, double frequency, int threads,
                               std::vector<field_sample> &samples);
}
