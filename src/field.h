/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Electromagnetic field evaluation at observation points.
 *
 * For every observation point the solver sums the contributions of all wire
 * segments carrying phasor currents I (peak amplitude, angular frequency w):
 *
 *   B = mu0/(4 pi) * sum I dl x R / |R|^3            (Biot-Savart)
 *   A = mu0/(4 pi) * sum I dl / |R|                  (magnetic vector potential)
 *   E = -j w A                                       (induced electric field)
 *   H = B / mu0
 *   S = 1/2 E x conj(H)                              (complex Poynting vector)
 *
 * Re{S} is the time-averaged power flow density (W/m^2), Im{S} the reactive part.
 * The instantaneous Poynting vector is
 *   s(t) = Re{E e^{jwt}} x Re{H e^{jwt}} = Re{S} + Re{1/2 E x H e^{j2wt}}.
 *
 * E only contains the part induced by the time-varying currents. The
 * electrostatic contribution of charges (e.g. on capacitor plates or due to the
 * distributed winding capacitance) is not modelled, consistent with the
 * quasi-static, lumped-circuit assumption of the inductance calculation.
 *
 * Two integration kernels are provided:
 *   midpoint - each segment is a current element at its midpoint (legacy -b kernel)
 *   exact    - closed-form integration over each straight segment, with the
 *              field of a round conductor of the wire's radius used for points
 *              inside the conductor so values stay finite near the wire.
 */
#pragma once

#include <vector>
#include "wire.h"

namespace ivsolver
{
    enum class kernel_type
    {
        midpoint,
        exact
    };

    struct field_sample
    {
        vector_3d_complex E; // V/m
        vector_3d_complex B; // T
        vector_3d_complex S; // W/m^2, complex Poynting vector 1/2 E x H*
        double min_distance; // distance to the nearest conductor axis (m)
    };

    // Magnetic flux density only (legacy -b behaviour, uses the wire currents).
    std::vector<vector_3d_complex> compute_b_field(const std::vector<IPTVisual::wire> &wires,
                                                   const std::vector<vector_3d> &points,
                                                   kernel_type kernel, int threads);

    // Full E/B/S evaluation. currents[i] is the phasor current of wires[i].
    std::vector<field_sample> compute_em_field(const std::vector<IPTVisual::wire> &wires,
                                               const std::vector<std::complex<double>> &currents,
                                               const std::vector<vector_3d> &points,
                                               double frequency, kernel_type kernel, int threads);
}
