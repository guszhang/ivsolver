/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Inductance matrix assembly and frequency-domain circuit solution.
 *
 * Every wire is one loop with a series impedance
 *   Z_k = R_k + j w (L_kk + L_ext,k) + 1/(j w C_k)          (C_k = 0: no capacitor)
 * and either a current source I_k ('i' directive) or a voltage source V_k
 * ('V' directive, default 0 V, i.e. a passive closed loop). The coupled loops obey
 *   V = (Z_diag + j w M) I
 * Rows of voltage-driven loops are solved for their currents, then the source
 * voltages of current-driven loops are recovered. Phasors are peak amplitudes,
 * so complex power is 1/2 V I*.
 */
#pragma once

#include <complex>
#include <vector>
#include "wire.h"

namespace ivsolver
{
    // Symmetric n x n inductance matrix in row-major order (henries).
    std::vector<double> compute_inductance_matrix(const std::vector<IPTVisual::wire> &wires, int threads);

    struct circuit_result
    {
        double frequency = 0;
        std::vector<std::complex<double>> current;        // A (peak)
        std::vector<std::complex<double>> voltage;        // source voltage, V (peak)
        std::vector<std::complex<double>> self_impedance; // Z_k including the self inductance
        std::vector<double> capacitance;                  // series capacitance used (F), 0 = none
        std::vector<std::complex<double>> source_power;   // 1/2 V I* (W, var)
        std::vector<double> resistive_loss;               // 1/2 |I|^2 R (W)
        std::vector<bool> current_driven;
    };

    circuit_result solve_circuit(const std::vector<IPTVisual::wire> &wires,
                                 const std::vector<double> &inductance_matrix,
                                 double frequency);
}
