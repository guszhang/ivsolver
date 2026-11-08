/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * High-level solver interface shared by the command-line programme and the
 * WebAssembly module. Inputs are the plain-text wire and observation definitions
 * documented in README.md; results can be written as the classic text files or as
 * a single JSON document.
 */
#pragma once

#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "circuit.h"
#include "core.h"
#include "field.h"
#include "wire.h"

#define IVSOLVER_VERSION "2.0.0"

namespace ivsolver
{
    enum class solve_mode
    {
        inductance, // -l
        b_field,    // -b (legacy, currents from 'i' directives)
        poynting    // -p (E, B and complex Poynting vector)
    };

    struct options
    {
        solve_mode mode = solve_mode::inductance;
        double frequency = -1;      // Hz; < 0 means "not given"
        int threads = 0;            // 0 = all hardware threads
        kernel_type kernel = kernel_type::exact;
        bool kernel_given = false;  // otherwise -b uses midpoint and -p uses exact
        bool include_trajectories = true;
        bool force_circuit = false; // -p: always compute L and solve the circuit
    };

    struct result
    {
        std::vector<IPTVisual::wire> wires;
        bool has_inductance = false;
        std::vector<double> inductance; // row-major n x n
        bool has_circuit = false;
        circuit_result circuit;
        std::vector<std::complex<double>> currents; // currents used for field evaluation
        std::vector<vector_3d> points;
        std::vector<vector_3d_complex> b_values; // -b
        std::vector<field_sample> field;         // -p
        kernel_type kernel = kernel_type::exact;
        size_t core_cells = 0;       // voxels in magnetic cores (0 = air only)
        double core_residual = 0;    // worst relative CG residual
        double geometry_ms = 0, inductance_ms = 0, field_ms = 0;
    };

    bool has_frequency(const options &opt);

    // obs_def may be null for the inductance mode.
    result run(std::istream &wire_def, std::istream *obs_def, const options &opt);

    void write_json(std::ostream &out, const result &res, const options &opt);
    void write_text_files(const std::string &prefix, const result &res, const options &opt);
    void write_json_error(std::ostream &out, const std::string &message);
}
