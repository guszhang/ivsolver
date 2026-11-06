/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Inductance matrix assembly (Neumann's formula, parallelised) and circuit solution.
 */
#include "circuit.h"
#include "parallel.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace
{
    using cplx = std::complex<double>;

    // Same kernel as IPTVisual::wire::inductance_with, restricted to rows
    // [begin, end) of wire a's segments.
    double partial_inductance(const IPTVisual::wire &a, const IPTVisual::wire &b, size_t begin, size_t end)
    {
        double L = 0;
        const double radius = a.section_radius;
        for (size_t i = begin; i < end; i++)
        {
            const segment &sa = a.segments[i];
            for (const segment &sb : b.segments)
            {
                if (sa.head.x == sb.head.x && sa.head.y == sb.head.y && sa.head.z == sb.head.z &&
                    sa.rear.x == sb.rear.x && sa.rear.y == sb.rear.y && sa.rear.z == sb.rear.z)
                {
                    const double l = sa.length;
                    L += 2 * u00 * (l * std::log((l + std::sqrt(l * l + radius * radius)) / radius) - std::sqrt(l * l + radius * radius) + l / 4 + radius);
                }
                else
                {
                    const double dx = sa.mid.x - sb.mid.x, dy = sa.mid.y - sb.mid.y, dz = sa.mid.z - sb.mid.z;
                    L += u00 * (sa.vec.x * sb.vec.x + sa.vec.y * sb.vec.y + sa.vec.z * sb.vec.z) / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
            }
        }
        return L;
    }

    // Gaussian elimination with partial pivoting; A is n x n row-major, solved in place.
    std::vector<cplx> solve_linear(std::vector<cplx> A, std::vector<cplx> b, size_t n)
    {
        for (size_t col = 0; col < n; col++)
        {
            size_t pivot = col;
            for (size_t r = col + 1; r < n; r++)
            {
                if (std::abs(A[r * n + col]) > std::abs(A[pivot * n + col]))
                {
                    pivot = r;
                }
            }
            if (std::abs(A[pivot * n + col]) == 0)
            {
                throw std::runtime_error("the circuit equations are singular; check that every voltage-driven loop has a non-zero impedance");
            }
            if (pivot != col)
            {
                for (size_t c = 0; c < n; c++)
                {
                    std::swap(A[col * n + c], A[pivot * n + c]);
                }
                std::swap(b[col], b[pivot]);
            }
            for (size_t r = col + 1; r < n; r++)
            {
                cplx factor = A[r * n + col] / A[col * n + col];
                if (factor == cplx(0))
                {
                    continue;
                }
                for (size_t c = col; c < n; c++)
                {
                    A[r * n + c] -= factor * A[col * n + c];
                }
                b[r] -= factor * b[col];
            }
        }
        std::vector<cplx> x(n);
        for (size_t i = n; i-- > 0;)
        {
            cplx sum = b[i];
            for (size_t c = i + 1; c < n; c++)
            {
                sum -= A[i * n + c] * x[c];
            }
            x[i] = sum / A[i * n + i];
        }
        return x;
    }
}

namespace ivsolver
{
    std::vector<double> compute_inductance_matrix(const std::vector<IPTVisual::wire> &wires, int threads)
    {
        const size_t n = wires.size();
        std::vector<double> matrix(n * n, 0.0);
        // Split every (i, j >= i) pair into chunks of rows of wire i so that a single
        // large winding still uses all threads. Partial sums are stored per chunk and
        // reduced in order, so the result does not depend on the thread count.
        constexpr size_t rows_per_chunk = 64;
        struct job
        {
            size_t i, j, begin, end;
        };
        std::vector<job> jobs;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = i; j < n; j++)
            {
                const size_t rows = wires[i].segments.size();
                for (size_t begin = 0; begin < rows; begin += rows_per_chunk)
                {
                    jobs.push_back({i, j, begin, std::min(rows, begin + rows_per_chunk)});
                }
            }
        }
        std::vector<double> partial(jobs.size(), 0.0);
        parallel_chunks(jobs.size(), 1, threads, [&](size_t begin, size_t end, size_t)
                        {
            for (size_t k = begin; k < end; k++)
            {
                const job &jb = jobs[k];
                partial[k] = partial_inductance(wires[jb.i], wires[jb.j], jb.begin, jb.end);
            } });
        for (size_t k = 0; k < jobs.size(); k++)
        {
            matrix[jobs[k].i * n + jobs[k].j] += partial[k];
        }
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = i + 1; j < n; j++)
            {
                matrix[j * n + i] = matrix[i * n + j];
            }
        }
        return matrix;
    }

    circuit_result solve_circuit(const std::vector<IPTVisual::wire> &wires,
                                 const std::vector<double> &inductance_matrix,
                                 double frequency)
    {
        const size_t n = wires.size();
        if (!(frequency >= 0))
        {
            throw std::runtime_error("the analysis frequency must not be negative");
        }
        const double omega = 2 * pi * frequency;
        circuit_result result;
        result.frequency = frequency;
        result.current.assign(n, 0);
        result.voltage.assign(n, 0);
        result.self_impedance.assign(n, 0);
        result.capacitance.assign(n, 0);
        result.source_power.assign(n, 0);
        result.resistive_loss.assign(n, 0);
        result.current_driven.assign(n, false);

        std::vector<cplx> Z(n * n);
        for (size_t i = 0; i < n; i++)
        {
            const auto &w = wires[i];
            const double L_loop = inductance_matrix[i * n + i] + w.external_inductance;
            double C = w.series_capacitance;
            if (w.tune_frequency >= 0)
            {
                const double f_tune = w.tune_frequency > 0 ? w.tune_frequency : frequency;
                const double w_tune = 2 * pi * f_tune;
                if (!(w_tune > 0) || !(L_loop > 0))
                {
                    throw std::runtime_error("wire " + std::to_string(i) + ": cannot tune a capacitor without a positive frequency and loop inductance");
                }
                C = 1.0 / (w_tune * w_tune * L_loop);
            }
            cplx z(w.series_resistance, omega * L_loop);
            if (C > 0)
            {
                if (omega == 0)
                {
                    throw std::runtime_error("wire " + std::to_string(i) + ": a series capacitor blocks DC; use a positive frequency");
                }
                z += cplx(0, -1.0 / (omega * C));
            }
            result.capacitance[i] = C > 0 ? C : 0;
            result.self_impedance[i] = z;
            for (size_t j = 0; j < n; j++)
            {
                Z[i * n + j] = (i == j) ? z : cplx(0, omega * inductance_matrix[i * n + j]);
            }
            result.current_driven[i] = w.current_defined;
            if (w.current_defined)
            {
                result.current[i] = w.current;
            }
            else
            {
                result.voltage[i] = w.voltage_defined ? w.voltage : cplx(0);
            }
        }

        std::vector<size_t> unknown;
        for (size_t i = 0; i < n; i++)
        {
            if (!result.current_driven[i])
            {
                unknown.push_back(i);
            }
        }
        const size_t m = unknown.size();
        if (m > 0)
        {
            std::vector<cplx> A(m * m), b(m);
            for (size_t r = 0; r < m; r++)
            {
                const size_t i = unknown[r];
                b[r] = result.voltage[i];
                for (size_t j = 0; j < n; j++)
                {
                    if (result.current_driven[j])
                    {
                        b[r] -= Z[i * n + j] * result.current[j];
                    }
                }
                for (size_t c = 0; c < m; c++)
                {
                    A[r * m + c] = Z[i * n + unknown[c]];
                }
            }
            std::vector<cplx> x = solve_linear(A, b, m);
            for (size_t r = 0; r < m; r++)
            {
                result.current[unknown[r]] = x[r];
            }
        }
        for (size_t i = 0; i < n; i++)
        {
            if (result.current_driven[i])
            {
                cplx v = 0;
                for (size_t j = 0; j < n; j++)
                {
                    v += Z[i * n + j] * result.current[j];
                }
                result.voltage[i] = v;
            }
            result.source_power[i] = 0.5 * result.voltage[i] * std::conj(result.current[i]);
            result.resistive_loss[i] = 0.5 * std::norm(result.current[i]) * wires[i].series_resistance;
        }
        return result;
    }
}
