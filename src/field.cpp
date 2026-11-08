/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Field kernels and observation-point evaluation. See field.h for the physics.
 */
#include "field.h"
#include "parallel.h"

#include <cmath>
#include <limits>

namespace
{
    using cplx = std::complex<double>;

    // Field of one wire per ampere of current. Geometry is real, so the per-wire
    // contributions are real vectors that are scaled by the complex current later.
    struct unit_field
    {
        vector_3d A; // vector potential per ampere (T m / A)
        vector_3d B; // flux density per ampere (T / A)
        double min_distance;
    };

    inline vector_3d cross(const vector_3d &a, const vector_3d &b)
    {
        return vector_3d(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }

    inline double dot(const vector_3d &a, const vector_3d &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    inline vector_3d_complex cross(const vector_3d_complex &a, const vector_3d_complex &b)
    {
        return vector_3d_complex(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }

    inline vector_3d_complex conj(const vector_3d_complex &a)
    {
        return vector_3d_complex(std::conj(a.x), std::conj(a.y), std::conj(a.z));
    }

    inline void accumulate(vector_3d_complex &target, const vector_3d &v, cplx scale)
    {
        target.x += v.x * scale;
        target.y += v.y * scale;
        target.z += v.z * scale;
    }

    // Current element at the segment midpoint (legacy kernel).
    unit_field midpoint_field(const IPTVisual::wire &w, const vector_3d &point, bool want_potential)
    {
        unit_field f;
        f.min_distance = std::numeric_limits<double>::infinity();
        for (const segment &seg : w.segments)
        {
            vector_3d R(point.x - seg.mid.x, point.y - seg.mid.y, point.z - seg.mid.z);
            double r = std::sqrt(dot(R, R));
            vector_3d dB = cross(seg.vec, R);
            double k3 = u00 / (r * r * r);
            f.B.x += dB.x * k3;
            f.B.y += dB.y * k3;
            f.B.z += dB.z * k3;
            if (want_potential)
            {
                double k1 = u00 / r;
                f.A.x += seg.vec.x * k1;
                f.A.y += seg.vec.y * k1;
                f.A.z += seg.vec.z * k1;
            }
            f.min_distance = std::min(f.min_distance, r);
        }
        return f;
    }

    // Closed-form integration over each straight segment. For a segment of
    // length L along unit vector u, with the observation point at axial
    // coordinate s (measured from the head) and perpendicular distance d:
    //   B = mu0/(4 pi) (u x r1) / d^2 * [ (L - s)/R2 + s/R1 ]
    //   A = mu0/(4 pi) u * [ asinh((L - s)/d) + asinh(s/d) ]
    // Inside the conductor (d < a) the distance is replaced by the radius a,
    // which gives the field of a uniformly distributed current: B grows
    // linearly with d and A stays finite.
    unit_field exact_field(const IPTVisual::wire &w, const vector_3d &point, bool want_potential)
    {
        unit_field f;
        f.min_distance = std::numeric_limits<double>::infinity();
        const double a2 = w.section_radius * w.section_radius;
        for (const segment &seg : w.segments)
        {
            const double L = seg.length;
            if (!(L > 0))
            {
                continue;
            }
            const vector_3d u(seg.vec.x / L, seg.vec.y / L, seg.vec.z / L);
            const vector_3d r1(point.x - seg.head.x, point.y - seg.head.y, point.z - seg.head.z);
            const double s = dot(r1, u);
            const double ls = L - s;
            const double rr = dot(r1, r1);
            const double d2 = std::max(0.0, rr - s * s);
            const double d2e = std::max(d2, a2);
            const double R1 = std::sqrt(s * s + d2e);
            const double R2 = std::sqrt(ls * ls + d2e);

            const double fb = ls / R2 + s / R1;
            const vector_3d ub = cross(u, r1);
            const double kb = u00 * fb / d2e;
            f.B.x += ub.x * kb;
            f.B.y += ub.y * kb;
            f.B.z += ub.z * kb;

            if (want_potential)
            {
                // asinh(s/d) = ln((s + R1)/d); use the conjugate form when s < 0
                // to avoid cancellation, likewise for the far end.
                const double t1 = s >= 0 ? s + R1 : d2e / (R1 - s);
                const double t2 = ls >= 0 ? ls + R2 : d2e / (R2 - ls);
                const double ka = u00 * std::log(t1 * t2 / d2e);
                f.A.x += u.x * ka;
                f.A.y += u.y * ka;
                f.A.z += u.z * ka;
            }

            const double sc = std::min(std::max(s, 0.0), L);
            const double dx = r1.x - u.x * sc, dy = r1.y - u.y * sc, dz = r1.z - u.z * sc;
            f.min_distance = std::min(f.min_distance, std::sqrt(dx * dx + dy * dy + dz * dz));
        }
        return f;
    }

    inline unit_field wire_field(const IPTVisual::wire &w, const vector_3d &point, ivsolver::kernel_type kernel, bool want_potential)
    {
        return kernel == ivsolver::kernel_type::exact ? exact_field(w, point, want_potential)
                                                      : midpoint_field(w, point, want_potential);
    }

    // Points per scheduling chunk.
    constexpr size_t point_chunk = 16;
}

namespace ivsolver
{
    vector_3d unit_b_field(const IPTVisual::wire &w, const vector_3d &point, kernel_type kernel)
    {
        if (w.segments.empty())
        {
            return vector_3d();
        }
        return wire_field(w, point, kernel, false).B;
    }

    std::vector<vector_3d_complex> compute_b_field(const std::vector<IPTVisual::wire> &wires,
                                                   const std::vector<vector_3d> &points,
                                                   kernel_type kernel, int threads)
    {
        std::vector<vector_3d_complex> values(points.size());
        parallel_chunks(points.size(), point_chunk, threads, [&](size_t begin, size_t end, size_t)
                        {
            for (size_t p = begin; p < end; p++)
            {
                vector_3d_complex B;
                for (const auto &w : wires)
                {
                    unit_field f = wire_field(w, points[p], kernel, false);
                    accumulate(B, f.B, w.current);
                }
                values[p] = B;
            } });
        return values;
    }

    std::vector<field_sample> compute_em_field(const std::vector<IPTVisual::wire> &wires,
                                               const std::vector<std::complex<double>> &currents,
                                               const std::vector<vector_3d> &points,
                                               double frequency, kernel_type kernel, int threads)
    {
        const double omega = 2 * pi * frequency;
        const cplx minus_j_omega(0, -omega);
        std::vector<field_sample> samples(points.size());
        parallel_chunks(points.size(), point_chunk, threads, [&](size_t begin, size_t end, size_t)
                        {
            for (size_t p = begin; p < end; p++)
            {
                vector_3d_complex A, B;
                double min_distance = std::numeric_limits<double>::infinity();
                for (size_t i = 0; i < wires.size(); i++)
                {
                    if (wires[i].segments.empty())
                    {
                        continue;
                    }
                    unit_field f = wire_field(wires[i], points[p], kernel, true);
                    accumulate(A, f.A, currents[i]);
                    accumulate(B, f.B, currents[i]);
                    min_distance = std::min(min_distance, f.min_distance);
                }
                field_sample &out = samples[p];
                out.B = B;
                out.E = vector_3d_complex(A.x * minus_j_omega, A.y * minus_j_omega, A.z * minus_j_omega);
                vector_3d_complex H(B.x / u0, B.y / u0, B.z / u0);
                vector_3d_complex S = cross(out.E, conj(H));
                out.S = vector_3d_complex(S.x * 0.5, S.y * 0.5, S.z * 0.5);
                out.min_distance = min_distance;
            } });
        return samples;
    }
}
