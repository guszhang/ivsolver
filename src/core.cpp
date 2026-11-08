/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * Linear magnetic cores by the volume integral method. See core.h.
 */
#include "core.h"
#include "parallel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace
{
    using cplx = std::complex<double>;
    using ivsolver::core_body;
    using ivsolver::core_model;
    using ivsolver::parallel_chunks;

    const double inv4pi = 1.0 / (4.0 * pi);
    constexpr double far_factor = 8.0; // dipole beyond this many half-sizes

    inline double half_axis(const vector_3d &h, int a) { return a == 0 ? h.x : (a == 1 ? h.y : h.z); }
    inline double comp(const vector_3d &v, int a) { return a == 0 ? v.x : (a == 1 ? v.y : v.z); }

    // ln(s + R) without cancellation for s < 0 (rho2 = R^2 - s^2).
    inline double ln_sr(double s, double R, double rho2)
    {
        return s >= 0 ? std::log(s + R) : std::log(rho2) - std::log(R - s);
    }

    // Field per unit magnetisation of a uniformly magnetised prism with half
    // sizes h and its centre at the origin, evaluated at d. T[3*a+b] = H_a / M_b.
    void prism_tensor(const vector_3d &d, const vector_3d &h, double *T)
    {
        const double scale = std::max({h.x, h.y, h.z});
        const double eps = 1e-12 * scale;
        for (int i = 0; i < 9; i++)
        {
            T[i] = 0;
        }
        for (int b = 0; b < 3; b++)
        {
            const int n = b, t1 = (b + 1) % 3, t2 = (b + 2) % 3;
            const double ht1 = half_axis(h, t1), ht2 = half_axis(h, t2), hn = half_axis(h, n);
            const double dt1 = comp(d, t1), dt2 = comp(d, t2), dn = comp(d, n);
            for (int face = -1; face <= 1; face += 2)
            {
                double w = dn - face * hn;
                if (std::fabs(w) < eps)
                {
                    w = eps;
                }
                const double us[2] = {dt1 + ht1, dt1 - ht1};
                const double vs[2] = {dt2 + ht2, dt2 - ht2};
                double e1 = 0, e2 = 0, en = 0;
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        const double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
                        const double u = us[i], v = vs[j];
                        const double R = std::sqrt(u * u + v * v + w * w);
                        e1 += sign * -ln_sr(v, R, u * u + w * w);
                        e2 += sign * -ln_sr(u, R, v * v + w * w);
                        en += sign * std::atan(u * v / (w * R));
                    }
                }
                // surface charge +M on the +n face, -M on the -n face
                const double s = face * inv4pi;
                T[3 * t1 + b] += s * e1;
                T[3 * t2 + b] += s * e2;
                T[3 * n + b] += s * en;
            }
        }
    }

    void dipole_tensor(const vector_3d &d, double volume, double *T)
    {
        const double r2 = d.x * d.x + d.y * d.y + d.z * d.z;
        const double r = std::sqrt(r2);
        const double n[3] = {d.x / r, d.y / r, d.z / r};
        const double k = volume * inv4pi / (r2 * r);
        for (int a = 0; a < 3; a++)
        {
            for (int b = 0; b < 3; b++)
            {
                T[3 * a + b] = k * (3 * n[a] * n[b] - (a == b ? 1.0 : 0.0));
            }
        }
    }

    // H tensor from a cell (half size h) at offset d = target - source.
    void cell_tensor(const vector_3d &d, const vector_3d &h, double *T)
    {
        const double hmax = std::max({h.x, h.y, h.z});
        const double r2 = d.x * d.x + d.y * d.y + d.z * d.z;
        if (r2 > far_factor * far_factor * hmax * hmax)
        {
            dipole_tensor(d, 8 * h.x * h.y * h.z, T);
        }
        else
        {
            prism_tensor(d, h, T);
        }
    }

    struct tensor_key
    {
        long long v[6];
        bool operator==(const tensor_key &o) const
        {
            for (int i = 0; i < 6; i++)
            {
                if (v[i] != o.v[i])
                    return false;
            }
            return true;
        }
    };
    struct tensor_key_hash
    {
        size_t operator()(const tensor_key &k) const
        {
            uint64_t h = 1469598103934665603ull;
            for (int i = 0; i < 6; i++)
            {
                h ^= static_cast<uint64_t>(k.v[i]) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
            }
            return static_cast<size_t>(h);
        }
    };

    struct interaction
    {
        size_t n = 0;
        std::vector<uint32_t> index; // n*n -> table entry
        std::vector<std::array<double, 9>> table;
    };

    interaction build_interaction(const core_model &c, int threads)
    {
        interaction it;
        const size_t n = c.centre.size();
        it.n = n;
        it.index.resize(n * n);
        double q = 1e300;
        for (const auto &h : c.half)
        {
            q = std::min({q, h.x, h.y, h.z});
        }
        q *= 1e-4;
        std::unordered_map<tensor_key, uint32_t, tensor_key_hash> seen;
        std::vector<std::pair<vector_3d, vector_3d>> reps; // offset, half size
        for (size_t v = 0; v < n; v++)
        {
            for (size_t u = 0; u < n; u++)
            {
                const vector_3d d = c.centre[v] - c.centre[u];
                const vector_3d &h = c.half[u];
                tensor_key k;
                k.v[0] = std::llround(d.x / q);
                k.v[1] = std::llround(d.y / q);
                k.v[2] = std::llround(d.z / q);
                k.v[3] = std::llround(h.x / q);
                k.v[4] = std::llround(h.y / q);
                k.v[5] = std::llround(h.z / q);
                auto found = seen.find(k);
                if (found == seen.end())
                {
                    found = seen.emplace(k, static_cast<uint32_t>(reps.size())).first;
                    reps.push_back({d, h});
                }
                it.index[v * n + u] = found->second;
            }
        }
        it.table.resize(reps.size());
        parallel_chunks(reps.size(), 64, threads, [&](size_t begin, size_t end, size_t)
                        {
            for (size_t k = begin; k < end; k++)
            {
                double T[9];
                cell_tensor(reps[k].first, reps[k].second, T);
                for (int i = 0; i < 9; i++)
                {
                    it.table[k][i] = T[i];
                }
            } });
        return it;
    }

    // out = (1/chi) x - G x
    void apply(const core_model &c, const interaction &it, const std::vector<double> &x,
               std::vector<double> &out, int threads)
    {
        const size_t n = it.n;
        parallel_chunks(n, 16, threads, [&](size_t begin, size_t end, size_t)
                        {
            for (size_t v = begin; v < end; v++)
            {
                double gx = 0, gy = 0, gz = 0;
                const uint32_t *row = &it.index[v * n];
                for (size_t u = 0; u < n; u++)
                {
                    const double *T = it.table[row[u]].data();
                    const double mx = x[3 * u], my = x[3 * u + 1], mz = x[3 * u + 2];
                    gx += T[0] * mx + T[1] * my + T[2] * mz;
                    gy += T[3] * mx + T[4] * my + T[5] * mz;
                    gz += T[6] * mx + T[7] * my + T[8] * mz;
                }
                const double inv = 1.0 / c.chi[v];
                out[3 * v] = x[3 * v] * inv - gx;
                out[3 * v + 1] = x[3 * v + 1] * inv - gy;
                out[3 * v + 2] = x[3 * v + 2] * inv - gz;
            } });
    }

    double dot(const std::vector<double> &a, const std::vector<double> &b)
    {
        double s = 0;
        for (size_t i = 0; i < a.size(); i++)
        {
            s += a[i] * b[i];
        }
        return s;
    }

    // Jacobi-preconditioned conjugate gradients for (1/chi - G) M = H0.
    void solve_cg(const core_model &c, const interaction &it, const std::vector<double> &diag,
                  const std::vector<double> &rhs, std::vector<double> &x, double &residual,
                  int &iterations, int threads)
    {
        const size_t len = rhs.size();
        x.assign(len, 0.0);
        std::vector<double> r = rhs, z(len), p(len), ap(len);
        const double bnorm = std::sqrt(dot(rhs, rhs));
        residual = 0;
        iterations = 0;
        if (bnorm == 0)
        {
            return;
        }
        for (size_t i = 0; i < len; i++)
        {
            z[i] = r[i] / diag[i];
        }
        p = z;
        double rz = dot(r, z);
        for (int k = 0; k < 2000; k++)
        {
            residual = std::sqrt(dot(r, r)) / bnorm;
            if (residual < 1e-9)
            {
                break;
            }
            apply(c, it, p, ap, threads);
            const double alpha = rz / dot(p, ap);
            for (size_t i = 0; i < len; i++)
            {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            for (size_t i = 0; i < len; i++)
            {
                z[i] = r[i] / diag[i];
            }
            const double rz_new = dot(r, z);
            const double beta = rz_new / rz;
            rz = rz_new;
            for (size_t i = 0; i < len; i++)
            {
                p[i] = z[i] + beta * p[i];
            }
            iterations = k + 1;
        }
        residual = std::sqrt(dot(r, r)) / bnorm;
    }

    void voxelise(const core_body &b, core_model &c)
    {
        const double h = b.voxel;
        auto add = [&](const vector_3d &centre, const vector_3d &half)
        {
            c.centre.push_back(centre);
            c.half.push_back(half);
            c.volume.push_back(8 * half.x * half.y * half.z);
            c.chi.push_back(b.mu_r - 1);
            if (c.centre.size() > ivsolver::max_core_cells)
            {
                throw std::runtime_error("too many core cells (limit " + std::to_string(ivsolver::max_core_cells) +
                                         "); increase the voxel size");
            }
        };
        if (b.kind == core_body::shape::box)
        {
            const double ext[3] = {b.p1.x - b.p0.x, b.p1.y - b.p0.y, b.p1.z - b.p0.z};
            int cnt[3];
            for (int a = 0; a < 3; a++)
            {
                if (!(ext[a] > 0))
                {
                    throw std::runtime_error("core box must have positive extent");
                }
                const double need = std::round(ext[a] / h);
                if (need > static_cast<double>(ivsolver::max_core_cells))
                {
                    throw std::runtime_error("too many core cells; increase the voxel size");
                }
                cnt[a] = std::max(1, static_cast<int>(need));
            }
            const vector_3d half(ext[0] / (2 * cnt[0]), ext[1] / (2 * cnt[1]), ext[2] / (2 * cnt[2]));
            for (int i = 0; i < cnt[0]; i++)
                for (int j = 0; j < cnt[1]; j++)
                    for (int k = 0; k < cnt[2]; k++)
                        add(vector_3d(b.p0.x + (2 * i + 1) * half.x, b.p0.y + (2 * j + 1) * half.y,
                                      b.p0.z + (2 * k + 1) * half.z),
                            half);
            return;
        }
        // sphere and cylinder: cubic lattice centred on p0, cells kept when their centre is inside.
        double reach[3];
        for (int a = 0; a < 3; a++)
        {
            reach[a] = b.radius;
        }
        if (b.kind == core_body::shape::cylinder)
        {
            reach[b.axis] = b.length / 2;
        }
        int cnt[3];
        for (int a = 0; a < 3; a++)
        {
            const double need = std::ceil(2 * reach[a] / h);
            if (need > static_cast<double>(ivsolver::max_core_cells))
            {
                throw std::runtime_error("too many core cells; increase the voxel size");
            }
            cnt[a] = std::max(1, static_cast<int>(need));
        }
        const vector_3d half(h / 2, h / 2, h / 2);
        for (int i = 0; i < cnt[0]; i++)
            for (int j = 0; j < cnt[1]; j++)
                for (int k = 0; k < cnt[2]; k++)
                {
                    const double o[3] = {(i - (cnt[0] - 1) / 2.0) * h, (j - (cnt[1] - 1) / 2.0) * h,
                                         (k - (cnt[2] - 1) / 2.0) * h};
                    bool inside;
                    if (b.kind == core_body::shape::sphere)
                    {
                        inside = o[0] * o[0] + o[1] * o[1] + o[2] * o[2] <= b.radius * b.radius;
                    }
                    else
                    {
                        const int a = b.axis, p = (a + 1) % 3, q = (a + 2) % 3;
                        inside = o[p] * o[p] + o[q] * o[q] <= b.radius * b.radius && std::fabs(o[a]) <= b.length / 2;
                    }
                    if (inside)
                    {
                        add(vector_3d(b.p0.x + o[0], b.p0.y + o[1], b.p0.z + o[2]), half);
                    }
                }
    }
}

namespace ivsolver
{
    std::vector<core_body> read_cores(std::istream &in)
    {
        std::string text, line;
        while (std::getline(in, line))
        {
            const size_t hash = line.find('#');
            text += (hash == std::string::npos ? line : line.substr(0, hash)) + "\n";
        }
        std::istringstream s(text);
        std::string keyword;
        std::vector<core_body> bodies;
        if (!(s >> keyword))
        {
            return bodies;
        }
        if (keyword != "cores")
        {
            throw std::runtime_error("expected 'cores' keyword");
        }
        int count;
        if (!(s >> count) || count < 0)
        {
            throw std::runtime_error("invalid core count");
        }
        for (int i = 0; i < count; i++)
        {
            const std::string tag = "core " + std::to_string(i);
            std::string kind;
            core_body b;
            auto need = [&](double &x, const char *what)
            {
                if (!(s >> x))
                    throw std::runtime_error(tag + ": missing " + what);
            };
            if (!(s >> kind))
                throw std::runtime_error(tag + ": missing shape");
            need(b.mu_r, "relative permeability");
            need(b.voxel, "voxel size");
            if (!(b.mu_r >= 1))
                throw std::runtime_error(tag + ": relative permeability must be at least 1");
            if (!(b.voxel > 0))
                throw std::runtime_error(tag + ": voxel size must be positive");
            if (kind == "box")
            {
                b.kind = core_body::shape::box;
                need(b.p0.x, "x0"); need(b.p0.y, "y0"); need(b.p0.z, "z0");
                need(b.p1.x, "x1"); need(b.p1.y, "y1"); need(b.p1.z, "z1");
                if (b.p1.x < b.p0.x) std::swap(b.p0.x, b.p1.x);
                if (b.p1.y < b.p0.y) std::swap(b.p0.y, b.p1.y);
                if (b.p1.z < b.p0.z) std::swap(b.p0.z, b.p1.z);
            }
            else if (kind == "sphere")
            {
                b.kind = core_body::shape::sphere;
                need(b.p0.x, "cx"); need(b.p0.y, "cy"); need(b.p0.z, "cz");
                need(b.radius, "radius");
            }
            else if (kind == "cylinder")
            {
                b.kind = core_body::shape::cylinder;
                need(b.p0.x, "cx"); need(b.p0.y, "cy"); need(b.p0.z, "cz");
                need(b.radius, "radius");
                need(b.length, "length");
                char ax;
                if (!(s >> ax) || (ax != 'x' && ax != 'y' && ax != 'z'))
                    throw std::runtime_error(tag + ": cylinder axis must be x, y or z");
                b.axis = ax - 'x';
            }
            else
            {
                throw std::runtime_error(tag + ": unknown shape '" + kind + "' (box, sphere, cylinder)");
            }
            if ((b.kind != core_body::shape::box) && !(b.radius > 0))
                throw std::runtime_error(tag + ": radius must be positive");
            bodies.push_back(b);
        }
        return bodies;
    }

    core_model build_core_model(const std::vector<core_body> &bodies, const std::vector<IPTVisual::wire> &wires,
                                kernel_type kernel, int threads)
    {
        core_model c;
        for (const auto &b : bodies)
        {
            if (b.mu_r > 1)
            {
                voxelise(b, c);
            }
        }
        if (c.empty())
        {
            return c;
        }
        const size_t n = c.centre.size();
        interaction it = build_interaction(c, threads);
        std::vector<double> diag(3 * n);
        for (size_t v = 0; v < n; v++)
        {
            const double *T = it.table[it.index[v * n + v]].data();
            for (int a = 0; a < 3; a++)
            {
                diag[3 * v + a] = 1.0 / c.chi[v] - T[4 * a];
            }
        }
        c.b0.resize(wires.size());
        c.m.resize(wires.size());
        c.residual.assign(wires.size(), 0.0);
        c.iterations.assign(wires.size(), 0);
        for (size_t w = 0; w < wires.size(); w++)
        {
            c.b0[w].resize(n);
            parallel_chunks(n, 16, threads, [&](size_t begin, size_t end, size_t)
                            {
                for (size_t v = begin; v < end; v++)
                {
                    c.b0[w][v] = unit_b_field(wires[w], c.centre[v], kernel);
                } });
            std::vector<double> rhs(3 * n), x;
            for (size_t v = 0; v < n; v++)
            {
                rhs[3 * v] = c.b0[w][v].x / u0;
                rhs[3 * v + 1] = c.b0[w][v].y / u0;
                rhs[3 * v + 2] = c.b0[w][v].z / u0;
            }
            solve_cg(c, it, diag, rhs, x, c.residual[w], c.iterations[w], threads);
            if (c.residual[w] > 1e-4)
            {
                throw std::runtime_error("core magnetisation solve did not converge (wire " + std::to_string(w) + ")");
            }
            c.m[w].resize(n);
            for (size_t v = 0; v < n; v++)
            {
                c.m[w][v] = vector_3d(x[3 * v], x[3 * v + 1], x[3 * v + 2]);
            }
        }
        return c;
    }

    void add_core_inductance(const core_model &c, std::vector<double> &L)
    {
        const size_t n = c.m.size();
        const size_t cells = c.centre.size();
        std::vector<double> d(n * n, 0.0);
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                double s = 0;
                for (size_t v = 0; v < cells; v++)
                {
                    const vector_3d &m = c.m[j][v], &b = c.b0[i][v];
                    s += c.volume[v] * (m.x * b.x + m.y * b.y + m.z * b.z);
                }
                d[i * n + j] = s;
            }
        }
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < n; j++)
            {
                L[i * n + j] += 0.5 * (d[i * n + j] + d[j * n + i]);
            }
        }
    }

    void core_fields(const core_model &c, const std::vector<cplx> &currents, const std::vector<vector_3d> &points,
                     int threads, std::vector<vector_3d_complex> &B, std::vector<vector_3d_complex> *A)
    {
        const size_t cells = c.centre.size();
        B.assign(points.size(), vector_3d_complex());
        if (A)
        {
            A->assign(points.size(), vector_3d_complex());
        }
        if (cells == 0)
        {
            return;
        }
        std::vector<cplx> mx(cells), my(cells), mz(cells);
        for (size_t w = 0; w < c.m.size(); w++)
        {
            for (size_t v = 0; v < cells; v++)
            {
                mx[v] += currents[w] * c.m[w][v].x;
                my[v] += currents[w] * c.m[w][v].y;
                mz[v] += currents[w] * c.m[w][v].z;
            }
        }
        parallel_chunks(points.size(), 8, threads, [&](size_t begin, size_t end, size_t)
                        {
            for (size_t p = begin; p < end; p++)
            {
                cplx hx = 0, hy = 0, hz = 0, ax = 0, ay = 0, az = 0;
                for (size_t v = 0; v < cells; v++)
                {
                    const vector_3d d = points[p] - c.centre[v];
                    double T[9];
                    cell_tensor(d, c.half[v], T);
                    hx += T[0] * mx[v] + T[1] * my[v] + T[2] * mz[v];
                    hy += T[3] * mx[v] + T[4] * my[v] + T[5] * mz[v];
                    hz += T[6] * mx[v] + T[7] * my[v] + T[8] * mz[v];
                    if (d.x >= -c.half[v].x && d.x < c.half[v].x && d.y >= -c.half[v].y && d.y < c.half[v].y &&
                        d.z >= -c.half[v].z && d.z < c.half[v].z) // half-open: a shared face belongs to one cell
                    {
                        hx += mx[v];
                        hy += my[v];
                        hz += mz[v];
                    }
                    if (A)
                    {
                        // point dipole potential, softened inside one cell size
                        const double hmax = std::max({c.half[v].x, c.half[v].y, c.half[v].z});
                        const double r = std::max(std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z), hmax);
                        const double k = u00 * c.volume[v] / (r * r * r);
                        ax += k * (my[v] * d.z - mz[v] * d.y);
                        ay += k * (mz[v] * d.x - mx[v] * d.z);
                        az += k * (mx[v] * d.y - my[v] * d.x);
                    }
                }
                B[p] = vector_3d_complex(hx * u0, hy * u0, hz * u0);
                if (A)
                {
                    (*A)[p] = vector_3d_complex(ax, ay, az);
                }
            } });
    }

    void apply_core_to_samples(const core_model &c, const std::vector<cplx> &currents,
                               const std::vector<vector_3d> &points, double frequency, int threads,
                               std::vector<field_sample> &samples)
    {
        if (c.empty())
        {
            return;
        }
        std::vector<vector_3d_complex> B, A;
        core_fields(c, currents, points, threads, B, &A);
        const cplx mjw(0, -2 * pi * frequency);
        for (size_t p = 0; p < samples.size(); p++)
        {
            field_sample &s = samples[p];
            s.B.x += B[p].x; s.B.y += B[p].y; s.B.z += B[p].z;
            s.E.x += A[p].x * mjw; s.E.y += A[p].y * mjw; s.E.z += A[p].z * mjw;
            // S = 1/2 E x conj(H), H = B / mu0
            const cplx hx = std::conj(s.B.x) / u0, hy = std::conj(s.B.y) / u0, hz = std::conj(s.B.z) / u0;
            s.S = vector_3d_complex(0.5 * (s.E.y * hz - s.E.z * hy), 0.5 * (s.E.z * hx - s.E.x * hz),
                                    0.5 * (s.E.x * hy - s.E.y * hx));
        }
    }
}
