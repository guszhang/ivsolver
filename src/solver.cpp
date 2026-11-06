/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * High-level solver driver and result writers.
 */
#include "solver.h"
#include "io.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace
{
    using clock_type = std::chrono::steady_clock;

    double ms_since(clock_type::time_point start)
    {
        return std::chrono::duration<double, std::milli>(clock_type::now() - start).count();
    }

    // Compact JSON number formatting; non-finite values become null.
    class json_numbers
    {
    public:
        explicit json_numbers(std::ostream &out) : out(out) {}

        void num(double v, int digits = 17)
        {
            if (!std::isfinite(v))
            {
                out << "null";
                return;
            }
            char buf[40];
            std::snprintf(buf, sizeof(buf), "%.*g", digits, v);
            out << buf;
        }

        void cplx(std::complex<double> c, int digits = 17)
        {
            out << '[';
            num(c.real(), digits);
            out << ',';
            num(c.imag(), digits);
            out << ']';
        }

        void cplx_array(const std::vector<std::complex<double>> &values)
        {
            out << '[';
            for (size_t i = 0; i < values.size(); i++)
            {
                if (i)
                    out << ',';
                cplx(values[i]);
            }
            out << ']';
        }

        // Flattened complex 3-vectors: [xr, xi, yr, yi, zr, zi, ...]
        template <typename Getter>
        void vector_field(size_t count, Getter get, int digits)
        {
            out << '[';
            for (size_t i = 0; i < count; i++)
            {
                const vector_3d_complex &v = get(i);
                if (i)
                    out << ',';
                num(v.x.real(), digits);
                out << ',';
                num(v.x.imag(), digits);
                out << ',';
                num(v.y.real(), digits);
                out << ',';
                num(v.y.imag(), digits);
                out << ',';
                num(v.z.real(), digits);
                out << ',';
                num(v.z.imag(), digits);
            }
            out << ']';
        }

        void points(const std::vector<vector_3d> &pts, int digits)
        {
            out << '[';
            for (size_t i = 0; i < pts.size(); i++)
            {
                if (i)
                    out << ',';
                num(pts[i].x, digits);
                out << ',';
                num(pts[i].y, digits);
                out << ',';
                num(pts[i].z, digits);
            }
            out << ']';
        }

    private:
        std::ostream &out;
    };

    std::string json_escape(const std::string &s)
    {
        std::string r;
        for (char c : s)
        {
            switch (c)
            {
            case '"':
                r += "\\\"";
                break;
            case '\\':
                r += "\\\\";
                break;
            case '\n':
                r += "\\n";
                break;
            case '\r':
                r += "\\r";
                break;
            case '\t':
                r += "\\t";
                break;
            default:
                if (static_cast<unsigned char>(c) < 0x20)
                {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    r += buf;
                }
                else
                {
                    r += c;
                }
            }
        }
        return r;
    }

    const char *mode_name(ivsolver::solve_mode m)
    {
        switch (m)
        {
        case ivsolver::solve_mode::inductance:
            return "inductance";
        case ivsolver::solve_mode::b_field:
            return "b_field";
        default:
            return "poynting";
        }
    }

    constexpr int field_digits = 10;
    constexpr int geometry_digits = 12;

    void write_circuit_file(const std::string &filename, const ivsolver::circuit_result &c)
    {
        std::ofstream out(filename);
        out << "# frequency_hz " << c.frequency << "\n";
        out << "# wire\tdriven_by\tI_re\tI_im\tV_re\tV_im\tZ_re\tZ_im\tC_series\tP_source\tQ_source\tP_resistive\n";
        out.precision(12);
        for (size_t i = 0; i < c.current.size(); i++)
        {
            out << i << '\t' << (c.current_driven[i] ? "current" : "voltage") << '\t'
                << c.current[i].real() << '\t' << c.current[i].imag() << '\t'
                << c.voltage[i].real() << '\t' << c.voltage[i].imag() << '\t'
                << c.self_impedance[i].real() << '\t' << c.self_impedance[i].imag() << '\t'
                << c.capacitance[i] << '\t'
                << c.source_power[i].real() << '\t' << c.source_power[i].imag() << '\t'
                << c.resistive_loss[i] << "\n";
        }
    }
}

namespace ivsolver
{
    bool has_frequency(const options &opt)
    {
        return std::isfinite(opt.frequency) && opt.frequency >= 0;
    }

    result run(std::istream &wire_def, std::istream *obs_def, const options &opt)
    {
        result res;
        auto t0 = clock_type::now();
        res.wires = read_wires(wire_def);
        if (opt.mode != solve_mode::inductance)
        {
            if (!obs_def)
            {
                throw std::runtime_error("an observation definition is required for field calculations");
            }
            IPTVisual::observations obs = read_observation_definition(*obs_def);
            res.points = std::move(obs.observation_points);
        }
        res.geometry_ms = ms_since(t0);
        const size_t n = res.wires.size();

        bool need_circuit = false;
        if (opt.mode == solve_mode::inductance)
        {
            need_circuit = has_frequency(opt);
        }
        else if (opt.mode == solve_mode::poynting)
        {
            if (!has_frequency(opt))
            {
                throw std::runtime_error("the Poynting vector calculation needs an analysis frequency (--freq)");
            }
            need_circuit = opt.force_circuit;
            for (const auto &w : res.wires)
            {
                if (!w.current_defined)
                {
                    need_circuit = true;
                }
            }
        }

        if (opt.mode == solve_mode::inductance || need_circuit)
        {
            auto t1 = clock_type::now();
            res.inductance = compute_inductance_matrix(res.wires, opt.threads);
            res.has_inductance = true;
            res.inductance_ms = ms_since(t1);
        }
        if (need_circuit)
        {
            res.circuit = solve_circuit(res.wires, res.inductance, opt.frequency);
            res.has_circuit = true;
        }

        res.currents.resize(n);
        for (size_t i = 0; i < n; i++)
        {
            res.currents[i] = res.has_circuit ? res.circuit.current[i] : res.wires[i].current;
        }

        auto t2 = clock_type::now();
        if (opt.mode == solve_mode::b_field)
        {
            res.kernel = opt.kernel_given ? opt.kernel : kernel_type::midpoint;
            res.b_values = compute_b_field(res.wires, res.points, res.kernel, opt.threads);
        }
        else if (opt.mode == solve_mode::poynting)
        {
            res.kernel = opt.kernel_given ? opt.kernel : kernel_type::exact;
            res.field = compute_em_field(res.wires, res.currents, res.points, opt.frequency, res.kernel, opt.threads);
        }
        res.field_ms = ms_since(t2);
        return res;
    }

    void write_json(std::ostream &out, const result &res, const options &opt)
    {
        json_numbers j(out);
        const size_t n = res.wires.size();
        out << "{\"solver\":\"ivsolver\",\"version\":\"" IVSOLVER_VERSION "\",\"mode\":\"" << mode_name(opt.mode) << "\",";
        out << "\"frequency\":";
        if (has_frequency(opt))
            j.num(opt.frequency);
        else
            out << "null";
        out << ",\"wires\":[";
        for (size_t i = 0; i < n; i++)
        {
            const auto &w = res.wires[i];
            if (i)
                out << ',';
            out << "{\"radius\":";
            j.num(w.section_radius);
            out << ",\"length\":";
            j.num(w.total_length);
            out << ",\"segments\":" << w.segments.size();
            out << ",\"current_defined\":" << (w.current_defined ? "true" : "false");
            if (opt.include_trajectories)
            {
                out << ",\"trajectory\":";
                j.points(w.trajectory, geometry_digits);
            }
            out << '}';
        }
        out << "],\"inductance\":";
        if (res.has_inductance)
        {
            out << '[';
            for (size_t i = 0; i < n; i++)
            {
                if (i)
                    out << ',';
                out << '[';
                for (size_t k = 0; k < n; k++)
                {
                    if (k)
                        out << ',';
                    j.num(res.inductance[i * n + k]);
                }
                out << ']';
            }
            out << ']';
        }
        else
        {
            out << "null";
        }
        out << ",\"circuit\":";
        if (res.has_circuit)
        {
            const auto &c = res.circuit;
            out << "{\"current\":";
            j.cplx_array(c.current);
            out << ",\"voltage\":";
            j.cplx_array(c.voltage);
            out << ",\"impedance\":";
            j.cplx_array(c.self_impedance);
            out << ",\"capacitance\":[";
            for (size_t i = 0; i < n; i++)
            {
                if (i)
                    out << ',';
                j.num(c.capacitance[i]);
            }
            out << "],\"source_power\":";
            j.cplx_array(c.source_power);
            out << ",\"resistive_loss\":[";
            for (size_t i = 0; i < n; i++)
            {
                if (i)
                    out << ',';
                j.num(c.resistive_loss[i]);
            }
            out << "],\"current_driven\":[";
            for (size_t i = 0; i < n; i++)
            {
                if (i)
                    out << ',';
                out << (c.current_driven[i] ? "true" : "false");
            }
            out << "]}";
        }
        else
        {
            out << "null";
        }
        out << ",\"currents\":";
        j.cplx_array(res.currents);
        out << ",\"field\":";
        if (opt.mode == solve_mode::inductance)
        {
            out << "null";
        }
        else
        {
            const size_t m = res.points.size();
            out << "{\"kernel\":\"" << (res.kernel == kernel_type::exact ? "exact" : "midpoint") << "\",\"count\":" << m << ",\"points\":";
            j.points(res.points, geometry_digits);
            if (opt.mode == solve_mode::b_field)
            {
                out << ",\"B\":";
                j.vector_field(m, [&](size_t i) -> const vector_3d_complex &
                               { return res.b_values[i]; }, field_digits);
            }
            else
            {
                out << ",\"E\":";
                j.vector_field(m, [&](size_t i) -> const vector_3d_complex &
                               { return res.field[i].E; }, field_digits);
                out << ",\"B\":";
                j.vector_field(m, [&](size_t i) -> const vector_3d_complex &
                               { return res.field[i].B; }, field_digits);
                out << ",\"S\":";
                j.vector_field(m, [&](size_t i) -> const vector_3d_complex &
                               { return res.field[i].S; }, field_digits);
                out << ",\"min_distance\":[";
                for (size_t i = 0; i < m; i++)
                {
                    if (i)
                        out << ',';
                    j.num(res.field[i].min_distance, 8);
                }
                out << ']';
            }
            out << '}';
        }
        out << ",\"timing_ms\":{\"geometry\":";
        j.num(res.geometry_ms, 6);
        out << ",\"inductance\":";
        j.num(res.inductance_ms, 6);
        out << ",\"field\":";
        j.num(res.field_ms, 6);
        out << "}}";
    }

    void write_json_error(std::ostream &out, const std::string &message)
    {
        out << "{\"solver\":\"ivsolver\",\"version\":\"" IVSOLVER_VERSION "\",\"error\":\"" << json_escape(message) << "\"}";
    }

    void write_text_files(const std::string &prefix, const result &res, const options &opt)
    {
        const size_t n = res.wires.size();
        if (res.has_inductance)
        {
            std::vector<double> lengths(n);
            for (size_t i = 0; i < n; i++)
            {
                lengths[i] = res.wires[i].total_length;
            }
            std::vector<double> matrix = res.inductance;
            write_l_matrix(prefix + ".lmat", static_cast<int>(n), matrix.data(), lengths.data());
        }
        if (res.has_circuit)
        {
            write_circuit_file(prefix + ".circuit", res.circuit);
        }
        if (opt.mode == solve_mode::b_field)
        {
            IPTVisual::observations obs;
            obs.observation_points = res.points;
            obs.observation_values = res.b_values;
            write_observation_values(prefix + ".obs", obs);
        }
        else if (opt.mode == solve_mode::poynting)
        {
            std::ofstream out(prefix + ".poynting");
            out << "# ivsolver " IVSOLVER_VERSION " complex Poynting vector, frequency_hz " << opt.frequency
                << ", kernel " << (res.kernel == kernel_type::exact ? "exact" : "midpoint") << "\n";
            out << "# x y z Ex_re Ex_im Ey_re Ey_im Ez_re Ez_im Bx_re Bx_im By_re By_im Bz_re Bz_im Sx_re Sx_im Sy_re Sy_im Sz_re Sz_im min_distance\n";
            out.precision(10);
            for (size_t i = 0; i < res.points.size(); i++)
            {
                const auto &p = res.points[i];
                const auto &f = res.field[i];
                out << p.x << '\t' << p.y << '\t' << p.z;
                for (const vector_3d_complex *v : {&f.E, &f.B, &f.S})
                {
                    out << '\t' << v->x.real() << '\t' << v->x.imag()
                        << '\t' << v->y.real() << '\t' << v->y.imag()
                        << '\t' << v->z.real() << '\t' << v->z.imag();
                }
                out << '\t' << f.min_distance << "\n";
            }
        }
        if (opt.include_trajectories)
        {
            for (size_t i = 0; i < n; i++)
            {
                write_trajectory_points(prefix + "_" + std::to_string(i) + ".trj", res.wires[i].trajectory);
            }
        }
    }
}
