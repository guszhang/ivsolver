/*
 * Copyright (c) 2026 Gus Zhang.
 *
 * WebAssembly entry point (built with Emscripten, see `make wasm`).
 *
 * JavaScript usage:
 *   const ivsolver = await createIVSolver();
 *   const json = ivsolver.solve(wireDefText, obsDefText, "p", 85000, "exact", true);
 *   const result = JSON.parse(json);   // same document as `ivsolver ... --json`
 *
 * mode: "l", "b" or "p"; frequency < 0 means "not given"; kernel "" selects the
 * mode's default. Errors are returned as {"error": "..."} instead of exceptions.
 */
#include <sstream>
#include <string>

#include <emscripten/bind.h>

#include "solver.h"

namespace
{
    std::string solve(const std::string &wire_def, const std::string &obs_def, const std::string &mode,
                      double frequency, const std::string &kernel, bool include_trajectories)
    {
        std::ostringstream out;
        try
        {
            ivsolver::options opt;
            if (mode == "l" || mode.empty())
                opt.mode = ivsolver::solve_mode::inductance;
            else if (mode == "b")
                opt.mode = ivsolver::solve_mode::b_field;
            else if (mode == "p")
                opt.mode = ivsolver::solve_mode::poynting;
            else
                throw std::runtime_error("unknown mode '" + mode + "'");
            opt.frequency = frequency;
            if (kernel == "exact")
            {
                opt.kernel = ivsolver::kernel_type::exact;
                opt.kernel_given = true;
            }
            else if (kernel == "midpoint")
            {
                opt.kernel = ivsolver::kernel_type::midpoint;
                opt.kernel_given = true;
            }
            else if (!kernel.empty())
            {
                throw std::runtime_error("unknown kernel '" + kernel + "'");
            }
            opt.include_trajectories = include_trajectories;
            opt.threads = 1;
            std::istringstream wire_in(wire_def);
            std::istringstream obs_in(obs_def);
            ivsolver::result res = ivsolver::run(wire_in, &obs_in, opt);
            ivsolver::write_json(out, res, opt);
        }
        catch (const std::exception &error)
        {
            out.str("");
            ivsolver::write_json_error(out, error.what());
        }
        return out.str();
    }

    std::string version()
    {
        return IVSOLVER_VERSION;
    }
}

EMSCRIPTEN_BINDINGS(ivsolver_module)
{
    emscripten::function("solve", &solve);
    emscripten::function("version", &version);
}
