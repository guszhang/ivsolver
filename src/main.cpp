/*
 * Copyright (c) 2024-2026 Gus Zhang.
 *
 * @file main.cpp
 * @brief Command-line entry point for IVSolver.
 *
 * Usage:
 *   ivsolver <wire_def> <output> [-l]                      inductance matrix (+ circuit with --freq)
 *   ivsolver <wire_def> <output> -b <observation_def>      magnetic flux density (legacy)
 *   ivsolver <wire_def> <output> -p <observation_def> --freq <Hz>
 *                                                          E, B and complex Poynting vector
 * Options:
 *   -f, --freq <Hz>        analysis frequency; enables the circuit solution
 *   -t, --threads <n>      worker threads (default: all hardware threads)
 *   --kernel <name>        field kernel: exact (default for -p) or midpoint (default for -b)
 *   --circuit              -p: solve the circuit even if every wire has an 'i' current
 *   --json                 write a single JSON document to <output>.json, or to stdout if <output> is '-'
 *   --no-trj               do not write trajectory files / trajectories in JSON
 *   -h, --help             show this help
 *   -v, --version          show the version
 * Use '-' as <wire_def> or <observation_def> to read from standard input.
 *
 * Author: Gus Zhang
 */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "solver.h"

namespace
{
    void print_usage(std::ostream &out)
    {
        out << "IVSolver " IVSOLVER_VERSION " - air-core winding inductance and field solver\n\n"
               "Usage:\n"
               "  ivsolver <wire_def> <output> [-l]                   inductance matrix (+ circuit with --freq)\n"
               "  ivsolver <wire_def> <output> -b <observation_def>   magnetic flux density\n"
               "  ivsolver <wire_def> <output> -p <observation_def> --freq <Hz>\n"
               "                                                      E, B and complex Poynting vector\n\n"
               "Options:\n"
               "  -f, --freq <Hz>     analysis frequency; enables the circuit solution\n"
               "  -t, --threads <n>   worker threads (default: all hardware threads)\n"
               "  --kernel <name>     exact (default for -p) or midpoint (default for -b)\n"
               "  --circuit           -p: solve the circuit even if every wire has an 'i' current\n"
               "  --json              write JSON to <output>.json, or to stdout if <output> is '-'\n"
               "  --no-trj            do not write trajectories\n"
               "  -h, --help          show this help\n"
               "  -v, --version       show the version\n\n"
               "Use '-' as <wire_def> or <observation_def> to read from standard input.\n";
    }

    double parse_double(const std::string &text, const std::string &option)
    {
        try
        {
            size_t used = 0;
            double v = std::stod(text, &used);
            if (used == text.size())
            {
                return v;
            }
        }
        catch (const std::exception &)
        {
        }
        throw std::runtime_error("invalid number '" + text + "' for " + option);
    }

    std::unique_ptr<std::istream> open_input(const std::string &name, const std::string &what)
    {
        if (name == "-")
        {
            auto buffer = std::make_unique<std::stringstream>();
            *buffer << std::cin.rdbuf();
            return buffer;
        }
        auto file = std::make_unique<std::ifstream>(name);
        if (!*file)
        {
            throw std::runtime_error("cannot open " + what + " '" + name + "'");
        }
        return file;
    }
}

int main(int argc, char *argv[])
{
    bool json = false;
    bool json_to_stdout = false;
    try
    {
        std::vector<std::string> positional;
        ivsolver::options opt;
        char mode_char = 0;
        for (int i = 1; i < argc; i++)
        {
            std::string arg = argv[i];
            auto value = [&](const std::string &name) -> std::string
            {
                if (i + 1 >= argc)
                {
                    throw std::runtime_error("missing value for " + name);
                }
                return argv[++i];
            };
            if (arg == "-h" || arg == "--help")
            {
                print_usage(std::cout);
                return 0;
            }
            else if (arg == "-v" || arg == "--version")
            {
                std::cout << "ivsolver " IVSOLVER_VERSION "\n";
                return 0;
            }
            else if (arg == "-l" || arg == "-b" || arg == "-p")
            {
                if (mode_char)
                {
                    throw std::runtime_error("only one of -l, -b and -p may be given");
                }
                mode_char = arg[1];
            }
            else if (arg == "-f" || arg == "--freq")
            {
                opt.frequency = parse_double(value(arg), arg);
                if (!(opt.frequency >= 0))
                {
                    throw std::runtime_error("the frequency must not be negative");
                }
            }
            else if (arg == "-t" || arg == "--threads")
            {
                opt.threads = static_cast<int>(parse_double(value(arg), arg));
            }
            else if (arg == "--kernel")
            {
                std::string k = value(arg);
                if (k == "exact")
                    opt.kernel = ivsolver::kernel_type::exact;
                else if (k == "midpoint")
                    opt.kernel = ivsolver::kernel_type::midpoint;
                else
                    throw std::runtime_error("unknown kernel '" + k + "' (use exact or midpoint)");
                opt.kernel_given = true;
            }
            else if (arg == "--circuit")
            {
                opt.force_circuit = true;
            }
            else if (arg == "--json")
            {
                json = true;
            }
            else if (arg == "--no-trj")
            {
                opt.include_trajectories = false;
            }
            else if (arg.size() > 1 && arg[0] == '-')
            {
                throw std::runtime_error("unknown option '" + arg + "'");
            }
            else
            {
                positional.push_back(arg);
            }
        }

        switch (mode_char)
        {
        case 0:
        case 'l':
            opt.mode = ivsolver::solve_mode::inductance;
            break;
        case 'b':
            opt.mode = ivsolver::solve_mode::b_field;
            break;
        case 'p':
            opt.mode = ivsolver::solve_mode::poynting;
            break;
        }
        const bool needs_obs = opt.mode != ivsolver::solve_mode::inductance;
        // The legacy form "ivsolver in out -l x" carried a dummy fourth argument.
        if (positional.size() < 2 || positional.size() > 3 || (needs_obs && positional.size() != 3))
        {
            print_usage(std::cerr);
            return 1;
        }
        const std::string &output = positional[1];
        json_to_stdout = json && output == "-";
        if (output == "-" && !json)
        {
            throw std::runtime_error("output '-' (stdout) requires --json");
        }
        if (positional[0] == "-" && needs_obs && positional[2] == "-")
        {
            throw std::runtime_error("only one input can be read from standard input");
        }

        auto wire_in = open_input(positional[0], "wire definition");
        std::unique_ptr<std::istream> obs_in;
        if (needs_obs)
        {
            obs_in = open_input(positional[2], "observation definition");
        }
        ivsolver::result res = ivsolver::run(*wire_in, obs_in.get(), opt);

        if (json)
        {
            if (json_to_stdout)
            {
                ivsolver::write_json(std::cout, res, opt);
                std::cout << std::endl;
            }
            else
            {
                std::ofstream out(output + ".json");
                ivsolver::write_json(out, res, opt);
                out << "\n";
            }
        }
        else
        {
            ivsolver::write_text_files(output, res, opt);
        }
        return 0;
    }
    catch (const std::exception &error)
    {
        if (json_to_stdout)
        {
            ivsolver::write_json_error(std::cout, error.what());
            std::cout << std::endl;
        }
        std::cerr << "error: " << error.what() << std::endl;
        return 1;
    }
}
