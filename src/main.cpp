/*
 * Copyright (c) 2024 Gus Zhang.
 *
 * @file main.cpp
 * @brief Main entry point for the IVSolver application.
 *
 * This program processes wire data to compute inductance matrices or magnetic field observations,
 * depending on the selected mode. It reads wire definitions from an input file, performs calculations
 * such as self and mutual inductance or magnetic field at observation points, and writes results to output files.
 * Supported modes:
 *   -l : Compute and output the inductance matrix and wire trajectories.
 *   -b : Compute and output the magnetic field at specified observation points and wire trajectories.
 *
 * Usage:
 *   ivsolver <input_file> <output_file> [-l|-b] [observation_file]
 *
 * Author: Gus Zhang
 */
#include <iostream>
#include "io.h"

int xy_to_index(int n, int x, int y)
{
    return y * n + x;
}

int index_to_x(int n, int index)
{
    return index % n;
}

int index_to_y(int n, int index)
{
    return index / n;
}

int main(int argc, char *argv[])
{
    std::string filename;
    IPTVisual::observations obs_B;
    if (argc != 3 && argc != 5)
    {
        std::cout << "argc error" << std::endl;
    }
    std::string filename_in(argv[1]);
    std::string filename_out(argv[2]);
    std::string mode = argc == 5 ? argv[3] : "-l";
    if ((mode.size()) != 2)
    {
        std::cout << "mode argument error" << std::endl;
        return 1;
    }
    char mode_char = mode[1];
    std::string filename_b_observed = argc == 5 ? argv[4] : "";
    std::vector<IPTVisual::wire> wires = read_wires(filename_in);
    int number_of_wires = wires.size();
    double *l_matrix = new double[number_of_wires * number_of_wires];
    double *len_array = new double[number_of_wires];
    switch (mode_char)
    {
    case 'l':
        for (int i = 0; i < number_of_wires; i++)
        {
            l_matrix[xy_to_index(number_of_wires, i, i)] = wires[i].inductance_with(wires[i]);
            len_array[i] = wires[i].total_length;
        }
        for (int i = 0; i < number_of_wires; i++)
        {
            for (int j = i + 1; j < number_of_wires; j++)
            {
                double l = wires[i].inductance_with(wires[j]);
                l_matrix[xy_to_index(number_of_wires, i, j)] = l;
                l_matrix[xy_to_index(number_of_wires, j, i)] = l;
            }
        }
        for (int i = 0; i < number_of_wires; i++)
        {
            filename = filename_out + "_" + std::to_string(i) + ".trj";
            write_trajectory_points(filename, wires[i].trajectory);
        }
        filename = filename_out + ".lmat";
        write_l_matrix(filename, number_of_wires, l_matrix, len_array);
        break;
    case 'b':
        obs_B = read_observation_definition(filename_b_observed);
        for (int i = 0; i < number_of_wires; i++)
        {
            for (int j = 0; j < obs_B.observation_points.size(); j++)
            {
                obs_B.observation_values[j] = obs_B.observation_values[j] + wires[i].b_field_at_point(obs_B.observation_points[j]);
            }
        }
        filename = filename_out + ".obs";
        write_observation_values(filename, obs_B);
        for (int i = 0; i < number_of_wires; i++)
        {
            filename = filename_out + "_" + std::to_string(i) + ".trj";
            write_trajectory_points(filename, wires[i].trajectory);
        }
        break;
    default:
        std::cout << "mode argument error" << std::endl;
        return 1;
    }
    delete[] l_matrix;
    delete[] len_array;
    return 0;
}