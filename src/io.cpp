
/*
 * Copyright (c) 2024 Gus Zhang.
 *
 * Summary:
 *   This file provides input/output utility functions for the IPTVisual project.
 *   It includes functions to read wire definitions from files, write trajectory points,
 *   serialize and deserialize L-matrix data, and handle observation grid definitions and values.
 *   The functions facilitate the conversion between file representations and in-memory
 *   data structures used for electromagnetic simulations and visualizations.
 *
 * Functions:
 *   - std::vector<IPTVisual::wire> read_wires(std::string filename):
 *       Reads wire definitions from a file and constructs a vector of wire objects.
 *
 *   - void write_trajectory_points(std::string filename, std::vector<vector_3d> trajectory):
 *       Writes a sequence of 3D trajectory points to a file.
 *
 *   - void write_l_matrix(std::string filename, int n, double *l_matrix, double *len_array):
 *       Serializes an L-matrix and associated length array to a file.
 *
 *   - IPTVisual::observations read_observation_definition(std::string filename):
 *       Reads observation grid definitions from a file and constructs an observations object.
 *
 *   - void write_observation_values(std::string filename, IPTVisual::observations obs):
 *       Writes observation points and their associated values to a file.
 */

#include "io.h"

std::vector<IPTVisual::wire> read_wires(std::string filename)
{
    std::vector<IPTVisual::wire> winding_array;
    std::ifstream file_in(filename);
    double section_radius;
    int number_of_windings;
    int number_of_sectors;
    file_in >> number_of_windings;

    for (int i = 0; i < number_of_windings; i++)
    {
        IPTVisual::wire *wire_obj = new IPTVisual::wire(section_radius);
        file_in >> section_radius >> number_of_sectors;
        wire_obj->section_radius = section_radius;
        for (int sector = 0; sector < number_of_sectors; sector++)
        {
            char mode;
            int number_of_points;
            std::vector<vector_3d> key_points;
            vector_3d vec;
            std::vector<double> length_array;
            double angle;
            double radius;
            file_in >> mode;
            switch (mode)
            {
            case 'l':
                file_in >> number_of_points;
                for (int j = 0; j < number_of_points; j++)
                {
                    file_in >> vec.x >> vec.y >> vec.z;
                    key_points.push_back(vec);
                }
                wire_obj->append_any_lines(key_points);
                break;
            case 's':
                file_in >> number_of_points;
                for (int j = 0; j < number_of_points; j++)
                {
                    file_in >> vec.x >> vec.y >> vec.z;
                    key_points.push_back(vec);
                }
                wire_obj->append_any_spiral(key_points);
                break;
            case 'g': // g for "saggy", with parameter "a" known
                double a;
                file_in >> number_of_points >> a;
                for (int j = 0; j < number_of_points; j++)
                {
                    file_in >> vec.x >> vec.y >> vec.z;
                    key_points.push_back(vec);
                }
                wire_obj->append_any_saggy_lines(key_points, a);
                break;
            case 'c': // c for "catenary", with length of each segment given as the fourth parameter in each line connecting to the previous point.
                double length;
                file_in >> number_of_points;
                for (int j = 0; j < number_of_points; j++)
                {
                    file_in >> vec.x >> vec.y >> vec.z >> length;
                    key_points.push_back(vec);
                    length_array.push_back(length);
                }
                wire_obj->append_any_saggy_lines_with_length(key_points, length_array);
                break;
            case 'x':
                file_in >> angle;
                wire_obj->rotate_x(angle / 180 * pi);
                break;
            case 'y':
                file_in >> angle;
                wire_obj->rotate_y(angle / 180 * pi);
                break;
            case 'z':
                file_in >> angle;
                wire_obj->rotate_z(angle / 180 * pi);
                break;
            case 'o':
                file_in >> vec.x >> vec.y >> vec.z;
                wire_obj->move_offset(vec);
                break;
            case 'n':
                file_in >> vec.x >> vec.y >> vec.z;
                wire_obj->align_normal(vec);
                break;
            case 'b':
                file_in >> radius;
                wire_obj->bend_circle_x(radius);
                break;
            case 'i':
                double real, imag;
                file_in >> real >> imag;
                wire_obj->current = std::complex<double>(real, imag);
            default:
                break;
            }
        }
        wire_obj->update_from_trajectory();
        winding_array.push_back(*wire_obj);
    }
    file_in.close();
    return winding_array;
}

void write_trajectory_points(std::string filename, std::vector<vector_3d> trajectory)
{
    std::ofstream file_out(filename);
    for (auto it = trajectory.begin(); it < trajectory.end(); it++)
    {
        file_out << it->x << '\t' << it->y << '\t' << it->z << std::endl;
    }
    file_out.close();
    return;
}

void write_l_matrix(std::string filename, int n, double *l_matrix, double *len_array)
{
    std::ofstream file_out(filename);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            file_out << l_matrix[j * n + i] << "\t";
        }
        file_out << "\r\n";
    }
    for (int i = 0; i < n; i++)
    {
        file_out << len_array[i] << "\t";
    }
    file_out.close();
    return;
}

IPTVisual::observations read_observation_definition(std::string filename)
{
    std::ifstream file_in(filename);
    double x_min, x_max, y_min, y_max, z_min, z_max;
    int x_num, y_num, z_num;
    file_in >> x_min >> x_max >> x_num >> y_min >> y_max >> y_num >> z_min >> z_max >> z_num;
    IPTVisual::observations obs(x_min, x_max, x_num, y_min, y_max, y_num, z_min, z_max, z_num);
    file_in.close();
    return obs;
}

void write_observation_values(std::string filename, IPTVisual::observations obs)
{
    std::ofstream file_out(filename);
    for (int i = 0; i < obs.observation_points.size(); i++)
    {
        file_out << obs.observation_points[i].x << '\t' << obs.observation_points[i].y << '\t' << obs.observation_points[i].z << '\t' << obs.observation_values[i].x.real() << '\t' << obs.observation_values[i].x.imag() << '\t' << obs.observation_values[i].y.real() << '\t' << obs.observation_values[i].y.imag() << '\t' << obs.observation_values[i].z.real() << '\t' << obs.observation_values[i].z.imag() << std::endl;
    }
    file_out.close();
    return;
}