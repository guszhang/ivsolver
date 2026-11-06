
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
 *   - void write_observation_values(std::string filename, const IPTVisual::observations &obs):
 *       Writes observation points and their associated values to a file.
 */

#include "io.h"
#include <sstream>

namespace
{
// Reads the next value and raises a descriptive error if the stream runs dry or
// holds something that is not a number.
template <typename T>
void expect(std::istream &in, T &value, const std::string &what)
{
    if (!(in >> value))
    {
        throw std::runtime_error("unexpected end of input or invalid value while reading " + what);
    }
}

void expect_point(std::istream &in, vector_3d &vec, const std::string &what)
{
    expect(in, vec.x, what);
    expect(in, vec.y, what);
    expect(in, vec.z, what);
}

void expect_count(std::istream &in, int &count, const std::string &what)
{
    expect(in, count, what);
    if (count < 0)
    {
        throw std::runtime_error("negative count while reading " + what);
    }
}

// Copies the stream into memory with '#' comments (to the end of the line) removed.
std::istringstream strip_comments(std::istream &in)
{
    std::string text, line;
    while (std::getline(in, line))
    {
        size_t hash = line.find('#');
        text += (hash == std::string::npos ? line : line.substr(0, hash));
        text += '\n';
    }
    return std::istringstream(text);
}

std::vector<vector_3d> read_points(std::istream &in, int count, const std::string &what)
{
    std::vector<vector_3d> key_points;
    vector_3d vec;
    for (int j = 0; j < count; j++)
    {
        expect_point(in, vec, what);
        key_points.push_back(vec);
    }
    return key_points;
}
}

std::vector<IPTVisual::wire> read_wires(std::istream &raw_in)
{
    std::istringstream file_in = strip_comments(raw_in);
    std::vector<IPTVisual::wire> winding_array;
    double section_radius;
    int number_of_windings;
    int number_of_sectors;
    expect(file_in, number_of_windings, "the number of wires");
    if (number_of_windings < 0)
    {
        throw std::runtime_error("the number of wires must not be negative");
    }

    for (int i = 0; i < number_of_windings; i++)
    {
        const std::string wire_tag = "wire " + std::to_string(i);
        expect(file_in, section_radius, wire_tag + " cross-section radius");
        expect_count(file_in, number_of_sectors, wire_tag + " directive count");
        if (!(section_radius > 0))
        {
            throw std::runtime_error(wire_tag + ": cross-section radius must be positive");
        }
        IPTVisual::wire wire_obj(section_radius);
        for (int sector = 0; sector < number_of_sectors; sector++)
        {
            char mode;
            int number_of_points;
            std::vector<vector_3d> key_points;
            vector_3d vec;
            std::vector<double> length_array;
            double angle;
            double radius;
            expect(file_in, mode, wire_tag + " directive " + std::to_string(sector));
            const std::string what = wire_tag + " directive " + std::to_string(sector) + " '" + mode + "'";
            switch (mode)
            {
            case 'l':
                expect_count(file_in, number_of_points, what);
                key_points = read_points(file_in, number_of_points, what);
                wire_obj.append_any_lines(key_points);
                break;
            case 's':
                expect_count(file_in, number_of_points, what);
                key_points = read_points(file_in, number_of_points, what);
                wire_obj.append_any_spiral(key_points);
                break;
            case 'a': // a for "arc", rotating around the z-direction axis at (axis_x, axis_y), in the given direction (+1 or -1).
            {
                double axis_x, axis_y;
                int direction;
                expect(file_in, axis_x, what);
                expect(file_in, axis_y, what);
                expect(file_in, direction, what);
                expect_count(file_in, number_of_points, what);
                key_points = read_points(file_in, number_of_points, what);
                wire_obj.append_any_arc(key_points, axis_x, axis_y, direction);
                break;
            }
            case 'g': // g for "saggy", with parameter "a" known
            {
                double a;
                expect_count(file_in, number_of_points, what);
                expect(file_in, a, what);
                key_points = read_points(file_in, number_of_points, what);
                wire_obj.append_any_saggy_lines(key_points, a);
                break;
            }
            case 'c': // c for "catenary", with length of each segment given as the fourth parameter in each line connecting to the previous point.
            {
                double length;
                expect_count(file_in, number_of_points, what);
                for (int j = 0; j < number_of_points; j++)
                {
                    expect_point(file_in, vec, what);
                    expect(file_in, length, what);
                    key_points.push_back(vec);
                    length_array.push_back(length);
                }
                wire_obj.append_any_saggy_lines_with_length(key_points, length_array);
                break;
            }
            case 'e':
            {
                vector_3d center, normal, x_direction;
                double major_axis_length, minor_axis_length, start_angle, end_angle;
                expect_point(file_in, center, what);
                expect_point(file_in, normal, what);
                expect_point(file_in, x_direction, what);
                expect(file_in, major_axis_length, what);
                expect(file_in, minor_axis_length, what);
                expect(file_in, start_angle, what);
                expect(file_in, end_angle, what);
                wire_obj.append_ellipse(center, normal, x_direction, major_axis_length, minor_axis_length, start_angle, end_angle);
                break;
            }
            case 'p':
            {
                vector_3d origin, normal, x_direction;
                double a, b, c, x1, x2;
                expect_point(file_in, origin, what);
                expect_point(file_in, normal, what);
                expect_point(file_in, x_direction, what);
                expect(file_in, a, what);
                expect(file_in, b, what);
                expect(file_in, c, what);
                expect(file_in, x1, what);
                expect(file_in, x2, what);
                wire_obj.append_parabola(origin, normal, x_direction, a, b, c, x1, x2);
                break;
            }
            case 'x':
                expect(file_in, angle, what);
                wire_obj.rotate_x(angle / 180 * pi);
                break;
            case 'y':
                expect(file_in, angle, what);
                wire_obj.rotate_y(angle / 180 * pi);
                break;
            case 'z':
                expect(file_in, angle, what);
                wire_obj.rotate_z(angle / 180 * pi);
                break;
            case 'o':
                expect_point(file_in, vec, what);
                wire_obj.move_offset(vec);
                break;
            case 'n':
                expect_point(file_in, vec, what);
                wire_obj.align_normal(vec);
                break;
            case 'b':
                expect(file_in, radius, what);
                wire_obj.bend_circle_x(radius);
                break;
            case 'i':
            {
                double real, imag;
                expect(file_in, real, what);
                expect(file_in, imag, what);
                wire_obj.current = std::complex<double>(real, imag);
                wire_obj.current_defined = true;
                break;
            }
            // Circuit directives (upper case). They describe the loop closed by
            // this wire and are only used when an analysis frequency is given.
            case 'V':
            {
                double real, imag;
                expect(file_in, real, what);
                expect(file_in, imag, what);
                wire_obj.voltage = std::complex<double>(real, imag);
                wire_obj.voltage_defined = true;
                break;
            }
            case 'R':
                expect(file_in, wire_obj.series_resistance, what);
                break;
            case 'C':
                expect(file_in, wire_obj.series_capacitance, what);
                break;
            case 'L':
                expect(file_in, wire_obj.external_inductance, what);
                break;
            case 'T':
                expect(file_in, wire_obj.tune_frequency, what);
                break;
            default:
                throw std::runtime_error(what + ": unknown directive");
            }
        }
        if (wire_obj.current_defined && wire_obj.voltage_defined)
        {
            throw std::runtime_error(wire_tag + ": a wire cannot have both a current source ('i') and a voltage source ('V')");
        }
        wire_obj.update_from_trajectory();
        winding_array.push_back(wire_obj);
    }
    return winding_array;
}

std::vector<IPTVisual::wire> read_wires(std::string filename)
{
    std::ifstream file_in(filename);
    if (!file_in)
    {
        throw std::runtime_error("cannot open wire definition file '" + filename + "'");
    }
    return read_wires(file_in);
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

// Two observation formats are accepted:
//   grid (legacy):  x_min x_max x_n / y_min y_max y_n / z_min z_max z_n
//                   optionally preceded by the keyword "grid".
//   point list:     "points" k, followed by k lines of "x y z".
IPTVisual::observations read_observation_definition(std::istream &raw_in)
{
    std::istringstream file_in = strip_comments(raw_in);
    std::string first;
    if (!(file_in >> first))
    {
        throw std::runtime_error("empty observation definition");
    }
    if (first == "points")
    {
        int count;
        expect_count(file_in, count, "the number of observation points");
        IPTVisual::observations obs;
        obs.observation_points = read_points(file_in, count, "observation points");
        obs.observation_values.assign(obs.observation_points.size(), vector_3d_complex());
        return obs;
    }
    double x_min, x_max, y_min, y_max, z_min, z_max;
    int x_num, y_num, z_num;
    if (first == "grid")
    {
        expect(file_in, x_min, "observation grid x range");
    }
    else
    {
        try
        {
            size_t used = 0;
            x_min = std::stod(first, &used);
            if (used != first.size())
            {
                throw std::invalid_argument(first);
            }
        }
        catch (const std::exception &)
        {
            throw std::runtime_error("observation definition must start with a number, 'grid' or 'points', found '" + first + "'");
        }
    }
    expect(file_in, x_max, "observation grid x range");
    expect_count(file_in, x_num, "observation grid x count");
    expect(file_in, y_min, "observation grid y range");
    expect(file_in, y_max, "observation grid y range");
    expect_count(file_in, y_num, "observation grid y count");
    expect(file_in, z_min, "observation grid z range");
    expect(file_in, z_max, "observation grid z range");
    expect_count(file_in, z_num, "observation grid z count");
    return IPTVisual::observations(x_min, x_max, x_num, y_min, y_max, y_num, z_min, z_max, z_num);
}

IPTVisual::observations read_observation_definition(std::string filename)
{
    std::ifstream file_in(filename);
    if (!file_in)
    {
        throw std::runtime_error("cannot open observation definition file '" + filename + "'");
    }
    return read_observation_definition(file_in);
}

void write_observation_values(std::string filename, const IPTVisual::observations &obs)
{
    std::ofstream file_out(filename);
    for (int i = 0; i < obs.observation_points.size(); i++)
    {
        file_out << obs.observation_points[i].x << '\t' << obs.observation_points[i].y << '\t' << obs.observation_points[i].z << '\t' << obs.observation_values[i].x.real() << '\t' << obs.observation_values[i].x.imag() << '\t' << obs.observation_values[i].y.real() << '\t' << obs.observation_values[i].y.imag() << '\t' << obs.observation_values[i].z.real() << '\t' << obs.observation_values[i].z.imag() << std::endl;
    }
    file_out.close();
    return;
}
