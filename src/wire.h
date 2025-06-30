// Copyright (c) 2024 Gus Zhang.
// 
// Summary:
//   This header defines the IPTVisual namespace, which contains the 'wire' and 'observations' classes.
//   - The 'wire' class models a 3D wire with geometric and electromagnetic properties, supporting operations
//     such as appending trajectories, geometric transformations, and inductance calculations.
//   - The 'observations' class manages a collection of observation points and their associated values in 3D space,
//     useful for field calculations and visualization.
//
//   Utility macros for numerical precision and iteration limits are also defined.
#pragma once

#define EPS 1e-6
#define MAX_ITER 100

#include <vector>
#include "arithmetic.h"

namespace IPTVisual
{
    class wire
    {
    public:
        wire(double section_radius);
        ~wire();

        double section_radius;
        std::vector<vector_3d> trajectory;
        std::vector<segment> segments;
        std::complex<double> current;
        double total_length;

        void append_any_spiral(const std::vector<vector_3d> key_points);
        void append_any_lines(const std::vector<vector_3d> key_points);
        void append_any_saggy_lines(const std::vector<vector_3d> key_points, const double a);
        void append_any_saggy_lines_with_length(const std::vector<vector_3d> key_points, const std::vector<double> length_array);
        void rotate_x(double angle);
        void rotate_y(double angle);
        void rotate_z(double angle);
        void move_offset(vector_3d offset);
        void align_normal(vector_3d normal);
        void bend_circle_x(double radius);

        void update_from_trajectory();

        double inductance_with(wire w);
        vector_3d_complex b_field_at_point(vector_3d point);

    private:
        double segment_self_inductance(double segment_length);
        void fill_points_line(vector_3d head, vector_3d rear);
        void fill_points_arc(vector_3d head, vector_3d rear);
        void fill_points_saggy_line(vector_3d head, vector_3d rear, const double a);
        double solve_a(const vector_3d head, const vector_3d rear, const double length);
        double solve_rx(const double a, const double dr, const double dz);
        bool model_is_valid;
    };

    class observations
    {
    public:
        observations(double x_min, double x_max, int x_num, double y_min, double y_max, int y_num, double z_min, double z_max, int z_num);
        observations();
        ~observations();
        std::vector<vector_3d> observation_points;
        std::vector<vector_3d_complex> observation_values;
    };
}