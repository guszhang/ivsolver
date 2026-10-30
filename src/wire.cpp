/*
 * Copyright (c) 2024 Gus Zhang.
 *
 * Summary:
 *   This file implements the IPTVisual::wire class and related functions for modeling, manipulating,
 *   and analyzing 3D wire trajectories and their electromagnetic properties. The wire class supports
 *   generating wire paths as straight lines, arcs, and sagging (catenary) lines, as well as
 *   calculating self and mutual inductance, magnetic field at observation points, and geometric
 *   transformations such as rotation, translation, and bending. The IPTVisual::observations class
 *   provides a grid of observation points for field calculations.
 *
 *   Key features:
 *     - Generation of wire trajectories using straight, arc, and catenary segments.
 *     - Calculation of segment self-inductance and mutual inductance between wires.
 *     - Magnetic field computation at arbitrary points in space.
 *     - Support for geometric transformations (rotation, translation, bending).
 *     - Observation grid management for field visualization or analysis.
 */
#include "wire.h"
#include <algorithm>

namespace
{
vector_3d normalized(vector_3d v)
{
    double length = mLen(v);
    return vector_3d(v.x / length, v.y / length, v.z / length);
}

void make_plane_basis(vector_3d normal, vector_3d x_direction, vector_3d &x_axis, vector_3d &y_axis)
{
    normal = normalized(normal);
    double projection = mDotProduct(x_direction, normal);
    x_direction.x -= projection * normal.x;
    x_direction.y -= projection * normal.y;
    x_direction.z -= projection * normal.z;
    x_axis = normalized(x_direction);
    y_axis = normalized(mCrossProduct(normal, x_axis));
}
}

// fill segments from <head>(non-inclusive) to <rear>(inclusive)
void IPTVisual::wire::fill_points_line(vector_3d head, vector_3d rear)
{
    double stroke_length = mLen(head, rear);
    int segment_count = std::floor(stroke_length / this->section_radius / 2.0);
    vector_3d increment;
    increment.x = ((rear.x - head.x) / stroke_length * this->section_radius * 2);
    increment.y = ((rear.y - head.y) / stroke_length * this->section_radius * 2);
    increment.z = ((rear.z - head.z) / stroke_length * this->section_radius * 2);
    vector_3d working_point = vector_3d(head);
    for (int i = 0; i < segment_count - 1; i++)
    {
        working_point.x += increment.x;
        working_point.y += increment.y;
        working_point.z += increment.z;
        this->trajectory.push_back(working_point);
    }
    this->trajectory.push_back(rear);
}

void IPTVisual::wire::fill_points_arc(vector_3d head, vector_3d rear)
{
    double r_head = sqrt(head.x * head.x + head.y * head.y);
    double r_rear = sqrt(rear.x * rear.x + rear.y * rear.y);
    double t_head = atan2(head.y, head.x);
    double t_rear = atan2(rear.y, rear.x);
    // atan2 returns (-pi, pi]; shift into [0, 2*pi) so that the arc below always
    // sweeps in the increasing-angle direction rather than taking the short route.
    if (t_head < 0)
    {
        t_head += 2 * pi;
    }
    if (t_rear < 0)
    {
        t_rear += 2 * pi;
    }
    if (t_head >= t_rear)
    {
        t_rear += 2 * pi;
    }
    int sc = (int)floor(sqrt((r_head - r_rear) * (r_head - r_rear) + (r_head * (t_rear - t_head)) * (r_head * (t_rear - t_head)) + (head.z - rear.z) * (head.z - rear.z)) / (2 * this->section_radius));
    if (sc < 1)
    {
        sc = 1;
    }
    double dt = (t_rear - t_head) / sc;
    double dz = (rear.z - head.z) / sc;
    double dr = (r_rear - r_head) / sc;
    vector_3d working_point = vector_3d(head);
    for (int i = 1; i < sc; i++)
    {
        working_point.x = (r_head + i * dr) * cos(t_head + i * dt);
        working_point.y = (r_head + i * dr) * sin(t_head + i * dt);
        working_point.z += dz;
        this->trajectory.push_back(working_point);
    }
    this->trajectory.push_back(rear);
}

// fill segments from <head>(non-inclusive) to <rear>(inclusive) along an arc that
// rotates around the z-direction axis passing through (axis_x, axis_y), rather than
// around the origin as in fill_points_arc. <direction> is +1 for the increasing-angle
// direction (the same growth direction as the spiral segments) and -1 for the opposite.
void IPTVisual::wire::fill_points_arc_about_axis(vector_3d head, vector_3d rear, const double axis_x, const double axis_y, const int direction)
{
    double head_x = head.x - axis_x;
    double head_y = head.y - axis_y;
    double rear_x = rear.x - axis_x;
    double rear_y = rear.y - axis_y;
    double r_head = sqrt(head_x * head_x + head_y * head_y);
    double r_rear = sqrt(rear_x * rear_x + rear_y * rear_y);
    double t_head = atan2(head_y, head_x);
    double t_rear = atan2(rear_y, rear_x);
    if (direction >= 0)
    {
        while (t_rear <= t_head)
        {
            t_rear += 2 * pi;
        }
    }
    else
    {
        while (t_rear >= t_head)
        {
            t_rear -= 2 * pi;
        }
    }
    int sc = (int)floor(sqrt((r_head - r_rear) * (r_head - r_rear) + (r_head * (t_rear - t_head)) * (r_head * (t_rear - t_head)) + (head.z - rear.z) * (head.z - rear.z)) / (2 * this->section_radius));
    if (sc < 1)
    {
        sc = 1;
    }
    double dt = (t_rear - t_head) / sc;
    double dz = (rear.z - head.z) / sc;
    double dr = (r_rear - r_head) / sc;
    vector_3d working_point = vector_3d(head);
    for (int i = 1; i < sc; i++)
    {
        working_point.x = axis_x + (r_head + i * dr) * cos(t_head + i * dt);
        working_point.y = axis_y + (r_head + i * dr) * sin(t_head + i * dt);
        working_point.z = head.z + i * dz;
        this->trajectory.push_back(working_point);
    }
    this->trajectory.push_back(rear);
}

double fx(double x, double a, double dr, double dz)
{
    return a * cosh((dr - x) / a) - a * cosh(-x / a) - dz;
}

double dfxdx(double x, double a, double dr)
{
    return -sinh((dr - x) / a) - sinh(x / a);
}

double dfxda(double x, double a, double dr)
{
    return (a * (cosh((dr - x) / a) - cosh(x / a)) - (dr - x) * sinh((dr - x) / a) + x * sinh(x / a)) / a;
}

double gx(double x, double a, double dr, double l)
{
    return a * sinh((dr - x) / a) - a * sinh(-x / a) - l;
}

double dgxdx(double x, double a, double dr)
{
    return -cosh((dr - x) / a) + cosh(x / a);
}

double dgxda(double x, double a, double dr)
{
    return (a * (sinh((dr - x) / a) + sinh(x / a)) - (dr - x) * cosh((dr - x) / a) - x * cosh(x / a)) / a;
}

double IPTVisual::wire::solve_rx(const double a, const double dr, const double dz)
{
    double x = 0;
    for (int i = 0; i < MAX_ITER; i++)
    {
        double dx = -fx(x, a, dr, dz) / dfxdx(x, a, dr);
        x += dx;
        if (fabs(dx) < EPS)
        {
            return x;
        }
    }
    return x;
}

double IPTVisual::wire::solve_a(const vector_3d head, const vector_3d rear, const double length)
{
    double a = 0.5;
    double dr = sqrt((head.x - rear.x) * (head.x - rear.x) + (head.y - rear.y) * (head.y - rear.y));
    double r0 = dr / 2;
    double dz = rear.z - head.z;
    double j11, j12, j21, j22, ji11, ji12, ji21, ji22, deti, f1, f2, da;
    for (int i = 0; i < MAX_ITER; i++)
    {
        // Inverse Jacobian Matrix
        j11 = dfxdx(r0, a, dr);
        j12 = dfxda(r0, a, dr);
        j21 = dgxdx(r0, a, dr);
        j22 = dgxda(r0, a, dr);
        deti = 1 / (j11 * j22 - j12 * j21);
        ji11 = j22 * deti;
        ji12 = -j12 * deti;
        ji21 = -j21 * deti;
        ji22 = j11 * deti;
        f1 = fx(r0, a, dr, dz);
        f2 = gx(r0, a, dr, length);
        r0 -= ji11 * f1 + ji12 * f2;
        da = ji21 * f1 + ji22 * f2;
        a -= da;
        if (fabs(da) < EPS)
        {
            return a;
        }
    }
    return a;
}

void IPTVisual::wire::fill_points_saggy_line(vector_3d head, vector_3d rear, const double a)
{
    // head and rear are not vertically aligned, checked in the caller.
    // interpolate and estimate total length.
    double dx = rear.x - head.x;
    double dy = rear.y - head.y;
    double dz = rear.z - head.z;
    double dr = sqrt(dx * dx + dy * dy);
    double theta = atan2(dy, dx);
    // saggy line equation
    // y = a*cosh((x-x0)/a)+y0
    double r00 = solve_rx(a, dr, dz);
    double z00 = -a * cosh(-r00 / a);
    double r_current = 0;
    while (1)
    {
        double slope = sinh((r_current - r00) / a);
        // seg = sqrt(dx*dx+dy*dy) = sqrt(dx*dx(1+slope*slope)))
        double drr = sqrt(this->section_radius * this->section_radius * 4.0 / (1 + slope * slope));
        r_current += drr;
        if (r_current >= dr)
        {
            break;
        }
        vector_3d working_point;
        working_point.x = head.x + r_current * cos(theta);
        working_point.y = head.y + r_current * sin(theta);
        working_point.z = head.z + a * cosh((r_current - r00) / a) + z00;
        // double distance = mLen(this->trajectory.back(), working_point);
        // std::cout << r_current << " " << distance << std::endl;
        this->trajectory.push_back(working_point);
    }
    this->trajectory.push_back(rear);
}

void IPTVisual::wire::append_any_lines(const std::vector<vector_3d> key_points)
{
    if (key_points.size() == 0)
    {
        return;
    }
    this->model_is_valid = false;
    if (this->trajectory.size() == 0)
    {
        this->trajectory.push_back(key_points[0]);
    }
    else
    {
        this->fill_points_line(this->trajectory.back(), key_points[0]);
    }

    for (int i = 1; i < key_points.size(); i++)
    {
        this->fill_points_line(key_points[i - 1], key_points[i]);
    }
}

void IPTVisual::wire::append_any_spiral(const std::vector<vector_3d> key_points)
{
    if (key_points.size() == 0)
    {
        return;
    }
    this->model_is_valid = false;
    if (this->trajectory.size() == 0)
    {
        this->trajectory.push_back(key_points[0]);
    }
    else
    {
        this->fill_points_arc(this->trajectory.back(), key_points[0]);
    }

    for (int i = 1; i < key_points.size(); i++)
    {
        this->fill_points_arc(key_points[i - 1], key_points[i]);
    }
}

void IPTVisual::wire::append_any_arc(const std::vector<vector_3d> key_points, const double axis_x, const double axis_y, const int direction)
{
    if (key_points.size() == 0)
    {
        return;
    }
    this->model_is_valid = false;
    if (this->trajectory.size() == 0)
    {
        this->trajectory.push_back(key_points[0]);
    }
    else
    {
        this->fill_points_arc_about_axis(this->trajectory.back(), key_points[0], axis_x, axis_y, direction);
    }

    for (int i = 1; i < key_points.size(); i++)
    {
        this->fill_points_arc_about_axis(key_points[i - 1], key_points[i], axis_x, axis_y, direction);
    }
}

void IPTVisual::wire::append_any_saggy_lines(const std::vector<vector_3d> key_points, const double a)
{
    if (key_points.size() == 0)
    {
        return;
    }
    this->model_is_valid = false;
    if (this->trajectory.size() == 0)
    {
        this->trajectory.push_back(key_points[0]);
    }
    else
    {
        if (this->trajectory.back().x == key_points[0].x && this->trajectory.back().y == key_points[0].y)
        {
            this->fill_points_line(this->trajectory.back(), key_points[0]);
        }
        else
        {
            this->fill_points_saggy_line(this->trajectory.back(), key_points[0], a);
        }
    }

    for (int i = 1; i < key_points.size(); i++)
    {
        if (key_points[i - 1].x == key_points[i].x && key_points[i - 1].y == key_points[i].y)
        {
            this->fill_points_line(key_points[i - 1], key_points[i]);
        }
        else
        {
            this->fill_points_saggy_line(key_points[i - 1], key_points[i], a);
        }
    }
}

void IPTVisual::wire::append_any_saggy_lines_with_length(const std::vector<vector_3d> key_points, const std::vector<double> length_array)
{
    if (key_points.size() == 0)
    {
        return;
    }
    this->model_is_valid = false;
    if (this->trajectory.size() == 0)
    {
        this->trajectory.push_back(key_points[0]);
    }
    else
    {
        if (this->trajectory.back().x == key_points[0].x && this->trajectory.back().y == key_points[0].y)
        {
            this->fill_points_line(this->trajectory.back(), key_points[0]);
        }
        else
        {
            double a = solve_a(this->trajectory.back(), key_points[0], length_array[0]);
            this->fill_points_saggy_line(this->trajectory.back(), key_points[0], a);
        }
    }

    for (int i = 1; i < key_points.size(); i++)
    {
        if (key_points[i - 1].x == key_points[i].x && key_points[i - 1].y == key_points[i].y)
        {
            this->fill_points_line(key_points[i - 1], key_points[i]);
        }
        else
        {
            double a = solve_a(key_points[i - 1], key_points[i], length_array[i]);
            this->fill_points_saggy_line(key_points[i - 1], key_points[i], a);
        }
    }
}

void IPTVisual::wire::append_ellipse(vector_3d center, vector_3d normal, vector_3d x_direction, double major_axis_length, double minor_axis_length, double start_angle, double end_angle)
{
    this->model_is_valid = false;
    vector_3d x_axis, y_axis;
    make_plane_basis(normal, x_direction, x_axis, y_axis);

    const double a = major_axis_length;
    const double b = minor_axis_length;
    const double theta0 = start_angle * pi / 180.0;
    const double theta1 = end_angle * pi / 180.0;
    const double theta_span = theta1 - theta0;
    auto point_at = [&](double theta)
    {
        return vector_3d(center.x + a * cos(theta) * x_axis.x + b * sin(theta) * y_axis.x,
                         center.y + a * cos(theta) * x_axis.y + b * sin(theta) * y_axis.y,
                         center.z + a * cos(theta) * x_axis.z + b * sin(theta) * y_axis.z);
    };

    vector_3d start = point_at(theta0);
    if (this->trajectory.empty())
    {
        this->trajectory.push_back(start);
    }
    else if (mLen(this->trajectory.back(), start) > 0)
    {
        this->fill_points_line(this->trajectory.back(), start);
    }

    if (theta_span == 0)
    {
        return;
    }
    double max_radius = std::max(a, b);
    double estimated_path = sqrt((a - b) * (a - b) + (max_radius * theta_span) * (max_radius * theta_span));
    int steps = std::max(1, static_cast<int>(floor(estimated_path / (2.0 * this->section_radius))));
    for (int i = 1; i <= steps; ++i)
    {
        this->trajectory.push_back(point_at(theta0 + theta_span * i / steps));
    }
}

void IPTVisual::wire::append_parabola(vector_3d origin, vector_3d normal, vector_3d x_direction, double a, double b, double c, double x1, double x2)
{
    this->model_is_valid = false;
    vector_3d x_axis, y_axis;
    make_plane_basis(normal, x_direction, x_axis, y_axis);
    auto point_at = [&](double x)
    {
        double y = a * x * x + b * x + c;
        return vector_3d(origin.x + x * x_axis.x + y * y_axis.x,
                         origin.y + x * x_axis.y + y * y_axis.y,
                         origin.z + x * x_axis.z + y * y_axis.z);
    };

    vector_3d start = point_at(x1);
    if (this->trajectory.empty())
    {
        this->trajectory.push_back(start);
    }
    else if (mLen(this->trajectory.back(), start) > 0)
    {
        this->fill_points_line(this->trajectory.back(), start);
    }

    if (x1 == x2)
    {
        return;
    }
    double max_slope = std::max(fabs(2 * a * x1 + b), fabs(2 * a * x2 + b));
    double estimated_length = fabs(x2 - x1) * sqrt(1 + max_slope * max_slope);
    int steps = std::max(1, static_cast<int>(ceil(estimated_length / (2.0 * this->section_radius))));
    for (int i = 1; i <= steps; ++i)
    {
        double x = x1 + (x2 - x1) * i / steps;
        this->trajectory.push_back(point_at(x));
    }
}

IPTVisual::wire::wire(double section_radius)
{
    this->section_radius = section_radius;
    this->trajectory.clear();
    this->model_is_valid = false;
}

IPTVisual::wire::~wire()
{
}

double IPTVisual::wire::segment_self_inductance(double segment_length)
{
    double l = segment_length; // alias, assume the compiler will optimise this automatically.
    double a = this->section_radius;
    return 2 * u00 * (l * log((l + sqrt(l * l + a * a)) / a) - sqrt(l * l + a * a) + l / 4 + a);
}

void IPTVisual::wire::update_from_trajectory()
{
    this->segments.clear();
    this->total_length = 0;
    if (this->trajectory.size() < 2)
    {
        this->model_is_valid = true;
        return;
    }
    for (int i = 0; i < this->trajectory.size() - 1; i++)
    {
        segment seg;
        seg.head = this->trajectory[i];
        seg.rear = this->trajectory[i + 1];
        seg.mid = mMid(this->trajectory[i], this->trajectory[i + 1]);
        seg.length = mLen(this->trajectory[i], this->trajectory[i + 1]);
        seg.vec = mVec(this->trajectory[i], this->trajectory[i + 1]);
        this->segments.push_back(seg);
        this->total_length += seg.length;
    }
    this->model_is_valid = true;
}

double IPTVisual::wire::inductance_with(const wire &w)
{
    if (!this->model_is_valid)
    {
        this->update_from_trajectory();
    }
    double L = 0;
    for (int i = 0; i < this->segments.size(); i++)
        for (int j = 0; j < w.segments.size(); j++)
            if ((this->segments[i].head.x == w.segments[j].head.x) &&
                (this->segments[i].head.y == w.segments[j].head.y) &&
                (this->segments[i].head.z == w.segments[j].head.z) &&
                (this->segments[i].rear.x == w.segments[j].rear.x) &&
                (this->segments[i].rear.y == w.segments[j].rear.y) &&
                (this->segments[i].rear.z == w.segments[j].rear.z))
                L += segment_self_inductance(this->segments[i].length);
            else
                L += u00 * mDotProduct(this->segments[i].vec, w.segments[j].vec) / mLen(this->segments[i].mid, w.segments[j].mid);
    return L;
}

void IPTVisual::wire::rotate_x(double angle)
{
    this->model_is_valid = false;
    matrix_33 Rx = mRxGen(angle);
    for (auto v = this->trajectory.begin(); v < this->trajectory.end(); v++)
    {
        *v = mMul(Rx, *v);
    }
}
void IPTVisual::wire::rotate_y(double angle)
{
    this->model_is_valid = false;
    matrix_33 Ry = mRyGen(angle);
    for (auto v = this->trajectory.begin(); v < this->trajectory.end(); v++)
    {
        *v = mMul(Ry, *v);
    }
}
void IPTVisual::wire::rotate_z(double angle)
{
    this->model_is_valid = false;
    matrix_33 Rz = mRzGen(angle);
    for (auto v = this->trajectory.begin(); v < this->trajectory.end(); v++)
    {
        *v = mMul(Rz, *v);
    }
}
void IPTVisual::wire::move_offset(vector_3d offset)
{
    this->model_is_valid = false;
    for (auto v = this->trajectory.begin(); v < this->trajectory.end(); v++)
    {
        v->x += offset.x;
        v->y += offset.y;
        v->z += offset.z;
    }
}
void IPTVisual::wire::align_normal(vector_3d normal)
{
    this->model_is_valid = false;
    matrix_33 R = mRGen(normal);
    for (auto v = this->trajectory.begin(); v < this->trajectory.end(); v++)
    {
        *v = mMul(R, *v);
    }
}
void IPTVisual::wire::bend_circle_x(double radius)
{
    this->model_is_valid = false;
    for (auto v = this->trajectory.begin(); v < this->trajectory.end(); v++)
    {
        double angle = pi / 2 - v->y / radius;
        double r = radius + v->z;
        double new_y = cos(angle) * r;
        double new_z = sin(angle) * r - radius;
        v->y = new_y;
        v->z = new_z;
    }
}

vector_3d_complex IPTVisual::wire::b_field_at_point(vector_3d point)
{
    vector_3d B_scalar = vector_3d(0, 0, 0), dB;
    vector_3d_complex B;
    double r;
    for (int i = 0; i < this->segments.size(); i++)
    {
        dB = mCrossProduct(this->segments[i].vec, mVec(this->segments[i].mid, point));
        r = mLen(this->segments[i].mid, point);
        B_scalar.x += dB.x / (r * r * r) * 1e-7;
        B_scalar.y += dB.y / (r * r * r) * 1e-7;
        B_scalar.z += dB.z / (r * r * r) * 1e-7;
    }
    B = B_scalar * this->current;
    return B;
}

IPTVisual::observations::observations(double x_min, double x_max, int x_num, double y_min, double y_max, int y_num, double z_min, double z_max, int z_num)
{
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;
    double z_range = z_max - z_min;
    double x_inc = x_range / (x_num - 1);
    if (x_num <= 1)
    {
        x_inc = 0;
    }
    double y_inc = y_range / (y_num - 1);
    if (y_num <= 1)
    {
        y_inc = 0;
    }
    double z_inc = z_range / (z_num - 1);
    if (z_num <= 1)
    {
        z_inc = 0;
    }
    for (int i = 0; i < x_num; i++)
    {
        for (int j = 0; j < y_num; j++)
        {
            for (int k = 0; k < z_num; k++)
            {
                vector_3d point;
                point.x = x_min + i * x_inc;
                point.y = y_min + j * y_inc;
                point.z = z_min + k * z_inc;
                this->observation_points.push_back(point);
            }
        }
    }
    this->observation_values.clear();
    this->observation_values.resize(this->observation_points.size(), vector_3d_complex());
}

IPTVisual::observations::observations()
{
}

IPTVisual::observations::~observations()
{
}
