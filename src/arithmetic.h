/*
 * arithmetic.h
 * 
 * Copyright (c) 2024 Gus Zhang.
 * 
 * This header defines fundamental arithmetic structures and operations for 3D vector and matrix algebra,
 * as well as utilities for handling complex numbers and color representations.
 * 
 * Structures:
 *   - vector_3d:      Represents a 3D vector with double precision.
 *   - vector_3d_complex: Represents a 3D vector with complex-valued components.
 *   - matrix_33:      Represents a 3x3 matrix with double precision.
 *   - segment:        Represents a line segment in 3D space, including endpoints, midpoint, direction, and length.
 *   - color:          Represents an RGB color with float components.
 * 
 * Operations:
 *   - Vector and matrix arithmetic (addition, subtraction, multiplication, division).
 *   - Generation of rotation matrices about principal axes.
 *   - Vector algebra (dot product, cross product, length, midpoint, direction vector).
 *   - Conversion utilities and logging/printing functions for vectors, matrices, and complex numbers.
 * 
 * Constants:
 *   - Mathematical and physical constants (pi, e, u0, e0).
 *   - Predefined color values (red, green, blue).
 * 
 * Dependencies:
 *   - <cmath>, <iostream>, <complex>, <string>
 * 
 * Usage:
 *   Include this header in projects requiring 3D vector/matrix operations, complex arithmetic, or color handling.
 */
#pragma once
#ifndef ARITHMETIC_H
#define ARITHMETIC_H
#include <cmath>
#include <iostream>
#include <complex>
#include <string>

// #define SHOW_DEBUG_INFO

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stuctures
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct vector_3d
{
	double x, y, z;
	vector_3d()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	vector_3d(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	vector_3d operator+(const vector_3d a)
	{
		return vector_3d(this->x + a.x, this->y + a.y, this->z + a.z);
	}
	vector_3d operator-(const vector_3d a)
	{
		return vector_3d(this->x - a.x, this->y - a.y, this->z - a.z);
	}
	vector_3d operator*(const double a)
	{
		return vector_3d(this->x * a, this->y * a, this->z * a);
	}
	vector_3d operator/(const double a)
	{
		return vector_3d(this->x / a, this->y / a, this->z / a);
	}

	vector_3d &operator+=(const vector_3d a)
	{
		this->x += a.x;
		this->y += a.y;
		this->z += a.z;
		return *this;
	}
	vector_3d &operator-=(const vector_3d a)
	{
		this->x -= a.x;
		this->y -= a.y;
		this->z -= a.z;
		return *this;
	}
};

struct vector_3d_complex
{
	std::complex<double> x, y, z;
	vector_3d_complex()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	vector_3d_complex(std::complex<double> x, std::complex<double> y, std::complex<double> z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	vector_3d_complex(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	vector_3d_complex(vector_3d a)
	{
		this->x = a.x;
		this->y = a.y;
		this->z = a.z;
	}
	vector_3d_complex operator+(const vector_3d_complex a)
	{
		return vector_3d_complex(this->x + a.x, this->y + a.y, this->z + a.z);
	}
	vector_3d_complex operator-(const vector_3d_complex a)
	{
		return vector_3d_complex(this->x - a.x, this->y - a.y, this->z - a.z);
	}
	vector_3d_complex operator*(const std::complex<double> a)
	{
		return vector_3d_complex(this->x * a, this->y * a, this->z * a);
	}
	vector_3d_complex operator/(const std::complex<double> a)
	{
		return vector_3d_complex(this->x / a, this->y / a, this->z / a);
	}
	vector_3d_complex operator+(const std::complex<double> a)
	{
		return vector_3d_complex(this->x + a, this->y + a, this->z + a);
	}
	vector_3d_complex operator-(const std::complex<double> a)
	{
		return vector_3d_complex(this->x - a, this->y - a, this->z - a);
	}
	vector_3d_complex operator*(const vector_3d a)
	{
		return vector_3d_complex(this->x * a.x, this->y * a.y, this->z * a.z);
	}
};

vector_3d_complex operator*(const vector_3d a, std::complex<double> b);


struct matrix_33
{
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	matrix_33()
	{
		this->a11 = 0;
		this->a12 = 0;
		this->a13 = 0;
		this->a21 = 0;
		this->a22 = 0;
		this->a23 = 0;
		this->a31 = 0;
		this->a32 = 0;
		this->a33 = 0;
	}
	matrix_33(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33)
	{
		this->a11 = a11;
		this->a12 = a12;
		this->a13 = a13;
		this->a21 = a21;
		this->a22 = a22;
		this->a23 = a23;
		this->a31 = a31;
		this->a32 = a32;
		this->a33 = a33;
	}
};

struct segment
{
	vector_3d head, rear, mid, vec;
	double length;
};

struct color
{
	float r, g, b;
	color()
	{
		this->r = 0;
		this->g = 0;
		this->b = 0;
	}
	color(float r, float g, float b)
	{
		this->r = r;
		this->g = g;
		this->b = b;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Operations
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

double rAlphaGen(vector_3d);
double rBetaGen(vector_3d);
double rGammaGen(vector_3d);

matrix_33 mRxGen(double);
matrix_33 mRyGen(double);
matrix_33 mRzGen(double);

matrix_33 mMul(matrix_33, matrix_33);
vector_3d mMul(matrix_33, vector_3d);
vector_3d mAdd(vector_3d, vector_3d);
double mDotProduct(vector_3d, vector_3d);
vector_3d mCrossProduct(vector_3d, vector_3d);

vector_3d mMid(vector_3d, vector_3d);
vector_3d mVec(vector_3d, vector_3d);
double mLen(vector_3d, vector_3d);
double mLen(vector_3d);

matrix_33 mRGen(vector_3d);

// // mvComplex functions removed - use cpp's complex library
// mvComplex cAdd(mvComplex, mvComplex);
// mvComplex cAdd(double, mvComplex);
// mvComplex cAdd(mvComplex, double);
// mvComplex cSub(mvComplex, mvComplex);
// mvComplex cSub(double, mvComplex);
// mvComplex cSub(mvComplex, double);
// mvComplex cMul(mvComplex, mvComplex);
// mvComplex cMul(double, mvComplex);
// mvComplex cMul(mvComplex, double);
// mvComplex cDiv(mvComplex, mvComplex);
// mvComplex cDiv(double, mvComplex);
// mvComplex cDiv(mvComplex, double);
// mvComplex cConj(mvComplex);
// double cAbs(mvComplex);
// double cAngle(mvComplex);
double cRadDeg(double);

void mLog(std::string, matrix_33);
void mLog(std::string, vector_3d);
void cLog(std::string, std::complex<double>);
void cPrint(std::complex<double>);
void cPrintPol(std::complex<double>);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constants
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
const double pi = 3.14159265358979323846264338327950288;
const double u00 = 1e-7;
const double u0 = 4 * pi * u00;
const double e = 2.718281828459045235360287471352662498;
const double e0 = 8.854187817e-12;

const color red = {1, 0, 0};
const color green = {0, 1, 0};
const color blue = {0, 0, 1};

#endif