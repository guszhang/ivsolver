/*
 * arithmetic.cpp
 * 
 * Copyright (c) 2024 Gus Zhang.
 * 
 * This file provides a set of utility functions and operations for 3D vector and 3x3 matrix arithmetic,
 * including rotation matrix generation, vector operations (addition, dot product, cross product, etc.),
 * and logging utilities for debugging. It also includes helper functions for working with complex numbers.
 * 
 * Functions:
 *  - rAlphaGen, rBetaGen, rGammaGen: Generate rotation angles (alpha, beta, gamma) from a 3D vector.
 *  - mRxGen, mRyGen, mRzGen: Generate rotation matrices for X, Y, Z axes.
 *  - mMul: Multiply two 3x3 matrices or a matrix and a vector.
 *  - mAdd: Add two 3D vectors.
 *  - mDotProduct: Compute the dot product of two 3D vectors.
 *  - mCrossProduct: Compute the cross product of two 3D vectors.
 *  - mRGen: Generate a rotation matrix from a 3D vector.
 *  - mMid: Compute the midpoint between two 3D vectors.
 *  - mVec: Compute the vector from one 3D point to another.
 *  - mLen: Compute the length between two 3D points or the magnitude of a vector.
 *  - mLog: Log the contents of a matrix or vector to the standard output.
 *  - operator*: Multiply a 3D vector by a complex number.
 *  - cRadDeg: Convert radians to degrees.
 *  - cLog, cPrint, cPrintPol: Logging and printing utilities for complex numbers.
 * 
 * Note: Some complex number arithmetic functions are commented out.
 * 
 * Dependencies:
 *  - arithmetic.h
 *  - <cmath>, <iostream>, <complex>, <string>
 * 
 * Usage:
 *  Include this file and its header in your project to perform 3D vector and matrix operations,
 *  as well as basic complex number manipulations and logging.
 */
#include "arithmetic.h"

double rAlphaGen(vector_3d V)
{
	double A = atan(V.y / V.x);
	if (V.x < 0)
		A = A + pi;
	return A;
};

double rBetaGen(vector_3d V)
{
	double MagV = sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
	double PrjV = sqrt(V.x * V.x + V.y * V.y);
	double B = asin(PrjV / MagV);
	if (V.z < 0)
		B = pi - B;
	return B;
};

double rGammaGen(vector_3d V)
{
	double G;
	if (V.z == 0)
	{
		G = -pi / 2;
		if (V.y < 0)
			G = G + pi;
	}
	else
		G = -atan(V.y / V.z);
	if (V.z < 0)
		G = G + pi;
	return G;
};

matrix_33 mRxGen(double gamma)
{
	matrix_33 R = {1, 0, 0, 0, cos(gamma), -sin(gamma), 0, sin(gamma), cos(gamma)};
	return R;
};

matrix_33 mRyGen(double beta)
{
	matrix_33 R = {cos(beta), 0, sin(beta), 0, 1, 0, -sin(beta), 0, cos(beta)};
	return R;
};

matrix_33 mRzGen(double alpha)
{
	matrix_33 R = {cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, 0, 0, 1};
	return R;
};

matrix_33 mMul(matrix_33 a, matrix_33 b)
{
	matrix_33 m = {b.a11 * a.a11 + b.a21 * a.a12 + b.a31 * a.a13, b.a12 * a.a11 + b.a22 * a.a12 + b.a32 * a.a13, b.a13 * a.a11 + b.a23 * a.a12 + b.a33 * a.a13,
				   b.a11 * a.a21 + b.a21 * a.a22 + b.a31 * a.a23, b.a12 * a.a21 + b.a22 * a.a22 + b.a32 * a.a23, b.a13 * a.a21 + b.a23 * a.a22 + b.a33 * a.a23,
				   b.a11 * a.a31 + b.a21 * a.a32 + b.a31 * a.a33, b.a12 * a.a31 + b.a22 * a.a32 + b.a32 * a.a33, b.a13 * a.a31 + b.a23 * a.a32 + b.a33 * a.a33};
	return m;
};

vector_3d mMul(matrix_33 a, vector_3d v)
{
	vector_3d o = {a.a11 * v.x + a.a12 * v.y + a.a13 * v.z,
				   a.a21 * v.x + a.a22 * v.y + a.a23 * v.z,
				   a.a31 * v.x + a.a32 * v.y + a.a33 * v.z};
	return o;
};

vector_3d mAdd(vector_3d a, vector_3d b)
{
	vector_3d c = {a.x + b.x, a.y + b.y, a.z + b.z};
	return c;
}

double mDotProduct(vector_3d a, vector_3d b)
{
	double c = a.x * b.x + a.y * b.y + a.z * b.z;
	return c;
}

vector_3d mCrossProduct(vector_3d a, vector_3d b)
{
	vector_3d c = {a.y * b.z - b.y * a.z, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
	return c;
}

matrix_33 mRGen(vector_3d v)
{
	matrix_33 R;
	if (v.x == 0)
		R = mRxGen(rGammaGen(v));
	else
		R = mMul(mRzGen(rAlphaGen(v)), mRyGen(rBetaGen(v)));
	return R;
};

vector_3d mMid(vector_3d v1, vector_3d v2)
{
	vector_3d v;
	v.x = (v1.x + v2.x) / 2;
	v.y = (v1.y + v2.y) / 2;
	v.z = (v1.z + v2.z) / 2;
	return v;
}

vector_3d mVec(vector_3d v1, vector_3d v2)
{
	vector_3d v;
	v.x = v2.x - v1.x;
	v.y = v2.y - v1.y;
	v.z = v2.z - v1.z;
	return v;
}

double mLen(vector_3d v1, vector_3d v2)
{
	double len;
	len = sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z));
	return len;
}

double mLen(vector_3d v)
{
	double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	return len;
}

void mLog(std::string m, matrix_33 a)
{
	std::cout << m << "=" << std::endl;
	std::cout << a.a11 << " " << a.a12 << " " << a.a13 << std::endl;
	std::cout << a.a21 << " " << a.a22 << " " << a.a23 << std::endl;
	std::cout << a.a31 << " " << a.a32 << " " << a.a33 << std::endl;
};

void mLog(std::string m, vector_3d a)
{
	std::cout << m << "=" << std::endl;
	std::cout << a.x << std::endl;
	std::cout << a.y << std::endl;
	std::cout << a.z << std::endl;
};

vector_3d_complex operator*(const vector_3d a, std::complex<double> b)
{
	return vector_3d_complex(a.x * b, a.y * b, a.z * b);
}

// mvComplex cAdd(mvComplex a, mvComplex b)
// {
// 	mvComplex c = {a.real + b.real, a.imag + b.imag};
// 	return c;
// }

// mvComplex cAdd(double a, mvComplex b)
// {
// 	mvComplex c = {a + b.real, b.imag};
// 	return c;
// }

// mvComplex cAdd(mvComplex a, double b)
// {
// 	mvComplex c = {a.real + b, a.imag};
// 	return c;
// }

// mvComplex cSub(mvComplex a, mvComplex b)
// {
// 	mvComplex c = {a.real - b.real, a.imag - b.imag};
// 	return c;
// }

// mvComplex cSub(double a, mvComplex b)
// {
// 	mvComplex c = {a - b.real, b.imag};
// 	return c;
// }

// mvComplex cSub(mvComplex a, double b)
// {
// 	mvComplex c = {a.real - b, a.imag};
// 	return c;
// }

// mvComplex cMul(mvComplex a, mvComplex b)
// {
// 	mvComplex c = {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
// 	return c;
// }

// mvComplex cMul(double a, mvComplex b)
// {
// 	mvComplex c = {a * b.real, a * b.imag};
// 	return c;
// }

// mvComplex cMul(mvComplex a, double b)
// {
// 	mvComplex c = {a.real * b, a.imag * b};
// 	return c;
// }

// mvComplex cDiv(mvComplex a, mvComplex b)
// {
// 	double d = b.real * b.real + b.imag * b.imag;
// 	mvComplex c = {(a.real * b.real + a.imag * b.imag) / d, (a.imag * b.real - a.real * b.imag) / d};
// 	return c;
// }

// mvComplex cDiv(double a, mvComplex b)
// {
// 	double d = b.real * b.real + b.imag * b.imag;
// 	mvComplex c = {a * b.real / d, -a * b.imag / d};
// 	return c;
// }

// mvComplex cDiv(mvComplex a, double b)
// {
// 	mvComplex c = {a.real / b, a.imag / b};
// 	return c;
// }

// mvComplex cConj(mvComplex a)
// {
// 	mvComplex c = {a.real, -a.imag};
// 	return c;
// }

// double cAbs(mvComplex a)
// {
// 	return sqrt(a.real * a.real + a.imag * a.imag);
// }

// double cAngle(mvComplex a)
// {
// 	if (a.real == 0)
// 	{
// 		if (a.imag > 0)
// 			return pi / 2;
// 		else
// 			return 3 * pi / 2;
// 	}
// 	double d = atan(a.imag / a.real);
// 	if (a.real < 0)
// 		d += pi;
// 	if (d > pi)
// 		d -= 2 * pi;
// 	return d;
// }

double cRadDeg(double a)
{
	return a / pi * 180;
}

void cLog(std::string str, std::complex<double> a)
{
	std::cout << str << "=" << std::endl
			  << a.real();
	if (a.imag() > 0)
		std::cout << "+";
	std::cout << a.imag() << "i" << std::endl;
}

void cPrint(std::complex<double> a)
{
	std::cout << a.real();
	if (a.imag() >= 0)
		std::cout << "+";
	std::cout << a.imag() << "i";
}

void cPrintPol(std::complex<double> a)
{
	std::cout << std::abs(a) * (1 + 1e-10) << " < " << std::arg(a) * (1 + 1e-10) << " (" << cRadDeg(std::arg(a)) * (1 + 1e-10) << "degree)" << std::endl;
}
