/*
 * Copyright (c) 2024 Gus Zhang.
 *
 * Summary:
 *   This header declares functions for reading and writing various data structures
 *   related to wire and observation definitions, as well as trajectory and matrix data,
 *   for use in the IPTVisual project. The functions provide interfaces for file I/O
 *   operations involving wires, observation definitions, trajectory points, L matrices,
 *   and observation values.
 *
 *   The readers are stream based so that the same parser serves files, standard input
 *   and in-memory strings (used by the WebAssembly build). Malformed input raises
 *   std::runtime_error with a description of the offending wire and directive.
 */
#pragma once

#include <iostream>
#include <fstream>
#include <stdexcept>
#include "wire.h"

std::vector<IPTVisual::wire> read_wires(std::istream &file_in);
std::vector<IPTVisual::wire> read_wires(std::string filename);
IPTVisual::observations read_observation_definition(std::istream &file_in);
IPTVisual::observations read_observation_definition(std::string filename);
void write_trajectory_points(std::string filename, std::vector<vector_3d> trajectory);
void write_l_matrix(std::string filename, int n, double *l_matrix, double *len_array);
void write_observation_values(std::string filename, const IPTVisual::observations &obs);
