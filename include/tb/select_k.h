#pragma once
#include <vector>
#include <Eigen/Dense>

std::vector<Eigen::Vector3d> gamma_grid(int nx, int ny, int nz);
std::vector<Eigen::Vector3d> generate_path_frac(std::vector<Eigen::Vector3d> sym_points, std::vector<int> pntsperseg);
