#pragma once
#include <vector>
#include <Eigen/Dense>

template <typename T>
struct KGridQuantity{
    int nx, ny, nz;
    std::vector<T> data;

    KGridQuantity(int nx, int ny, int nz, T init = T{})
        : nx(nx), ny(ny), nz(nz), data(nx*ny*nz, init)
        {}

    T& operator()(int i, int j, int k){
        return data[(i*ny + j)*nz+ k];
    }
    const T& operator()(int i, int j, int k) const {
        return data[(i*ny + j)*nz+ k];
    }
};

KGridQuantity<Eigen::Vector3d> gamma_grid(int nx, int ny, int nz);
KGridQuantity<Eigen::Vector3d> generate_path_frac(std::vector<Eigen::Vector3d> sym_points, std::vector<int> pntsperseg);
KGridQuantity<Eigen::Vector3d> mp_grid(int nx, int ny, int nz);
