#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <vector>
#include <numbers>
#include <tb/model.h>
#include <cmath>
#include <complex>

using namespace std::complex_literals;
using cplx = std::complex<double>;

Eigen::MatrixXcd TightBindingModel::get_Hk(const Eigen::Vector3d& k, bool cell_gauge) const {
    Eigen::MatrixXcd H_k {Eigen::MatrixXcd::Zero(norb(), norb())};
    if (cell_gauge){
        for (const auto& h: H_hops){
            H_k(h.i, h.j) += std::exp(1i* 2.0 * std::numbers::pi * h.R.cast<double>().dot(k)) * h.t;
        }
    }
    else{
        const Eigen::Vector3d rec_lat_cart {B_.transpose() * k};
        for (const auto& h: H_hops){
            const Eigen::Vector3d lat_cart {A_.transpose() * h.R.cast<double>()};
            const Eigen::Vector3d d {lat_cart + tau_[h.j] - tau_[h.i]};
            H_k(h.i, h.j) += std::exp(1i* rec_lat_cart.dot(d)) * h.t;
        }
    }
    double herm_delta {(H_k - H_k.adjoint()).cwiseAbs().maxCoeff()};
    if (herm_delta > 1e-14){
        std::cout << "Hamiltonian not hermitian: ";
        std::cout << herm_delta << '\n';
    }
    return H_k;
}


std::vector<Eigen::VectorXd> calculate_bands(std::vector<Eigen::Vector3d> kfrac, TightBindingModel model, bool cell_gauge){
    std::cout << "Calculate Bands" << std::endl;
    std::vector<Eigen::VectorXd> eigvals {};
    for (const auto& k: kfrac){
        Eigen::MatrixXcd Hk {model.get_Hk(k, cell_gauge)};
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
        eigvals.push_back(es.eigenvalues());
    }
    return eigvals;
}