#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <vector>
#include <numbers>
#include <tb/model.h>
#include <cmath>
#include <complex>

using namespace std::complex_literals;

Eigen::MatrixXcd TightBindingModel::get_Hk(const Eigen::Vector3d& k) const {
    Eigen::MatrixXcd H_k {Eigen::MatrixXcd::Zero(norb(), norb())};
    for (const auto& h: hops_){
        const Eigen::Vector3d lat_cart {A_ *h.R.cast<double>()};
        const Eigen::Vector3d rec_lat_cart {B_ * k};
        const Eigen::Vector3d d {lat_cart + tau_[h.j] - tau_[h.i]};
        H_k(h.i, h.j) += std::exp(1i* rec_lat_cart.dot(lat_cart + tau_[h.i] - tau_[h.j])) * h.t;
    }
    return H_k;
}

std::vector<Eigen::VectorXd> calculate_bands(std::vector<Eigen::Vector3d> kfrac, TightBindingModel model){
    std::vector<Eigen::VectorXd> eigvals {};
    for (const auto k: kfrac){
        Eigen::MatrixXcd Hk {model.get_Hk(k)};
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
        eigvals.push_back(es.eigenvalues());
    }
    return eigvals;
}