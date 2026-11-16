#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numbers>
#include <tb/model.h>
#include <cmath>
#include <complex>
#include <limits>

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

std::array<Eigen::MatrixXcd,3> TightBindingModel::get_Hk_deriv(const Eigen::Vector3d& k, bool cell_gauge) const {
    std::array<Eigen::MatrixXcd, 3> gradHk {};
    if (cell_gauge){
        for (int i=0; i<3; ++i){
            for (const auto& h: H_hops){
                const Eigen::Vector3d lat_cart {A_.transpose() * h.R.cast<double>()};
                gradHk[i](h.i, h.j) += 1i * lat_cart(i) * std::exp(1i * 2.0 * std::numbers::pi * h.R.cast<double>().dot(k)) * h.t;
            }
        }
    }
    else{
        const Eigen::Vector3d rec_lat_cart {B_.transpose() * k};
        for (int i=0; i<3; ++i){
            for (const auto& h: H_hops){
                const Eigen::Vector3d lat_cart {A_.transpose() * h.R.cast<double>()};
                const Eigen::Vector3d d {lat_cart + tau_[h.j] - tau_[h.i]};
                gradHk[i](h.i, h.j) += 1i * d(i) * std::exp(1i* rec_lat_cart.dot(d)) * h.t;
            }
        }
    }
    return gradHk;
}

bool test_rec_lat(Eigen::Matrix3d lat, Eigen::Matrix3d rec_lat){
    Eigen::Matrix3d result {(lat * rec_lat.transpose())};
    return result.isApprox(2 * std::numbers::pi * Eigen::Matrix3d::Identity());
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

std::vector<double> energy_range(std::vector<Eigen::VectorXd> eigvals, int n_energies){
    std::vector<double> energy_au(n_energies, 0.0);
    double emax = -std::numeric_limits<double>::infinity();
    double emin = std::numeric_limits<double>::infinity();
    for (const auto v: eigvals){
        emin = std::min(emin, v.minCoeff());
        emax = std::max(emax, v.maxCoeff());
    }
    emin -= 0.1;
    emax += 0.1;
    double de = (emax - emin)/(n_energies-1);
    for (int i =0; i < n_energies; i++){
        energy_au[i] = emin +i * de;
    }
    return energy_au;
}

std::vector<double> calculate_dos(const std::vector<Eigen::VectorXd>& eigvals, const std::vector<double>& energy_au, double broadening_au){
    std::cout<< "Calculate DOS" << std::endl;
    int n_k = static_cast<int>(eigvals.size());
    int n_bands = eigvals[0].size();
    std::vector<double> dos(energy_au.size(), 0.0);
    double norm_const {2.0/(n_k * broadening_au * std::sqrt(2 * std::numbers::pi))};
    for (std::size_t i =0; i < energy_au.size(); i++){
        for (int j =0; j < n_k; j++){
            for (int k=0; k < n_bands; k++){
                dos[i] += norm_const * std::exp(-0.5 * std::pow((energy_au[i] - eigvals[j](k))/broadening_au, 2));
            }
        }
    }
    return dos;
}
