#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numbers>
#include <tb/model.h>
#include <tb/select_k.h>
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

KGridQuantity<Eigen::VectorXd> calculate_bands(KGridQuantity<Eigen::Vector3d>& kfrac, TightBindingModel model, bool cell_gauge){
    std::cout << "Calculate Bands" << std::endl;
    KGridQuantity<Eigen::VectorXd> eigvals(kfrac.nx, kfrac.ny, kfrac.nz);
    #pragma omp parallel for schedule(static)
    for (std::size_t i=0; i<kfrac.data.size(); ++i){
        Eigen::MatrixXcd Hk {model.get_Hk(kfrac.data[i], cell_gauge)};
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hk);
        eigvals.data[i] = es.eigenvalues();
    }
    return eigvals;
}

std::vector<double> get_energy_range(KGridQuantity<Eigen::VectorXd> eigvals, int n_energies){
    std::vector<double> energy_au(n_energies, 0.0);
    double emax = -std::numeric_limits<double>::infinity();
    double emin = std::numeric_limits<double>::infinity();
    for (const auto v: eigvals.data){
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

// density of states per unit cell
Eigen::MatrixXd calculate_gauss_dos(const KGridQuantity<Eigen::VectorXd>& eigvals, const std::vector<double>& erange_au, double broadening_au){
    std::cout<< "Calculate DOS" << '\n';
    int n_k = static_cast<int>(eigvals.data.size());
    int n_bands = eigvals.data[0].size();
    int n_e = static_cast<int>(erange_au.size());
    Eigen::MatrixXd dos {Eigen::MatrixXd::Zero(n_e, n_bands)};
    double norm_const {2.0/(n_k * broadening_au * std::sqrt(2 * std::numbers::pi))};
    std::cout << "Number of k-points: " << '\n';
    std::cout << n_k << '\n';
    #pragma omp parallel for schedule(static)
    for (int i =0; i < n_e; i++){
        for (int j =0; j < n_k; j++){
            for (int k=0; k < n_bands; k++){
                dos(i, k) += norm_const * std::exp(-0.5 * std::pow((erange_au[i] - eigvals.data[j](k))/broadening_au, 2));
            }
        }
    }
    return dos;
}

Eigen::MatrixXd calculate_tetra_dos(const KGridQuantity<Eigen::VectorXd>& eigvals, const std::vector<double>& erange_au, const Eigen::MatrixXd& rec_lat){
    std::cout << "Calculate DOS" << '\n';
    int n_bands = eigvals.data[0].size();
    int n_e = static_cast<int>(erange_au.size());
    Eigen::MatrixXd dos {Eigen::MatrixXd::Zero(n_e, n_bands)};
    std::vector<int> single_points {
        (eigvals.nx == 1) ? 0 : 1,
        (eigvals.ny == 1) ? 0 : 1,
        (eigvals.nz == 1) ? 0 : 1,
    };
    int dimension {0};
    for (auto i: single_points){
        dimension += i;
    }
    if (dimension == 3){
        std::cout << "Applying 3D lattice tetrahedron method" << '\n';
        std::vector<Eigen::Vector3i> corners {
            {0, 0, 0}, // 0
            {1, 0, 0}, // 1
            {0, 1, 0}, // 2
            {1, 1, 0}, // 3
            {0, 0, 1}, // 4
            {1, 0, 1}, // 5
            {0, 1, 1}, // 6
            {1, 1, 1}, // 7
        };
        std::vector<std::vector<int>> tetrahedra {
            {0, 1, 2, 5},
            {0, 2, 4, 5},
            {1, 2, 3, 5},
            {2, 3, 5, 7},
            {2, 5, 6, 7},
            {2, 4, 5, 6},
        };
        double volume_bz {rec_lat.determinant()};
        double volume_tetrahedron {volume_bz /(6.0 * eigvals.nx * eigvals.ny * eigvals.nz)};//assume uniformly spaced grid points
        double prefactor {6 * volume_tetrahedron/ volume_bz};
        #pragma omp parallel for schedule(static)
        for (int ie=0; ie < n_e; ie++){
            for (int i=0; i<eigvals.nx; i++){
                for (int j=0; j<eigvals.ny; j++){
                    for (int k=0; k<eigvals.nz; k++){
                        for (auto t: tetrahedra){
                            for (int b=0; b<n_bands; b++){
                                std::vector<double> e_c {}; // energies at corners
                                for (auto c: t){
                                    e_c.push_back(eigvals((i+corners[c](0)) % eigvals.nx, (j+corners[c](1)) % eigvals.ny, (k+corners[c](2)) % eigvals.nz)(b));
                                }
                                std::ranges::sort(e_c);
                                if (e_c[0]<= erange_au[ie] && erange_au[ie] < e_c[1]){
                                    dos(ie, b) += prefactor * std::pow(erange_au[ie] - e_c[0], 2) / ((e_c[1] - e_c[0]) * (e_c[2] - e_c[0])* (e_c[3] - e_c[0]));
                                }
                                else if (e_c[1] <= erange_au[ie] && erange_au[ie] < e_c[2]){
                                    dos(ie, b) += prefactor * (e_c[1] - e_c[0] + 2 * (erange_au[ie] - e_c[1]) - (e_c[2] - e_c[0] + e_c[3] - e_c[1])*std::pow(erange_au[ie] - e_c[1], 2)/((e_c[2] - e_c[1]) * (e_c[3] - e_c[1])))/((e_c[2] - e_c[0]) * (e_c[3] - e_c[0]));
                                }
                                else if (e_c[2] <= erange_au[ie] && erange_au[ie] < e_c[3]){
                                    dos(ie, b) += prefactor * std::pow(e_c[3] - erange_au[ie], 2) / ((e_c[3] - e_c[0]) * (e_c[3] - e_c[1]) * (e_c[3] - e_c[2]));
                                }
                            }
                        }
                    }
                }
            }
        }
        return dos;
    }
    else if (dimension == 2){
        std::cout << "Applying 2D lattice tetrahedron method" << '\n';
        std::vector<Eigen::Vector3i> corners {
            {0, 0, 0}, // 0
            {1, 0, 0}, // 1
            {0, 1, 0}, // 2
            {1, 1, 0}, // 3
        };
        std::vector<std::vector<int>> triangles {
            {0, 1, 3},
            {0, 2, 3},
        };
        double volume_bz {rec_lat.topLeftCorner<2,2>().determinant()};
        double volume_tetrahedron {volume_bz * 0.5/(eigvals.nx * eigvals.ny)};//assume uniformly spaced grid points
        double prefactor {4 * volume_tetrahedron/ volume_bz};
        #pragma omp parallel for schedule(static)
        for (int ie=0; ie < n_e; ie++){
            for (int i=0; i<eigvals.nx; i++){
                for (int j=0; j<eigvals.ny; j++){
                    for (auto t: triangles){
                        for (int b=0; b<n_bands; b++){
                            std::vector<double> e_vertices {};
                            for (auto c: t){
                                e_vertices.push_back(eigvals((i+corners[c](0)) % eigvals.nx, (j+corners[c](1)) % eigvals.ny, 0)(b));
                            }
                            std::ranges::sort(e_vertices);
                            if (e_vertices[0]<= erange_au[ie] && erange_au[ie] < e_vertices[1]){
                                dos(ie, b) += prefactor * (erange_au[ie] - e_vertices[0])/((e_vertices[1] - e_vertices[0]) * (e_vertices[2] - e_vertices[0]));
                            }
                            else if (e_vertices[1] <= erange_au[ie] && erange_au[ie] < e_vertices[2]){
                                dos(ie, b) += prefactor * (erange_au[ie] - e_vertices[2])/((e_vertices[1] - e_vertices[2]) * (e_vertices[2] - e_vertices[0]));
                            }
                        }
                    }

                }
            }
        }
        return dos;
    }
    else if (dimension == 1){
        std::cout << "Applying 3D lattice tetrahedron method" << '\n';
        std::cout << "Not implemented yet" << '\n';
        return dos;
    }
    else if (dimension == 0){
        std::cout << "Calculating DOS of nonperiodic system" << '\n';
        std::cout << "Not implemented yet" << '\n';
        return dos;
    }
    return dos;
}

