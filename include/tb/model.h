#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <tb/select_k.h>
#include <cassert>

using cplx = std::complex<double>;


struct Hopping {
	int i,j;
	Eigen::Vector3i R;
	cplx t;
};

struct HoppingVec {
	int i,j;
	Eigen::Vector3i R;
	Eigen::Vector3cd t;
};

bool test_rec_lat(Eigen::Matrix3d lat, Eigen::Matrix3d rec_lat);

int norb(const std::vector<Hopping>& hops);

class TightBindingModel{
	public:
		TightBindingModel(Eigen::Matrix3d lattice, std::vector<Hopping> hoppings, bool orthogonal): 
				A_(lattice),
				H_hops(hoppings),
				orthogonal(orthogonal)
				{
                    B_ = A_.inverse().transpose()* 2 * std::numbers::pi; 
    				std::vector<Eigen::Vector3d> pos(norb(hoppings), Eigen::Vector3d::Zero());
					tau_ = pos;
					dim_basis = norb(hoppings);
					assert(test_rec_lat(A_, B_));
				}


		void set_positions(std::vector<Eigen::Vector3d> positions){
			tau_ = positions;
		}
		void set_r_hoppings(std::vector<HoppingVec> r_hoppings){
			r_hops = r_hoppings;
		}
		const std::vector<Hopping>& get_H_hops()const{
			return H_hops;
		}
		const std::vector<Hopping>& get_S_hops()const{
			return S_hops;
		}
		bool is_orthogonal()const{
			return orthogonal;
		}
		void set_overlap(std::vector<Hopping> overlap){
			S_hops = overlap;
			orthogonal = false;
		}
		const int get_dimension()const{
			return dim_basis;
		}
		std::array<Eigen::MatrixXcd,3> get_Hk_deriv(const Eigen::Vector3d& k, bool cell_gauge) const;
		static TightBindingModel from_file(const std::string& path);
        const Eigen::Matrix3d& get_rec_lat()const{
			return B_;
		}
		Eigen::MatrixXcd ft_to_k(const std::vector<Hopping>& hops, const Eigen::Vector3d& k, bool cell_gauge) const;
		const Eigen::Matrix3d& get_lattice(){
			return A_;
		}
		std::vector<Eigen::Vector3i> get_lat_points(){
			std::vector<Eigen::Vector3i> unique_latpoints {};
			for (auto h: H_hops){
				if (std::ranges::find(unique_latpoints, h.R) == unique_latpoints.end()){
					unique_latpoints.push_back(h.R);
				}
			}
			return unique_latpoints;
		}

	private:
		Eigen::Matrix3d A_ {}; // lattice vec, vectors as rows
		std::vector<Eigen::Vector3d> tau_ {}; // orbital centers, repeats for several orbitals at same center
		std::vector<Hopping> H_hops {};
		std::vector<Hopping> S_hops {};
		std::vector<HoppingVec> r_hops {};
		int dim_basis {};
		Eigen::Matrix3d B_ {}; // rec lattice vectors, vectors as rows
		bool orthogonal {};
};

KGridQuantity<Eigen::VectorXd> calculate_bands(KGridQuantity<Eigen::Vector3d>& kfrac, const TightBindingModel& model, bool cell_gauge);
std::vector<double> get_energy_range(KGridQuantity<Eigen::VectorXd> eigvals, int n_energies = 1000);
Eigen::MatrixXd calculate_gauss_dos(const KGridQuantity<Eigen::VectorXd>& eigvals, const std::vector<double>& e_range_au, double broadening_au);
Eigen::MatrixXd calculate_tetra_dos(const KGridQuantity<Eigen::VectorXd>& eigvals, const std::vector<double>& erange_au, const Eigen::MatrixXd& rec_lat);



