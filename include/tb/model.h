#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>

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

class TightBindingModel{
	public:
		TightBindingModel(Eigen::Matrix3d lattice, std::vector<Hopping> hoppings): 
				A_(lattice),
				H_hops(hoppings)
				{
                    B_ = A_.inverse().transpose()* 2 * std::numbers::pi; 
    				std::vector<Eigen::Vector3d> pos(norb(), Eigen::Vector3d::Zero());
					tau_ = pos;
				}

		int norb() const {
			int max = 0;
			for (auto h: H_hops){
				max = std::max(std::max(h.i, h.j), max);
			}
			return max + 1;
		}

		void set_positions(std::vector<Eigen::Vector3d> positions){
			tau_ = positions;
		}
		void set_r_hoppings(std::vector<HoppingVec> r_hoppings){
			r_hops = r_hoppings;
		}
		Eigen::MatrixXcd get_Hk(const Eigen::Vector3d& k, bool cell_gauge) const;
		// void velocity(const Eigen::Vector3d& k, int alpha, Eigen::MatrixXcd& dH);
		static TightBindingModel from_file(const std::string& path);
        Eigen::Matrix3d get_rec_lat(){return B_;};
		Eigen::Matrix3d get_lattice(){return A_;};
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
		std::vector<HoppingVec> r_hops {};
		int dim_basis {};
		Eigen::Matrix3d B_ {}; // rec lattice vectors, vectors as rows
};

std::vector<Eigen::VectorXd> calculate_bands(std::vector<Eigen::Vector3d> kfrac, TightBindingModel model, bool cell_gauge);






