#pragma once
#include <Eigen/Dense>

using cplx = std::complex<double>;

struct Hopping {
	int i,j;
	Eigen::Vector3i R;
	cplx t;
};

class TightBindingModel{
	public:
		TightBindingModel(Eigen::Matrix3d lattice,
				std::vector<Eigen::Vector3d> positions,
				std::vector<Hopping> hoppings
				)
				: A_(lattice),
				hops_(hoppings),
				tau_(positions)
				{
					dim_basis = norb();
                    B_ = A_.inverse().transpose()* 2 * std::numbers::pi; 
				}

		int norb() const {
			return static_cast<int>(tau_.size());
		}

		void hamiltonian(const Eigen::Vector3d& k, Eigen::MatrixXcd& H) const;
		void velocity(const Eigen::Vector3d& k, int alpha, Eigen::MatrixXcd& dH) const;
		static TightBindingModel from_file(const std::string& path);
        Eigen::Matrix3d get_rec_lat(){return B_;};

	private:
		Eigen::Matrix3d A_ {}; // lattice vec
		std::vector<Eigen::Vector3d> tau_ {}; // orbital centers, repeats for several orbitals at same center
		std::vector<Hopping> hops_ {};
		int dim_basis {};
		Eigen::Matrix3d B_ {}; // rec lattice vectors
};