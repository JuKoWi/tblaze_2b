#include <Eigen/Dense>
#include <complex>
#include <vector>

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
				);
		int norb() const {
			return static_cast<int>(tau_.size());
		}
		Eigen::Matrix3d reciprocal() const;

		void hamiltonian(const Eigen::Vector3d& k, Eigen::MatrixXcd& H) const;
		void velocity(const Eigen::Vector3d& k, int alpha, Eigen::MatrixXcd& dH) const;
		static TightBindingModel from_file(const std::string& path);

	private:
		Eigen::Marix3d A_;
		std::vector<Eigen::Vectro3d> tau_;
		std::vector<Hoppint> hops_;
};

