//calculate band structure
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <tb/model.h>
#include <tb/select_k.h>
#include <tb/read_and_write.h>
#include <tb/utils.h>
#include <complex>
#include <iostream>
#include <cmath>


int main(){
	auto grid {generate_mp_grid_frac(3,2,1)};
	std::vector<Eigen::Vector3d> sympoints {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {2.0/3.0, 1.0/3.0, 0.0}, {0.0, 0.0, 0.0}};
	std::vector<int> ppseg {100, 100, 100};
	auto path {generate_path_frac(sympoints, ppseg)};
	// for (const auto& k : path) std::cout << k.transpose() << "\n";
	// for (const auto& k : grid) std::cout << k.transpose() << "\n";
	double a_cc {Constants::ANGSTROM_TO_BOHR * 1.42};
	double a {std::sqrt(3) * a_cc};
	Eigen::Matrix3d ucell {{a, 0.0, 0.0}, {a* 0.5, 1.5 * a_cc, 0.0}, {0.0, 0.0, 10.0 * Constants::ANGSTROM_TO_BOHR}};
	Hopping hop1 {0,1, {0, 0, 0}, -2.8};
	Hopping hop2 {1,0, {0, 0, 0}, -2.8};
	Hopping hop3 {0,1, {-1, 0, 0}, -2.8};
	Hopping hop4 {0,1, {0, -1, 0}, -2.8};
	Hopping hop5 {1,0, {1, 0, 0}, -2.8};
	Hopping hop6 {1,0, {0, 1, 0}, -2.8};
	std::vector<Hopping> hoppings {hop1, hop2, hop3, hop4, hop5, hop6};
	std::cout << "Number of hoppings"<< '\n';
	std::cout << hoppings.size() << std::endl;
	std::vector<Eigen::Vector3d> pos {{0.0, 0.0, 0.0}, {a*0.5, a_cc*0.5, 0.0}};
	TightBindingModel model(ucell, pos, hoppings);
	std::cout << "Reciprocal lattice: " << '\n';
	std::cout << model.get_rec_lat() << std::endl;
	std::vector<Eigen::VectorXd> bands {calculate_bands(path, model)};
	write_bands_path(path, bands);
}
