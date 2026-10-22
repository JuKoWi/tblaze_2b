//calculate band structure
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <tb/model.h>
#include <tb/select_k.h>
#include <complex>
#include <iostream>

int main(){
	auto grid {generate_mp_grid_frac(3,2,1)};
	std::vector<Eigen::Vector3d> sympoints {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
	std::vector<int> ppseg {4};
	auto path {generate_path_frac(sympoints, ppseg)};
	for (const auto& k : path) std::cout << k.transpose() << "\n";
	for (const auto& k : grid) std::cout << k.transpose() << "\n";
	Eigen::Matrix3d ucell {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
	Eigen::Vector3i nearest {1, 0, 0};
	Hopping hop {1,1, nearest, -2.8};
	std::vector<Hopping> hoppings {hop};
	std::vector<Eigen::Vector3d> pos {{0.0, 0.0, 0.0}};
	TightBindingModel model(ucell, pos, hoppings);
	std::cout << model.get_rec_lat() << std::endl;
}
