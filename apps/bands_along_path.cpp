//calculate band structure
#include <Eigen/Dense>
#include <tb/model.h>
#include <tb/select_k.h>
#include <tb/read_and_write.h>
#include <tb/utils.h>
#include <iostream>


int main(){
	auto grid {generate_mp_grid_frac(3,2,1)};
	std::vector<Eigen::Vector3d> sympoints {{0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}, {1.0/3.0, 1.0/3.0, 0.0}, {0.0, 0.0, 0.0}};
	std::vector<int> ppseg {100, 100, 100};
	auto path {generate_path_frac(sympoints, ppseg)};

	std::string fname {};
	std::cout << "Waiting for seedname filename" << '\n';
	std::cin >> fname;
	auto tb = from_seedname(fname);
	std::cout << "Lattice in angstrom: " << '\n';
	std::cout << tb.get_lattice()/Constants::ANGSTROM_TO_BOHR << std::endl;
	std::cout << "Reciprocal lattice: " << '\n';
	std::cout << tb.get_rec_lat()  << std::endl;
	std::cout << "Dimension: " << '\n';
	std::cout << tb.norb() << '\n';
	std::cout << "Lattice points: " << '\n';
	for (auto& p: tb.get_lat_points()){
		std::cout << p.transpose() << '\n';
	}
	std::vector<Eigen::VectorXd> bands_w90 {calculate_bands(path, tb, true)};
	write_bands_path(path, bands_w90);
}
