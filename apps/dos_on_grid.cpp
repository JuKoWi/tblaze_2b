#include <Eigen/Dense>
#include <tb/model.h>
#include <tb/select_k.h>
#include <tb/read_and_write.h>
#include <tb/utils.h>
#include <iostream>


int main(){
	auto k_grid {gamma_grid(100,100,1)};

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
	std::vector<Eigen::VectorXd> bands {calculate_bands(k_grid, tb, true)};
	std::vector<double> e_range {energy_range(bands)};
	int n_k = int(k_grid.size());
	std::vector<double> dos {calculate_dos(bands, e_range, 0.001)};
	write_dos(e_range, dos);
}
