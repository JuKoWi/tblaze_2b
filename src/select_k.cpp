#include <vector>
#include <Eigen/Dense>

// generate Monkhorst Pack grid in fractional coordinates
std::vector<Eigen::Vector3d> generate_mp_grid_frac(int nx, int ny, int nz){
	std::vector<Eigen::Vector3d> grid {};
	for (int i = 0; i<nx; i++){
		for (int j = 0; j<ny; j++){
			for (int k=0; k<nz; k++){
				double kx {static_cast<double>(i) * 1/nx};
				double ky {static_cast<double>(j) * 1/ny};
				double kz {static_cast<double>(k) * 1/nz};
				Eigen::Vector3d kpoint {kx, ky, kz};
				grid.push_back(kpoint);
			}
		}
	}
	return grid;
}

std::vector<Eigen::Vector3d> generate_path_frac(std::vector<Eigen::Vector3d> sym_points, std::vector<int> pntsperseg){
	if (sym_points.size() != pntsperseg.size()+1 || sym_points.size() <2){
		throw std::invalid_argument("number of symmetry points must be one more than number of segments that connect the points");
	}
	std::vector<Eigen::Vector3d> path {};
	for (std::size_t i=0; i<pntsperseg.size(); i++){
		Eigen::Vector3d start {sym_points[i]};
		Eigen::Vector3d diff = sym_points[i+1] - start;
		int n = pntsperseg[i];
		for (int j=0; j<n; ++j){
			Eigen::Vector3d point {start + j * diff / pntsperseg[i]};
			path.push_back(point);
		}
	}
	path.push_back(sym_points.back());
	return path;
}

	
