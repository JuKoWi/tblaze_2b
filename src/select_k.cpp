#include <vector>
#include <Eigen/Dense>
#include <tb/select_k.h>

//generate a gamma-centered uniform grid
KGridQuantity<Eigen::Vector3d> gamma_grid(int nx, int ny, int nz){
	KGridQuantity<Eigen::Vector3d> grid(nx, ny, nz, Eigen::Vector3d::Zero());
	for (int i = 0; i<nx; i++){
		for (int j = 0; j<ny; j++){
			for (int k=0; k<nz; k++){
				double kx {static_cast<double>(i) * 1/nx};
				double ky {static_cast<double>(j) * 1/ny};
				double kz {static_cast<double>(k) * 1/nz};
				grid(i,j,k) = Eigen::Vector3d(kx, ky, kz);
			}
		}
	}
	return grid;
}

KGridQuantity<Eigen::Vector3d> generate_path_frac(std::vector<Eigen::Vector3d> sym_points, std::vector<int> pntsperseg){
	if (sym_points.size() != pntsperseg.size()+1 || sym_points.size() <2){
		throw std::invalid_argument("number of symmetry points must be one more than number of segments that connect the points");
	}
	int points_total {1}; //also include end point
	for (auto p: pntsperseg){
		points_total += p;
	}
	KGridQuantity<Eigen::Vector3d> path(points_total, 1, 1, Eigen::Vector3d::Zero());
	int count {0};
	for (std::size_t i=0; i<pntsperseg.size(); i++){
		Eigen::Vector3d start {sym_points[i]};
		Eigen::Vector3d diff = sym_points[i+1] - start;
		int n = pntsperseg[i];
		for (int j=0; j<n; ++j){
			Eigen::Vector3d point {start + static_cast<double>(j) * diff / pntsperseg[i]};
			path(count, 0, 0) = point;
			count +=1;
		}
	}
	path.data[count] = sym_points.back();
	return path;
}

	
KGridQuantity<Eigen::Vector3d> mp_grid(int nx, int ny, int nz){
	KGridQuantity<Eigen::Vector3d> grid(nx, ny, nz, Eigen::Vector3d::Zero());
	for (int i = 0; i<nx; i++){
		for (int j = 0; j<ny; j++){
			for (int k=0; k<nz; k++){
				double kx {(static_cast<double>(i) + (1 - nx)/2.0)/nx};
				double ky {(static_cast<double>(j) + (1 - ny)/2.0)/ny};
				double kz {(static_cast<double>(k) + (1 - nz)/2.0)/nz};
				Eigen::Vector3d kpoint {kx, ky, kz};
				grid(i,j,k) = kpoint;
			}
		}
	}
	return grid;
}