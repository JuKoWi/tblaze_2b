#pragma once
#include <fstream>
#include <Eigen/Dense>
#include <tb/select_k.h>

void write_bands_path(KGridQuantity<Eigen::Vector3d> kfrac, KGridQuantity<Eigen::VectorXd> eigvals);
TightBindingModel from_seedname(std::string fname);
std::istringstream read_line(std::ifstream& infile);
void write_dos(std::vector<double> E_au, Eigen::MatrixXd dos);
