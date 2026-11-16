#pragma once
#include <fstream>
#include <Eigen/Dense>

void write_bands_path(std::vector<Eigen::Vector3d> kfrac, std::vector<Eigen::VectorXd> eigvals);
TightBindingModel from_seedname(std::string fname);
std::istringstream read_line(std::ifstream& infile);
void write_dos(std::vector<double> E_au, std::vector<double> dos);
