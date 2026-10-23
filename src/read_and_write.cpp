#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <tb/model.h>
#include <string>

void write_bands_path(std::vector<Eigen::Vector3d> kfrac, std::vector<Eigen::VectorXd> eigvals){
    if (kfrac.size() != eigvals.size()){
        throw std::invalid_argument("Only one set of eigenvalues per k-point allowed");
    }
    std::ofstream resultfile("bands.txt");
    if (resultfile.is_open()){
        for (int i; i<kfrac.size(); i++){
            resultfile << kfrac[i].transpose() << '\t' << eigvals[i].transpose() << '\n';
        } 
        resultfile.close();
    }
    std::cout << "Finished writing bandfile" << '\n';
}

// TightBindingModel from_seedname(std::string fname){
//     Eigen::Matrix3d lattice {};
//     std::ifstream inputfile(fname);
//     if (!inputfile.is_open()){
//         throw std::invalid_argument("Loading seedname file failed");
//     }
//     std::string line {};
//     std::getline(inputfile, line);
    
//     while (std::getline(inputfile, line)){
//         line
//     }
// }