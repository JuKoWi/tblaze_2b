#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <tb/model.h>
#include <tb/read_and_write.h>
#include <tb/utils.h>
#include <string>
#include <cmath>
#include <sstream>

void write_bands_path(std::vector<Eigen::Vector3d> kfrac, std::vector<Eigen::VectorXd> eigvals){
    if (kfrac.size() != eigvals.size()){
        throw std::invalid_argument("Only one set of eigenvalues per k-point allowed");
    }
    std::vector<Eigen::VectorXd> eigvals_eV {};
    for (auto v: eigvals){
        eigvals_eV.push_back(v * Constants::HARTREE_TO_EV);
    }
    std::ofstream resultfile("bands.txt");
    if (resultfile.is_open()){
        for (int i=0; i<kfrac.size(); i++){
            resultfile << kfrac[i].transpose() << '\t' << eigvals_eV[i].transpose() << '\n';
        } 
        resultfile.close();
    }
    std::cout << "Finished writing bandfile" << '\n';
}

std::istringstream read_line(std::ifstream& infile){
    std::string line {};
    std::getline(infile, line);
    return std::istringstream(line);
}

// parse orthogonal Wannier90 type seedname_tb.dat file
TightBindingModel from_seedname(std::string fname){
    std::ifstream infile(fname);
    if (!infile.is_open()){
        throw std::invalid_argument("Loading seedname file failed");
    }
    std::string header {};
    std::getline(infile, header);
    Eigen::Matrix3d lattice {};
    std::string line {};
    double x, y, z;
    for (int i=0; i<3; i++){
        auto iss {read_line(infile)};
        iss >> x >> y >> z;
        lattice.row(i) << x, y, z;
    }
    lattice *= Constants::ANGSTROM_TO_BOHR;
    int norbs {};
    infile >> norbs;
    int n_latpoints {};
    infile >> n_latpoints;
    std::vector<int> degeneracies(n_latpoints);
    std::cout << "Read degeneracies: " << '\n';
    for (int i = 0; i < n_latpoints; ++i){
        infile >> degeneracies[i];
        std::cout << degeneracies[i];
        if ((i+1)%15 == 0){std::cout << '\n';}
    }
    std::cout << '\n';
    infile >> std::ws;
    // read_line(infile);
    std::vector<Eigen::Vector3i> lat_points {};
    std::vector<Hopping> hoppingsH {};
    std::vector<HoppingVec> hoppings_r {};
    for (int i=0; i < n_latpoints; ++i){
        Eigen::Vector3i R {}; 
        auto iss {read_line(infile)};
        iss >> R(0) >> R(1) >> R(2);
        for (int j=0; j < std::pow(norbs,2); ++j){
            int m,n;
            double ReH, ImH;
            auto iss {read_line(infile)};
            iss >> m >> n >> ReH >> ImH;
            if (m - 1 < 0 || m - 1 >= norbs || n - 1 < 0 || n - 1 >= norbs) {
                throw std::runtime_error("Error parsing orbital indices");
            }
            hoppingsH.push_back({m-1, n-1, R, cplx(ReH/Constants::HARTREE_TO_EV, ImH/Constants::HARTREE_TO_EV)/double(degeneracies[i])});
        }
        infile >> std::ws;
    }
    for (int i=0; i < n_latpoints; ++i){
        Eigen::Vector3i R {}; 
        auto iss {read_line(infile)};
        iss >> R(0) >> R(1) >> R(2);
        for (int j=0; j < std::pow(norbs,2); ++j){
            int m,n;
            double Rerx, Imrx, Rery, Imry, Rerz, Imrz;
            auto iss {read_line(infile)};
            iss >> m >> n >> Rerx >> Imrx >> Rery >> Imry >> Rerz >> Imrz;
            Eigen::Vector3cd r(cplx(Rerx,Imrx), cplx(Rery, Imry), cplx(Rerz, Imrz));
            hoppings_r.push_back({m-1, n-1, R, r * Constants::ANGSTROM_TO_BOHR/double(degeneracies[i])});
        }
        infile >> std::ws;
    }
    TightBindingModel model(lattice, hoppingsH);
    model.set_r_hoppings(hoppings_r);
    std::cout << "Finished reading seedname" << std::endl;
    return model;
}
