
#ifndef FERSKEN_SRC_PATHABUNDANCES_HPP
#define FERSKEN_SRC_PATHABUNDANCES_HPP

#include <vector>

#include <Eigen/Dense>

#include "path.hpp"
#include "utils.hpp"

using namespace std;


struct Abundances {
        
    Eigen::RowVectorXd confidence;
    Eigen::RowVectorXd expression;

    double read_count;
    	
   	Abundances() {}
    Abundances(const uint32_t num_components, const bool init_zero = false) {

        if (init_zero) {

            confidence = Eigen::RowVectorXd::Zero(1, num_components);
            expression = Eigen::RowVectorXd::Zero(1, num_components);

        } else {

            confidence = Eigen::RowVectorXd::Constant(1, num_components, 1);
            expression = Eigen::RowVectorXd::Constant(1, num_components, 1 / static_cast<float>(num_components));
        }

        read_count = 0;
    }
};

struct PathAbundances {
        
    vector<Path> paths;
    Abundances abundances;

    PathAbundances(const vector<Path> & paths_in, const bool add_noise, const bool init_zero = false) : paths(paths_in){

        abundances = Abundances(paths.size() + static_cast<uint32_t>(add_noise), init_zero);
    }
};


#endif
