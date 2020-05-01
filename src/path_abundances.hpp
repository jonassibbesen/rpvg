
#ifndef RPVG_SRC_PATHABUNDANCES_HPP
#define RPVG_SRC_PATHABUNDANCES_HPP

#include <vector>

#include <Eigen/Dense>

#include "path_info.hpp"
#include "utils.hpp"

using namespace std;


struct Abundances {
        
    Eigen::RowVectorXd confidence;
    Eigen::RowVectorXd expression;

    double read_count;
    	
   	Abundances() {}
    Abundances(const uint32_t num_components, const bool init_zero) {

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

struct PathAbundances : public Abundances {
        
    vector<PathInfo> paths;

    PathAbundances(const vector<PathInfo> & paths_in, const bool add_noise, const bool init_zero) : Abundances(paths_in.size() + static_cast<uint32_t>(add_noise), init_zero), paths(paths_in) {}
};


#endif
