
#ifndef FERSKEN_SRC_PATHABUNDANCES_HPP
#define FERSKEN_SRC_PATHABUNDANCES_HPP

#include <vector>

#include <Eigen/Dense>

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
        
    vector<string> names;
    vector<uint32_t> lengths;
    vector<double> effective_lengths;
       
    Abundances abundances;

    PathAbundances(const uint32_t num_components) {

        names.reserve(num_components);
        lengths.reserve(num_components);
        effective_lengths.reserve(num_components);
    }
};


#endif
