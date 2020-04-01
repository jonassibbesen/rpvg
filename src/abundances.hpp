
#ifndef RPVG_SRC_ABUNDANCES_HPP
#define RPVG_SRC_ABUNDANCES_HPP

#include <vector>

#include <Eigen/Dense>

#include "utils.hpp"

using namespace std;


struct Abundances {
        
    Eigen::RowVectorXd confidence;
    Eigen::RowVectorXd expression;

    Abundances(const uint32_t num_components) {

        confidence = Eigen::RowVectorXd::Constant(1, num_components + 1, 1);
        expression = Eigen::RowVectorXd::Constant(1, num_components + 1, 1 / static_cast<float>(num_components + 1));
    }
};


#endif
