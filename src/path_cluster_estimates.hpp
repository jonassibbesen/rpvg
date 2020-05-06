
#ifndef RPVG_SRC_PATHCLUSTERESTIMATES_HPP
#define RPVG_SRC_PATHCLUSTERESTIMATES_HPP

#include <vector>

#include <Eigen/Dense>

#include "utils.hpp"

using namespace std;


struct PathInfo {
        
    string name;
    string origin;
    uint32_t length;
    double effective_length;
    
    PathInfo() {

        name = "";
        origin = "";
        length = 0;
        effective_length = 0;
    }
};

struct Likelihoods {
        
    Eigen::RowVectorXd likelihoods;
    bool is_log; 
        
    vector<vector<uint32_t> > groups;

    Likelihoods() {}
    Likelihoods(const uint32_t num_components, const uint32_t group_size, const bool init_log) : is_log(init_log) {

        assert(group_size <= 2);
        groups.reserve((num_components * (num_components - 1) / 2) * (group_size - 1) + num_components);

        if (group_size == 1) {
            
            for (uint32_t i = 0; i < num_components; ++i) {

                groups.emplace_back(vector<uint32_t>({i}));
            }

        } else {

            for (uint32_t i = 0; i < num_components; ++i) {

                for (uint32_t j = i; j < num_components; ++j) {

                    groups.emplace_back(vector<uint32_t>({i, j}));
                }
            }
        }

        if (init_log) {

            likelihoods = Eigen::RowVectorXd::Constant(1, groups.size(), numeric_limits<double>::lowest());

        } else {

            likelihoods = Eigen::RowVectorXd::Zero(1, groups.size());
        }
    }
};

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

struct PathClusterEstimates {

    vector<PathInfo> paths;
    Abundances abundances;
    Likelihoods likelihoods;
};


#endif
