
#ifndef RPVG_SRC_PATHCLUSTERESTIMATES_HPP
#define RPVG_SRC_PATHCLUSTERESTIMATES_HPP

#include <vector>

#include <Eigen/Dense>

#include "utils.hpp"

using namespace std;


struct PathInfo {
        
    string name;
    string origin;
    uint32_t count;
    uint32_t length;
    double effective_length;
    
    PathInfo() {

        name = "";
        origin = "";
        count = 1;
        length = 0;
        effective_length = 0;
    }
};

struct PathClusterEstimates {

    vector<PathInfo> paths;

    Eigen::RowVectorXd posteriors;
    Eigen::RowVectorXd abundances;

    double read_count;

    vector<vector<uint32_t> > path_groups;

    PathClusterEstimates() {

        read_count = 0;
    }

    void generateGroupsRecursive(const uint32_t num_components, const uint32_t group_size, vector<uint32_t> cur_group) {

        assert(cur_group.size() <= group_size);

        if (cur_group.size() < group_size) {

            uint32_t start_idx = 0;

            if (!cur_group.empty()) {

                start_idx = cur_group.back();
            }

            for (uint32_t i = start_idx; i < num_components; ++i) {

                vector<uint32_t> new_group = cur_group;
                new_group.push_back(i);
                generateGroupsRecursive(num_components, group_size, new_group);
            }

        } else {

            path_groups.emplace_back(cur_group);
        }
    }

    void initEstimates(uint32_t num_components, const uint32_t group_size, const bool init_zero) {

        if (group_size > 0) {

            generateGroupsRecursive(num_components, group_size, vector<uint>());
            num_components = path_groups.size();
        }

        if (init_zero) {

            posteriors = Eigen::RowVectorXd::Zero(1, num_components);
            abundances = Eigen::RowVectorXd::Zero(1, num_components);

        } else {

            posteriors = Eigen::RowVectorXd::Constant(num_components, 1);
            abundances = Eigen::RowVectorXd::Constant(num_components, 1 / static_cast<float>(num_components));
        }

        read_count = 0;
    }
};


#endif
