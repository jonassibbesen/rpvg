
#ifndef RPVG_SRC_PATHCLUSTERESTIMATES_HPP
#define RPVG_SRC_PATHCLUSTERESTIMATES_HPP

#include <vector>

#include "Eigen/Dense"
#include "sparsepp/spp.h"

#include "utils.hpp"

using namespace std;


struct PathInfo {
        
    string name;
    uint32_t group_id;

    uint32_t source_count;
    spp::sparse_hash_set<uint32_t> source_ids;

    uint32_t length;
    double effective_length;
    
    PathInfo(const string & name_in) : name(name_in) {

        group_id = 0;
        source_count = 1;
        length = 0;
        effective_length = 0;
    }
};

struct CountSamples {
    
    vector<uint32_t> path_ids;
    vector<vector<double> > samples;
};

struct PathClusterEstimates {

    vector<PathInfo> paths;

    vector<double> posteriors;
    vector<vector<uint32_t> > path_group_sets;

    Eigen::RowVectorXd abundances;

    vector<CountSamples> gibbs_read_count_samples;

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

            path_group_sets.emplace_back(cur_group);
        }
    }

    void initEstimates(uint32_t num_components, const uint32_t group_size, const bool init_zero) {

        if (group_size > 0) {

            generateGroupsRecursive(num_components, group_size, vector<uint>());
            num_components = path_group_sets.size();
        
            if (init_zero) {

                posteriors = vector<double>(num_components, 0);

            } else {

                posteriors = vector<double>(num_components, 1 / static_cast<float>(num_components));
            }
        }

        if (init_zero) {

            abundances = Eigen::RowVectorXd::Zero(1, num_components);

        } else {

            abundances = Eigen::RowVectorXd::Constant(1, num_components, 1 / static_cast<float>(num_components));
        }
    }
};


#endif
