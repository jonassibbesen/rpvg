
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

    vector<double> noise_samples;
    vector<double> abundance_samples;

    CountSamples() {}
};

struct PathClusterEstimates {

    vector<PathInfo> paths;

    vector<vector<uint32_t> > path_group_sets;

    vector<double> posteriors;
    vector<double> abundances;

    double noise_count;
    double total_count;

    vector<CountSamples> gibbs_read_count_samples;

    PathClusterEstimates() {

        noise_count = 0;
        total_count = 0;
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

            path_group_sets.emplace_back(cur_group);
        }
    }

    void resetEstimates(uint32_t num_components, const uint32_t group_size) {

        path_group_sets.clear();

        posteriors.clear();
        abundances.clear();

        noise_count = 0;
        total_count = 0; 
        
        gibbs_read_count_samples.clear();

        if (group_size > 0) {

            generateGroupsRecursive(num_components, group_size, vector<uint>());
        
            posteriors = vector<double>(path_group_sets.size(), 0);
            abundances = vector<double>(path_group_sets.size() * group_size, 0);
        }
    }
};


#endif
