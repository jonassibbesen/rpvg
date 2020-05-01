
#ifndef RPVG_SRC_PATHLIKELIHOODS_HPP
#define RPVG_SRC_PATHLIKELIHOODS_HPP

#include <vector>
#include <limits>

#include <Eigen/Dense>

#include "path_info.hpp"
#include "utils.hpp"

using namespace std;


struct Likelihoods {
        
    Eigen::RowVectorXd likelihoods;
    bool is_log; 
    	
   	Likelihoods() {}
    Likelihoods(const uint32_t num_components, const bool init_log) : is_log(init_log) {

        if (init_log) {

            likelihoods = Eigen::RowVectorXd::Constant(1, num_components, numeric_limits<double>::lowest());

        } else {

            likelihoods = Eigen::RowVectorXd::Zero(1, num_components);
        }
    }
};

struct PathLikelihoods : Likelihoods {
        
    vector<PathInfo> paths;

    PathLikelihoods(const vector<PathInfo> & paths_in, const bool init_log) : Likelihoods(paths_in.size(), init_log), paths(paths_in) {
    }
};

struct PathComboLikelihoods : Likelihoods {
   
    vector<PathInfo> paths;     
    vector<vector<uint32_t> > path_combos;

    PathComboLikelihoods(const vector<PathInfo> & paths_in, const uint32_t combo_size, const bool init_log) : Likelihoods((paths_in.size() * (paths_in.size() - 1) / 2) * (combo_size - 1) + paths_in.size(), init_log), paths(paths_in) {

        assert(combo_size == 2);
        path_combos.reserve((paths.size() * (paths.size() - 1) / 2) * (combo_size - 1) + paths.size());

        if (combo_size == 1) {
            
            for (uint32_t i = 0; i < paths.size(); ++i) {

                path_combos.emplace_back(vector<uint32_t>({i}));
            }

        } else {

            for (uint32_t i = 0; i < paths.size(); ++i) {

                for (uint32_t j = i; j < paths.size(); ++j) {

                    path_combos.emplace_back(vector<uint32_t>({i, j}));
                }
            }
        }
    }
};


#endif
