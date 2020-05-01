
#ifndef RPVG_SRC_PATHABUNDANCES_HPP
#define RPVG_SRC_PATHABUNDANCES_HPP

#include <vector>
#include <limits>

#include <Eigen/Dense>

#include "path_info.hpp"
#include "utils.hpp"

using namespace std;


struct Likelihoods {
        
    Eigen::RowVectorXd likelihoods;
    	
   	Likelihoods() {}
    Likelihoods(const uint32_t num_components, const bool init_log) {

        if (init_log) {

            likelihoods = Eigen::RowVectorXd::Constant(1, num_components, numeric_limits<double>::lowest());

        } else {

            likelihoods = Eigen::RowVectorXd::Zero(1, num_components);
        }
    }
};

struct PathLikelihoods {
        
    vector<PathInfo> paths;
    Likelihoods likelihoods;

    PathAbundances(const vector<PathInfo> & paths_in, const bool add_noise, const bool init_log) : paths(paths_in) {

        abundances = Abundances(paths.size() + static_cast<uint32_t>(add_noise), init_log);
    }
};

struct PloidyPathLikelihoods {
        
    vector<PathInfo> paths;
    vector<vector<uint32_t> > path_combo;

    Likelihoods likelihoods;

    PathAbundances(const vector<PathInfo> & paths_in, const uint32_t ploidy, const bool init_log) : paths(paths_in) {

        assert(ploidy >= 1 && ploidy <= 2);
        path_combo.reserve(paths.size() * (paths.size() - 1) / 2 + paths.size());

        if (ploidy == 1) {
            
            for (uint32_t i = 0; i < paths.size(); ++i) {

                path_combo.emplace_back(vector<uint32_t>({i}));
            }

        } else {

            for (uint32_t i = 0; i < paths.size(); ++i) {

                for (uint32_t j = i; j < paths.size(); ++j) {
        
                    path_combo.emplace_back(vector<uint32_t>({i, j}));
                }
            }
        }

        likelihoods = Likelihoods(path_combo.size(), init_log);
    }
};


#endif
