
#ifndef RPVG_SRC_ABUNDANCEESTIMATOR_HPP
#define RPVG_SRC_ABUNDANCEESTIMATOR_HPP

#include <vector>

#include <Eigen/Dense>

#include "abundances.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


namespace Eigen {

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXd;
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXd;
}

class AbundanceEstimator {

    public:

        AbundanceEstimator(const double min_abundance_in);
        virtual ~AbundanceEstimator() {};

        virtual Abundances inferClusterAbundance(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) = 0;

    protected: 

        const double min_abundance;
};

class EMAbundanceEstimator : public AbundanceEstimator {

    public:

        EMAbundanceEstimator(const double min_abundance_in, const uint32_t max_em_iteration_in, const double stop_em_count_diff_in);
        ~EMAbundanceEstimator() {};

        Abundances inferClusterAbundance(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths);

    private:
    
        const uint32_t max_em_iteration;
        const double stop_em_count_diff;
};


#endif
