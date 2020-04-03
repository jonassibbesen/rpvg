
#ifndef RPVG_SRC_PATHABUNDANCEESTIMATOR_HPP
#define RPVG_SRC_PATHABUNDANCEESTIMATOR_HPP

#include <vector>

#include <Eigen/Dense>

#include "path_abundances.hpp"
#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


namespace Eigen {

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor> ColVectorXd;
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMatrixXd;
}

class PathAbundanceEstimator {

    public:

        PathAbundanceEstimator(const double min_abundance_in);
        virtual ~PathAbundanceEstimator() {};

        virtual Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const = 0;

    protected: 

        const double min_abundance;

        void removeNoiseAndRenormalize(Abundances * abundances) const;
};

class EMPathAbundanceEstimator : public PathAbundanceEstimator {

    public:

        EMPathAbundanceEstimator(const double min_abundance_in, const uint32_t max_em_iteration_in);
        ~EMPathAbundanceEstimator() {};

        Abundances inferPathClusterAbundances(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const uint32_t num_paths) const;

    private:
    
        const uint32_t max_em_iteration;
};


#endif
