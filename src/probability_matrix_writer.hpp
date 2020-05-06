
#ifndef RPVG_SRC_PROBABILITYMATRIXWRITER_HPP
#define RPVG_SRC_PROBABILITYMATRIXWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "read_path_probabilities.hpp"
#include "path_cluster_estimates.hpp"
#include "utils.hpp"

using namespace std;


class ProbabilityMatrixWriter {

    public: 
    	
    	ProbabilityMatrixWriter(const bool use_stdout_in, const string filename, const double prob_precision_in);
    	~ProbabilityMatrixWriter();

    	void lockWriter();
    	void unlockWriter();

        void writeReadPathProbabilityCluster(const vector<ReadPathProbabilities> & cluster_probs, const vector<PathInfo> & cluster_paths);

    private:

    	const bool use_stdout;
        const double prob_precision;
        const uint32_t prob_precision_digits;

	    ofstream writer_file;
    	ostream * writer_stream;

    	mutex writer_mutex;

        void writeCollapsedProbabilities(const vector<pair<double, vector<uint32_t> > > & collpased_probs, const bool write_zero);
};


#endif
