
#ifndef RPVG_SRC_PATHESTIMATESWRITER_HPP
#define RPVG_SRC_PATHESTIMATESWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "path_cluster_estimates.hpp"
#include "utils.hpp"

using namespace std;


class PathEstimatesWriter {

    public: 
    	
    	PathEstimatesWriter(const bool use_stdout_in, const string filename);
    	~PathEstimatesWriter();

        void writeThreadedPathClusterLikelihoods(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const uint32_t ploidy);
        void writeThreadedPathClusterAbundances(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates);

    private:

    	const bool use_stdout;

	    ofstream writer_file;
    	ostream * writer_stream;
};


#endif
