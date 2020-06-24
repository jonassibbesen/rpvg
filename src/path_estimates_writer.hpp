
#ifndef RPVG_SRC_PATHESTIMATESWRITER_HPP
#define RPVG_SRC_PATHESTIMATESWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "gbwt/gbwt.h"

#include "path_cluster_estimates.hpp"
#include "paths_index.hpp"
#include "utils.hpp"

using namespace std;


class PathEstimatesWriter {

    public: 
    	
    	PathEstimatesWriter(const bool use_stdout_in, const string filename);
    	~PathEstimatesWriter();

        void writePathClusterPosteriors(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const uint32_t ploidy);
        void writePathClusterAbundances(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates);

        void writeBestPathGroupsPosteriors(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const PathsIndex & paths_index, const bool is_call_traversals);
        void writeMarginalPathPosteriors(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const PathsIndex & paths_index, const bool is_call_traversals);

    private:

        pair<string, uint32_t> gbwtPathToGAFPath(const PathInfo & path, const PathsIndex & paths_index, const bool is_call_traversals);
        void writeGAFline(const string & path_name, const pair<string, uint32_t> & gaf_path, const double path_posterior, const uint32_t read_count);

    	const bool use_stdout;

	    ofstream writer_file;
    	ostream * writer_stream;
};


#endif
