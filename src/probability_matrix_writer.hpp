
#ifndef RPVG_SRC_PROBABILITYMATRIXWRITER_HPP
#define RPVG_SRC_PROBABILITYMATRIXWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "read_path_probabilities.hpp"
#include "path.hpp"
#include "utils.hpp"

using namespace std;


class ProbabilityMatrixWriter {

    public: 
    	
    	ProbabilityMatrixWriter(const bool use_stdout_in, const string filename, const double precision_in);
    	~ProbabilityMatrixWriter();

    	void lockWriter();
    	void unlockWriter();

        void writeReadPathProbabilityCluster(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths);

        bool collapseReadPathProbabilities(const ReadPathProbabilities & cluster_probs_1, const ReadPathProbabilities & cluster_probs_2) const;

    private:

    	const bool use_stdout;
        const double precision;
        const uint32_t precision_digits;

	    ofstream writer_file;
    	ostream * writer_stream;

    	mutex writer_mutex;

        void writeCollapsedProbabilities(const vector<pair<double, vector<uint32_t> > > & collpased_probs, const bool write_zero);
};


#endif
