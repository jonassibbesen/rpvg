
#ifndef FERSKEN_SRC_PROBABILITYMATRIXWRITER_HPP
#define FERSKEN_SRC_PROBABILITYMATRIXWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class ProbabilityMatrixWriter {

    public: 
    	
    	ProbabilityMatrixWriter(const bool use_stdout_in, const string filename, const double precision_in);
    	~ProbabilityMatrixWriter();

    	void lockWriter();
    	void unlockWriter();

        void writeReadPathProbabilityCluster(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<string> & path_names, const vector<uint32_t> & path_lengths, const vector<double> & effective_path_lengths);

    private:

    	const bool use_stdout;
        const double precision;
        const uint32_t precision_digits;

	    ofstream writer_file;
    	ostream * writer_stream;

    	mutex writer_mutex;

        bool collapseReadPathProbabilities(const ReadPathProbabilities & cluster_probs_1, const ReadPathProbabilities & cluster_probs_2) const;
        void writeCollapsedProbabilities(const vector<pair<double, vector<uint32_t> > > & collpased_probs, const bool write_zero);
};


#endif
