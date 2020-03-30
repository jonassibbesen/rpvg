
#ifndef RPVG_SRC_PROBABILITYMATRIXWRITER_HPP
#define RPVG_SRC_PROBABILITYMATRIXWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "read_path_probabilities.hpp"
#include "utils.hpp"

using namespace std;


class ProbabilityMatrixWriter {

    public: 
    	
    	ProbabilityMatrixWriter(const bool use_stdout_in, const string filename);
    	~ProbabilityMatrixWriter();

    	void lockWriter();
    	void unlockWriter();

    	void writeReadPathProbabilities(const vector<pair<ReadPathProbabilities, uint32_t> > & read_path_probs, const vector<string> & path_names);

    private:

    	const bool use_stdout;

	    ofstream writer_file;
    	ostream * writer_stream;

    	mutex writer_mutex;
};


#endif
