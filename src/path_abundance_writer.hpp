
#ifndef RPVG_SRC_PATHABUNDANCEWRITER_HPP
#define RPVG_SRC_PATHABUNDANCEWRITER_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>

#include "path_abundances.hpp"
#include "utils.hpp"

using namespace std;


class PathAbundanceWriter {

    public: 
    	
    	PathAbundanceWriter(const bool use_stdout_in, const string filename, const double precision);
    	~PathAbundanceWriter();

        void writeThreadedPathClusterAbundances(const vector<vector<PathAbundances> > & threaded_path_cluster_abundances);

    private:

    	const bool use_stdout;
        const uint32_t precision_digits;

	    ofstream writer_file;
    	ostream * writer_stream;
};


#endif
