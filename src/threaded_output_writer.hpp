
#ifndef RPVG_SRC_THREADEDOUTPUTWRITER_HPP
#define RPVG_SRC_THREADEDOUTPUTWRITER_HPP

#include <iostream>
#include <fstream>
#include <string>

#include "htslib/bgzf.h"
#include "htslib/hts.h"

#include "producer_consumer_queue.hpp"
#include "read_path_probabilities.hpp"
#include "path_cluster_estimates.hpp"
#include "utils.hpp"

using namespace std;

class ThreadedOutputWriter {

    public: 

        ThreadedOutputWriter(const string & filename, const string & compression_mode, const uint32_t num_threads);
        virtual ~ThreadedOutputWriter() {};

        void close();

    protected:

        ProducerConsumerQueue<stringstream *> * output_queue;

    private:

        BGZF * writer_stream;
        thread writing_thread; 

        void write();
};

class ProbabilityClusterWriter : public ThreadedOutputWriter {

    public: 
        
        ProbabilityClusterWriter(const string filename_prefix, const uint32_t num_threads, const double prob_precision_in);
        ~ProbabilityClusterWriter() {};

        void addCluster(const vector<ReadPathProbabilities> & cluster_probs, const vector<PathInfo> & cluster_paths);

    private:

        const double prob_precision;
        const uint32_t prob_precision_digits;

        void addCollapsedProbabilities(stringstream * prob_out_sstream, const ReadPathProbabilities & probs);
};

class GibbsSamplesWriter : public ThreadedOutputWriter {

    public: 
        
        GibbsSamplesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t num_gibbs_samples_in);
        ~GibbsSamplesWriter() {};

        void addSamples(const pair<uint32_t, PathClusterEstimates> & path_cluster_estimate);

    private:

        const uint32_t num_gibbs_samples; 
};

class PosteriorEstimatesWriter : public ThreadedOutputWriter {

    public: 
        
        PosteriorEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double min_posterior_in);
        ~PosteriorEstimatesWriter() {};

        void addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates);

    private:

        const uint32_t ploidy;
        const double min_posterior;
};

class AbundanceEstimatesWriter : public ThreadedOutputWriter {

    public: 
    	
    	AbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const double total_transcript_count_in);
    	~AbundanceEstimatesWriter() {};

        void addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates);

    private:

        const double total_transcript_count;
};


#endif
