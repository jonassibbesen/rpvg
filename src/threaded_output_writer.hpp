
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

        void addCluster(const vector<ReadPathProbabilities> & read_path_cluster_probs, const vector<PathInfo> & cluster_paths);

    private:

        const double prob_precision;
        const uint32_t prob_precision_digits;
};

class ReadCountGibbsSamplesWriter : public ThreadedOutputWriter {

    public: 
        
        ReadCountGibbsSamplesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t num_gibbs_samples_in);
        ~ReadCountGibbsSamplesWriter() {};

        void addSamples(const pair<uint32_t, PathClusterEstimates> & path_cluster_estimate);
        void addNoiseTranscript(const uint32_t unaligned_read_count);

    private:

        const uint32_t num_gibbs_samples;

        vector<double> noise_counts;
};

class JointHaplotypeEstimatesWriter : public ThreadedOutputWriter {

    public: 
        
        JointHaplotypeEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double min_posterior_in);
        ~JointHaplotypeEstimatesWriter() {};

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
        void addNoiseTranscript(const uint32_t unaligned_read_count);

    private:

        const double total_transcript_count;

        double noise_count;
};

class HaplotypeAbundanceEstimatesWriter : public ThreadedOutputWriter {

    public: 
    	
    	HaplotypeAbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double total_transcript_count_in);
    	~HaplotypeAbundanceEstimatesWriter() {};

        void addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates);
        void addNoiseTranscript(const uint32_t unaligned_read_count);

    private:

        const uint32_t ploidy;
        const double total_transcript_count;

        double noise_count;        
};

class JointHaplotypeAbundanceEstimatesWriter : public ThreadedOutputWriter {

    public: 
        
        JointHaplotypeAbundanceEstimatesWriter(const string filename_prefix, const uint32_t num_threads, const uint32_t ploidy_in, const double min_posterior_in, const double total_transcript_count_in);
        ~JointHaplotypeAbundanceEstimatesWriter() {};

        void addEstimates(const vector<pair<uint32_t, PathClusterEstimates> > & path_cluster_estimates);
        void addNoiseTranscript(const uint32_t unaligned_read_count);

    private:

        const uint32_t ploidy;
        const double min_posterior;
        const double total_transcript_count;

        vector<double> noise_counts;
};

#endif
