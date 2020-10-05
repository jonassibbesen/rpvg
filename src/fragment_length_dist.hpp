
#ifndef RPVG_SRC_FRAGMENTLENGTHDIST_HPP
#define RPVG_SRC_FRAGMENTLENGTHDIST_HPP

#include <iostream>
#include <fstream>

#include "vg/io/basic_stream.hpp"

using namespace std;


class FragmentLengthDist {

    public: 
    	
        FragmentLengthDist();
        FragmentLengthDist(const double mean_in, const double sd_in);
        FragmentLengthDist(istream * alignments_istream, const bool is_multipath);
        FragmentLengthDist(const vector<uint32_t> & frag_length_counts);

        double mean() const;
        double sd() const;

        bool isValid() const;
        uint32_t maxLength() const;
        double logProb(const uint32_t value) const;

        bool parseAlignment(const vg::Alignment & alignment);
        bool parseMultipathAlignment(const vg::MultipathAlignment & alignment);

    private:
        
        double mean_;
        double sd_;
        double max_length_;

        vector<double> log_prob_buffer;

        void setMaxLength();
        void setLogProbBuffer(const uint32_t size); 
};


#endif
