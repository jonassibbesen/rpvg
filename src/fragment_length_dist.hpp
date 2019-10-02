
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

        double mean() const;
        double sd() const;

        bool isValid() const;
        int32_t maxLength() const;
        double logProb(const int32_t value) const;

        bool parseAlignment(const vg::Alignment & alignment);
        bool parseMultipathAlignment(const vg::MultipathAlignment & alignment);

    private:

        double mean_;
        double sd_;
};


#endif
