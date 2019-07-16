
#ifndef VGPROB_FRAGMENTLENGTHDIST_HPP
#define VGPROB_FRAGMENTLENGTHDIST_HPP

#include <iostream>
#include <fstream>

using namespace std;


class FragmentLengthDist {

    public: 
    	
        FragmentLengthDist();
        FragmentLengthDist(const double mean_in, const double sd_in);
        FragmentLengthDist(ifstream * alignments_istream, const bool is_multipath);

        double mean() const;
        double sd() const;

        bool isValid() const;
        int32_t maxLength() const;
        double logProb(const int32_t value) const;

    private:

        double mean_;
        double sd_;
};

#endif

