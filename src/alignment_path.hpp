
#ifndef RPVG_SRC_ALIGNMENTPATH_HPP
#define RPVG_SRC_ALIGNMENTPATH_HPP

#include <vector>

#include "gbwt/gbwt.h"
#include "vg/io/basic_stream.hpp"
#include "utils.hpp"

using namespace std;

class AlignmentSearchPath;

class AlignmentPath {

    public: 
        
        AlignmentPath(const int32_t seq_length_in, const double mapq_prob_in, const int32_t score_sum_in, const vector<gbwt::size_type> & ids_in);
        AlignmentPath(const AlignmentSearchPath & align_path_in, const vector<gbwt::size_type> & ids_in);

        const int32_t seq_length;
        const double mapq_prob;
        const int32_t score_sum;
        
        const vector<gbwt::size_type> ids;
};

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs);

ostream & operator<<(ostream & os, const AlignmentPath & align_path);
ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_paths);

class AlignmentSearchPath {

    public: 
    
        AlignmentSearchPath();

        vector<gbwt::node_type> path;
        int32_t path_end_pos;
        int32_t seq_end_offset;

        gbwt::SearchState search;

        int32_t seq_length;

        vector<int32_t> mapqs;
        vector<int32_t> scores;

        double mapqProb() const;
        int32_t scoreSum() const;

        bool complete() const;
};

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path);
ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path);

#endif
