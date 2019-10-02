
#ifndef RPVG_SRC_ALIGNMENTPATH_HPP
#define RPVG_SRC_ALIGNMENTPATH_HPP

#include <vector>

#include "gbwt/gbwt.h"
#include "vg/io/basic_stream.hpp"
#include "utils.hpp"

using namespace std;


class AlignmentPath {

    public: 
    
        AlignmentPath();

        vector<gbwt::node_type> path;
        int32_t path_end_pos;
        int32_t seq_end_offset;

        gbwt::SearchState search;
        vector<gbwt::size_type> ids;

        int32_t seq_length;

        vector<int32_t> mapqs;
        vector<int32_t> scores;

        int32_t mapqMin() const;
        double mapqProb() const;
        int32_t scoreSum() const;
};

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs);

ostream & operator<<(ostream & os, const AlignmentPath & align_path);
ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_paths);

#endif
