
#ifndef VGPROB_ALIGNMENTPATH_HPP
#define VGPROB_ALIGNMENTPATH_HPP

#include <vector>

#include <gbwt/gbwt.h>
#include <vg/io/basic_stream.hpp>

#include "utils.hpp"

using namespace std;


class AlignmentPath {

    public: 
    
        AlignmentPath();

        gbwt::SearchState path; 
        vector<gbwt::size_type> path_ids;
        
        int32_t node_length;
        int32_t seq_length;

        pair<int32_t, int32_t> scores;
        pair<int32_t, int32_t> mapqs;

        void extent_align_path(const vg::Path & path, const uint32_t & node_offset, const gbwt::GBWT & paths_index);
};

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs);
ostream& operator<<(ostream& os, const vector<AlignmentPath> & align_paths);

#endif

