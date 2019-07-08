
#ifndef VGPROB_ALIGNMENTPATH_HPP
#define VGPROB_ALIGNMENTPATH_HPP

#include <vector>

#include <gbwt/gbwt.h>
#include <vg/io/basic_stream.hpp>

using namespace std;


class AlignmentPath {

    public: 
    
        AlignmentPath();

        gbwt::SearchState path;
        vector<gbwt::size_type> path_ids;

        int32_t node_length;
        int32_t seq_length;

        vector<int32_t> mapqs;
        vector<int32_t> scores;

        int32_t mapqMin() const;
        double mapqProb() const;
        int32_t scoreSum() const;
};

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs);

ostream& operator<<(ostream& os, const vector<int32_t> & values);
ostream& operator<<(ostream& os, const AlignmentPath & align_path);
ostream& operator<<(ostream& os, const vector<AlignmentPath> & align_paths);

#endif

