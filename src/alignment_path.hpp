
#ifndef RPVG_SRC_ALIGNMENTPATH_HPP
#define RPVG_SRC_ALIGNMENTPATH_HPP

#include <vector>

#include "gbwt/gbwt.h"
#include "sparsepp/spp.h"
#include "vg/io/basic_stream.hpp"

#include "paths_index.hpp"
#include "utils.hpp"

using namespace std;


class AlignmentSearchPath;

class AlignmentPath {

    public: 
        
        AlignmentPath(const AlignmentSearchPath & align_search_path);

        uint32_t seq_length;
        uint32_t mapq_comb;
        uint32_t score_sum;
        
        gbwt::SearchState search;

        static vector<AlignmentPath> alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths);
};

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs);
bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs);

ostream & operator<<(ostream & os, const AlignmentPath & align_path);
ostream & operator<<(ostream & os, const vector<AlignmentPath> & align_paths);

namespace std {

    template<> 
    struct hash<vector<AlignmentPath> >
    {
        size_t operator()(vector<AlignmentPath> const & align_paths) const
        {
            size_t seed = 0;

            for (auto & align_path: align_paths) {

                spp::hash_combine(seed, align_path.seq_length);
                spp::hash_combine(seed, align_path.mapq_comb);
                spp::hash_combine(seed, align_path.score_sum);

                spp::hash_combine(seed, align_path.search.node);
                spp::hash_combine(seed, align_path.search.range.first);
                spp::hash_combine(seed, align_path.search.range.second);
            }

            return seed;
        }
    };
}

class AlignmentSearchPath {

    public: 
    
        AlignmentSearchPath();

        vector<gbwt::node_type> path;
        uint32_t path_end_pos;

        uint32_t seq_start_offset;
        uint32_t seq_end_offset;

        gbwt::SearchState search;

        uint32_t seq_length;

        vector<uint32_t> mapqs;
        vector<uint32_t> scores;

        double mapqProb() const;
        uint32_t mapqComb() const;
        uint32_t scoreSum() const;

        bool complete() const;
};

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path);
ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path);


#endif
