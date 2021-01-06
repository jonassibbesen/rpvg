
#ifndef RPVG_SRC_ALIGNMENTPATH_HPP
#define RPVG_SRC_ALIGNMENTPATH_HPP

#include <vector>

#include "gbwt/gbwt.h"
#include "sparsepp/spp.h"
#include "vg/io/basic_stream.hpp"

#include "utils.hpp"

using namespace std;


class AlignmentSearchPath;

class AlignmentPath {

    public: 
        
        AlignmentPath(const uint32_t frag_length_in, const uint32_t min_mapq_in, const uint32_t score_sum_in, const bool is_multimap_in, const gbwt::SearchState & search_state_in);
        AlignmentPath(const AlignmentSearchPath & align_path_in, const bool is_multimap_in);

        uint32_t frag_length;
        uint32_t min_mapq;
        uint32_t score_sum;

        bool is_multimap;
        gbwt::SearchState search_state;

        static vector<AlignmentPath> alignmentSearchPathsToAlignmentPaths(const vector<AlignmentSearchPath> & align_search_paths, const uint32_t max_score_diff, const bool is_multimap);
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

                spp::hash_combine(seed, align_path.frag_length);
                spp::hash_combine(seed, align_path.min_mapq);
                spp::hash_combine(seed, align_path.score_sum);
                spp::hash_combine(seed, align_path.is_multimap);
                spp::hash_combine(seed, align_path.search_state.node);
                spp::hash_combine(seed, align_path.search_state.range.first);
                spp::hash_combine(seed, align_path.search_state.range.second);
            }

            return seed;
        }
    };
}

struct ReadAlignmentStats {

    uint32_t mapq;
    int32_t score;
    uint32_t length;

    int32_t left_softclip_length;
    int32_t right_softclip_length;

    ReadAlignmentStats() {

        mapq = 0;
        score = 0;
        length = 0;

        left_softclip_length = -1;
        right_softclip_length = -1; 
    }
};

class AlignmentSearchPath {

    public: 
    
        AlignmentSearchPath();

        vector<gbwt::node_type> path;
        uint32_t path_end_idx;

        uint32_t path_start_offset;        
        uint32_t path_end_offset;

        gbwt::SearchState search_state;

        int32_t insert_length;

        vector<ReadAlignmentStats> read_stats;

        uint32_t fragmentLength() const;

        uint32_t minMappingQuality() const;
        uint32_t scoreSum() const;

        double minBestScoreFraction() const;
        double maxSoftclipFraction() const;

        bool isComplete() const;
};

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path);
ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path);


#endif
