
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
        
        AlignmentPath(const pair<gbwt::SearchState, gbwt::size_type> & gbwt_search_in, const bool is_multimap_in, const uint32_t frag_length_in, const uint32_t min_mapq_in, const uint32_t score_sum_in);
        AlignmentPath(const AlignmentSearchPath & align_path_in, const bool is_multimap_in);

        pair<gbwt::SearchState, gbwt::size_type> gbwt_search;
        bool is_multimap;

        uint32_t frag_length;
        uint32_t min_mapq;
        uint32_t score_sum;

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

                spp::hash_combine(seed, align_path.gbwt_search.first.node);
                spp::hash_combine(seed, align_path.gbwt_search.first.range.first);
                spp::hash_combine(seed, align_path.gbwt_search.first.range.second);
                spp::hash_combine(seed, align_path.gbwt_search.second);
                spp::hash_combine(seed, align_path.is_multimap);
                spp::hash_combine(seed, align_path.frag_length);
                spp::hash_combine(seed, align_path.min_mapq);
                spp::hash_combine(seed, align_path.score_sum);
            }

            return seed;
        }
    };
}


class InternalAlignment {

    public: 

        InternalAlignment();

        bool is_internal;
        uint32_t penalty;

        uint32_t offset;
        uint32_t max_offset;
};

bool operator==(const InternalAlignment & lhs, const InternalAlignment & rhs);
bool operator!=(const InternalAlignment & lhs, const InternalAlignment & rhs);
bool operator<(const InternalAlignment & lhs, const InternalAlignment & rhs);

ostream & operator<<(ostream & os, const InternalAlignment & read_align_stats);


class AlignmentStats {

    public: 

        AlignmentStats();

        uint32_t mapq;
        int32_t score;
        uint32_t length;

        uint32_t left_softclip_length;
        uint32_t right_softclip_length;

        InternalAlignment internal_start;
        InternalAlignment internal_end;

        gbwt::node_type internal_end_next_node;

        void updateLeftSoftClipLength(const vg::Path & path);
        void updateRightSoftClipLength(const vg::Path & path);

        int32_t adjustedScore() const; 

        uint32_t clippedOffsetLeftBases() const;
        uint32_t clippedOffsetRightBases() const;
        uint32_t clippedOffsetTotalBases() const;
};

bool operator==(const AlignmentStats & lhs, const AlignmentStats & rhs);
bool operator!=(const AlignmentStats & lhs, const AlignmentStats & rhs);
bool operator<(const AlignmentStats & lhs, const AlignmentStats & rhs);

ostream & operator<<(ostream & os, const AlignmentStats & read_align_stats);


class AlignmentSearchPath {

    public: 
    
        AlignmentSearchPath();

        vector<gbwt::node_type> path;
        pair<gbwt::SearchState, gbwt::size_type> gbwt_search;

        uint32_t start_offset;
        uint32_t end_offset;

        int32_t insert_length;

        vector<AlignmentStats> read_align_stats;

        uint32_t fragmentLength() const;
        uint32_t minMappingQuality() const;
        uint32_t scoreSum() const;

        double minOptimalScoreFraction(const vector<int32_t> & optimal_align_scores) const;
        double maxSoftclipFraction() const;

        bool isEmpty() const;
        void clear();
};

bool operator==(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs);
bool operator!=(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs);
bool operator<(const AlignmentSearchPath & lhs, const AlignmentSearchPath & rhs);

ostream & operator<<(ostream & os, const AlignmentSearchPath & align_search_path);
ostream & operator<<(ostream & os, const vector<AlignmentSearchPath> & align_search_path);


#endif
