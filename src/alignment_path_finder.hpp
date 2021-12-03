
#ifndef RPVG_SRC_ALIGNMENTPATHFINDER_HPP
#define RPVG_SRC_ALIGNMENTPATHFINDER_HPP

#include <vector>
#include <map>

#include "vg/io/basic_stream.hpp"
#include "paths_index.hpp"
#include "alignment_path.hpp"

using namespace std;


template<class AlignmentType> 
class AlignmentPathFinder {

    public: 
    
       	AlignmentPathFinder(const PathsIndex & paths_index_in, const string library_type_in, const bool score_not_qual_in, const bool use_allelic_mapq_in, const uint32_t max_pair_frag_length_in, const uint32_t max_partial_offset_in, const bool est_missing_noise_prob_in, const int32_t max_score_diff_in, const double min_best_score_filter_in);

		vector<AlignmentPath> findAlignmentPaths(const AlignmentType & alignment) const;
		vector<AlignmentPath> findPairedAlignmentPaths(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;

	private:

       	const PathsIndex & paths_index;
       	const string library_type;

       	const bool score_not_qual;
       	const bool use_allelic_mapq;

       	const uint32_t max_pair_frag_length;
       	const uint32_t max_partial_offset;
       	const bool est_missing_noise_prob;

       	const int32_t max_score_diff;
       	const double min_best_score_filter;

		bool alignmentHasPath(const vg::Alignment & alignment) const;
		bool alignmentHasPath(const vg::MultipathAlignment & alignment) const;
		
       	bool alignmentStartInGraph(const AlignmentType & alignment) const;

       	int32_t alignmentScore(const char & quality) const;
		int32_t alignmentScore(const string & quality, const uint32_t & start_offset, const uint32_t & length) const;

       	int32_t optimalAlignmentScore(const string & quality, const uint32_t seq_length) const;
		int32_t optimalAlignmentScore(const vg::Alignment & alignment) const;
		int32_t optimalAlignmentScore(const vg::MultipathAlignment & alignment) const;

		uint32_t mappingQuality(const AlignmentType & alignment) const;

		vector<AlignmentSearchPath> extendAlignmentSearchPath(const AlignmentSearchPath & align_search_path, const vg::Alignment & alignment) const;

		void extendAlignmentSearchPath(vector<AlignmentSearchPath> * align_search_paths, const vg::Path & path, const bool is_first_path, const bool is_last_path, const string & quality, const uint32_t seq_length, const bool add_internal_start) const;
		void extendAlignmentSearchPath(AlignmentSearchPath * align_search_path, const vg::Mapping & mapping) const;

		vector<AlignmentSearchPath> extendAlignmentSearchPath(const AlignmentSearchPath & align_search_path, const vg::MultipathAlignment & alignment) const;
		void extendAlignmentSearchPaths(vector<AlignmentSearchPath> * align_search_paths, const AlignmentSearchPath & init_align_search_path, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const uint32_t start_subpath_idx, const string & quality, const uint32_t seq_length, spp::sparse_hash_map<pair<uint32_t, uint32_t>, int32_t> * internal_node_subpaths, int32_t * best_align_score, const bool has_right_bonus) const;
		
		void findAlignmentSearchPaths(vector<AlignmentSearchPath> * align_search_paths, const AlignmentType & alignment) const;
		void findPairedAlignmentSearchPaths(vector<AlignmentSearchPath> * paired_align_search_paths, const AlignmentType & start_alignment, const AlignmentType & end_alignment) const;

		void mergeAlignmentSearchPath(AlignmentSearchPath * main_align_search_path, uint32_t main_path_start_idx, const AlignmentSearchPath & second_align_search_path) const;

		vector<gbwt::node_type> getAlignmentStartNodes(const vg::Alignment & alignment) const;
		vector<gbwt::node_type> getAlignmentStartNodes(const vg::MultipathAlignment & alignment) const;

		vector<uint32_t> getAlignmentStartSoftclipLengths(const vg::MultipathAlignment & alignment) const;
		vector<uint32_t> getAlignmentEndSoftclipLengths(const vg::MultipathAlignment & alignment) const;

		bool isAlignmentDisconnected(const vg::Alignment & alignment) const;
		bool isAlignmentDisconnected(const vg::MultipathAlignment & alignment) const;

		bool filterAlignmentSearchPaths(const vector<AlignmentSearchPath> & align_search_paths, const vector<int32_t> & optimal_align_scores) const;
};

namespace std {

    template<> 
    struct hash<pair<uint32_t, uint32_t> >
    {
        size_t operator()(const pair<uint32_t, uint32_t> & values) const
        {
            size_t seed = 0;
            
            spp::hash_combine(seed, values.first);
            spp::hash_combine(seed, values.second);

            return seed;
        }
    };
}


#endif
