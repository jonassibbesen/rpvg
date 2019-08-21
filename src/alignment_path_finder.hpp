
#ifndef VGPROB_SRC_ALIGNMENTPATHFINDER_HPP
#define VGPROB_SRC_ALIGNMENTPATHFINDER_HPP

#include <vector>
#include <map>

#include "gbwt/gbwt.h"
#include "vg/io/basic_stream.hpp"
#include "alignment_path.hpp"

using namespace std;


template<class AlignmentType> 
class AlignmentPathFinder {

    public: 
    
       	AlignmentPathFinder(const vg::Graph & graph, const gbwt::GBWT & paths_index_in, const int32_t max_pair_seq_length_in);
       	void setMaxPairSeqLength(const int32_t max_pair_seq_length_in);

		vector<AlignmentPath> findAlignmentPaths(const AlignmentType & alignment) const;
		vector<AlignmentPath> findAlignmentPathsIds(const AlignmentType & alignment) const;

		vector<AlignmentPath> findPairedAlignmentPaths(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;
		vector<AlignmentPath> findPairedAlignmentPathsIds(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;

	private:

        vector<int32_t> node_seq_lengths;
       	const gbwt::GBWT & paths_index;

       	int32_t max_pair_seq_length;

		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const vector<gbwt::node_type> & stop_nodes = vector<gbwt::node_type>()) const;
		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const int32_t sub_path_start_idx, const vector<gbwt::node_type> & stop_nodes = vector<gbwt::node_type>()) const;
		void extendAlignmentPath(AlignmentPath * align_path, const vg::Path & path, const vector<gbwt::node_type> & stop_nodes = vector<gbwt::node_type>()) const;

		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const vector<gbwt::node_type> & stop_nodes = vector<gbwt::node_type>()) const;
		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const int32_t sub_path_start_idx, const vector<gbwt::node_type> & stop_nodes = vector<gbwt::node_type>()) const;
		void extendAlignmentPaths(vector<AlignmentPath> * align_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & subpaths, const int32_t sub_path_start_idx, const vector<gbwt::node_type> & stop_nodes = vector<gbwt::node_type>()) const;

		void pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const AlignmentType & end_alignment) const;

		vector<gbwt::node_type> getAlignmentStartNodes(const vg::Alignment & alignment) const;
		vector<gbwt::node_type> getAlignmentStartNodes(const vg::MultipathAlignment & alignment) const;

		map<gbwt::node_type, int32_t> getAlignmentStartNodesIndex(const vg::Alignment & alignment) const;
		map<gbwt::node_type, int32_t> getAlignmentStartNodesIndex(const vg::MultipathAlignment & alignment) const;

		void printDebug(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;
};


#endif
