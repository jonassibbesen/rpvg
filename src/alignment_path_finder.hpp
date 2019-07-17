
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
    
       	AlignmentPathFinder(const vg::Graph & graph, const gbwt::GBWT & paths_index_in, const uint32_t & max_pair_distance_in);

		vector<AlignmentPath> findAlignmentPaths(const AlignmentType & alignment) const;
		vector<AlignmentPath> findAlignmentPathsIds(const AlignmentType & alignment) const;

		vector<AlignmentPath> findPairedAlignmentPaths(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;
		vector<AlignmentPath> findPairedAlignmentPathsIds(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;

	private:

        vector<uint32_t> node_seq_lengths;
       	const gbwt::GBWT & paths_index;

       	uint32_t max_pair_distance;

		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment) const;
		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::Alignment & alignment, const pair<uint32_t, uint32_t> offset) const;
		void extendAlignmentPath(AlignmentPath * align_path, const vg::Path & path, const uint32_t offset) const;

		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment) const;
		vector<AlignmentPath> extendAlignmentPath(const AlignmentPath & align_path, const vg::MultipathAlignment & alignment, const pair<uint32_t, uint32_t> offset) const;
		void extendAlignmentPaths(vector<AlignmentPath> * align_paths, const google::protobuf::RepeatedPtrField<vg::Subpath> & sub_path, const pair<uint32_t, uint32_t> offset) const;

		void pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const AlignmentType & end_alignment) const;

		multimap<gbwt::node_type, pair<int32_t, int32_t> > getAlignmentNodeIndex(const vg::Alignment & alignment) const;
		multimap<gbwt::node_type, pair<int32_t, int32_t> > getAlignmentNodeIndex(const vg::MultipathAlignment & alignment) const;	

		vg::Mapping getMapping(const vg::Alignment & alignment, const pair<uint32_t, uint32_t> offset) const;
		vg::Mapping getMapping(const vg::MultipathAlignment & alignment, const pair<uint32_t, uint32_t> offset) const;

		void printDebug(const AlignmentType & alignment_1, const AlignmentType & alignment_2) const;
};


#endif
