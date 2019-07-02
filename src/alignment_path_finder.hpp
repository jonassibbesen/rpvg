
#ifndef VGPROB_ALIGNMENTPATHFINDER_HPP
#define VGPROB_ALIGNMENTPATHFINDER_HPP

#include <vector>
#include <map>

#include <gbwt/gbwt.h>
#include <vg/io/basic_stream.hpp>

#include "alignment_path.hpp"

using namespace std;

class AlignmentPathFinder {

    public: 
    
       	AlignmentPathFinder(const vg::Graph & graph, const gbwt::GBWT & paths_index_in);

        vector<uint32_t> node_seq_lengths;
       	const gbwt::GBWT & paths_index;

		AlignmentPath findAlignmentPath(const vg::Alignment & alignment) const;
		AlignmentPath findAlignmentPathIds(const vg::Alignment & alignment) const;

		vector<AlignmentPath> findPairedAlignmentPaths(const vg::Alignment & alignment_1, const vg::Alignment & alignment_2, const int32_t max_pair_distance) const;
		vector<AlignmentPath> findPairedAlignmentPathsIds(const vg::Alignment & alignment_1, const vg::Alignment & alignment_2, const int32_t max_pair_distance) const;

	protected:

		void extendAlignmentPath(AlignmentPath * align_path, const vg::Path & path, const uint32_t & node_offset) const;
		void pairAlignmentPaths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const vg::Alignment & end_alignment, const int32_t max_pair_distance) const;
		multimap<gbwt::node_type, pair<int32_t, int32_t> > getAlignmentNodeIndex(const vg::Alignment & alignment) const;
};

#endif

