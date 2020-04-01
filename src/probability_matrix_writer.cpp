
#include "probability_matrix_writer.hpp"

#include <iomanip>


ProbabilityMatrixWriter::ProbabilityMatrixWriter(const bool use_stdout_in, const string filename, const double precision_in) : use_stdout(use_stdout_in), precision(precision_in) {

    streambuf * writer_buffer;

    if (use_stdout) {

        writer_buffer = cout.rdbuf();
    
    } else {

        writer_file.open(filename);
        assert(writer_file.is_open());

        writer_buffer = writer_file.rdbuf();
    }

    writer_stream = new ostream(writer_buffer);

    *writer_stream << fixed << setprecision(ceil(-1 * log10(precision)));
}

ProbabilityMatrixWriter::~ProbabilityMatrixWriter() {

	delete writer_stream;

    if (!use_stdout) {

        writer_file.close();
    }
}

void ProbabilityMatrixWriter::lockWriter() {

	writer_mutex.lock();
}

void ProbabilityMatrixWriter::unlockWriter() {

	writer_mutex.unlock();
}

bool ProbabilityMatrixWriter::collapseReadPathProbabilities(const ReadPathProbabilities & cluster_probs_1, const ReadPathProbabilities & cluster_probs_2) {

    assert(cluster_probs_1.read_path_probs.size() == cluster_probs_2.read_path_probs.size());

    if (abs(cluster_probs_1.noise_prob - cluster_probs_2.noise_prob) < precision) {

        for (size_t i = 0; i < cluster_probs_1.read_path_probs.size(); ++i) {

            if (abs(cluster_probs_1.read_path_probs.at(i) - cluster_probs_2.read_path_probs.at(i)) >= precision) {

                return false;
            }
        }

        return true;
    } 

    return false;
}

void ProbabilityMatrixWriter::writeReadPathProbabilityCluster(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<string> & path_names, const vector<double> & path_lengths) {

    if (!cluster_probs.empty()) {

        assert(!path_names.empty());
        assert(path_names.size() == path_lengths.size());

    	assert(cluster_probs.front().first.read_path_probs.size() == path_names.size());

        *writer_stream << "#" << endl;
        *writer_stream << path_names.front() << "," << path_lengths.front();

        for (size_t i = 1; i < path_names.size(); ++i) {

            *writer_stream << " " << path_names.at(i) << "," << path_lengths.at(i);
        }

        *writer_stream << endl;

        uint32_t read_count = cluster_probs.front().second;
        uint32_t prev_unique_probs_idx = 0;

        for (size_t i = 1; i < cluster_probs.size(); ++i) {

            if (collapseReadPathProbabilities(cluster_probs.at(prev_unique_probs_idx).first, cluster_probs.at(i).first)) {

                read_count += cluster_probs.at(i).second;
            
            } else {

                *writer_stream << read_count << " " << cluster_probs.at(prev_unique_probs_idx).first.getCollapsedProbabilityString(precision) << endl;
                read_count = cluster_probs.at(i).second;
                prev_unique_probs_idx = i;
            }
        }

        *writer_stream << read_count << " " << cluster_probs.at(prev_unique_probs_idx).first.getCollapsedProbabilityString(precision) << endl;
    }
}

