
#include "probability_matrix_writer.hpp"

#include <iomanip>


ProbabilityMatrixWriter::ProbabilityMatrixWriter(const bool use_stdout_in, const string filename, const double precision_in) : use_stdout(use_stdout_in), precision(precision_in), num_digits(ceil(-1 * log10(precision))) {

    streambuf * writer_buffer;

    if (use_stdout) {

        writer_buffer = cout.rdbuf();
    
    } else {

        writer_file.open(filename);
        assert(writer_file.is_open());

        writer_buffer = writer_file.rdbuf();
    }

    writer_stream = new ostream(writer_buffer);

    *writer_stream << fixed;
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

bool ProbabilityMatrixWriter::collapseReadPathProbabilities(const ReadPathProbabilities & cluster_probs_1, const ReadPathProbabilities & cluster_probs_2) const {

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

void ProbabilityMatrixWriter::writeCollapsedProbabilities(const vector<pair<double, vector<uint32_t> > > & collpased_probs) {

    for (auto & prob: collpased_probs) {

        *writer_stream << " " << prob.first << ":";
        bool is_first = true;

        for (auto idx: prob.second) {

            if (is_first) {

                *writer_stream << idx;
                is_first = false;

            } else {

                *writer_stream << "," << idx;
            }
        }
    }
}

void ProbabilityMatrixWriter::writeReadPathProbabilityCluster(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<string> & path_names, const vector<double> & path_lengths) {

    if (!cluster_probs.empty()) {

        assert(!path_names.empty());
        assert(path_names.size() == path_lengths.size());

    	assert(cluster_probs.front().first.read_path_probs.size() == path_names.size());

        *writer_stream << "#" << endl << setprecision(3);
        *writer_stream << path_names.front() << "," << path_lengths.front();

        for (size_t i = 1; i < path_names.size(); ++i) {

            *writer_stream << " " << path_names.at(i) << "," << path_lengths.at(i);
        }

        *writer_stream << endl << setprecision(num_digits);

        uint32_t read_count = cluster_probs.front().second;
        uint32_t prev_unique_probs_idx = 0;

        for (size_t i = 1; i < cluster_probs.size(); ++i) {

            if (collapseReadPathProbabilities(cluster_probs.at(prev_unique_probs_idx).first, cluster_probs.at(i).first)) {

                read_count += cluster_probs.at(i).second;
            
            } else {

                *writer_stream << read_count << " " << cluster_probs.at(prev_unique_probs_idx).first.noise_prob << " ";
                writeCollapsedProbabilities(cluster_probs.at(prev_unique_probs_idx).first.collapsedProbabilities(precision));
                *writer_stream << endl;

                read_count = cluster_probs.at(i).second;
                prev_unique_probs_idx = i;
            }
        }

        *writer_stream << read_count << " " << cluster_probs.at(prev_unique_probs_idx).first.noise_prob << " ";
        writeCollapsedProbabilities(cluster_probs.at(prev_unique_probs_idx).first.collapsedProbabilities(precision));
        *writer_stream << endl;
    }
}

