
#include "probability_matrix_writer.hpp"

#include <iomanip>


ProbabilityMatrixWriter::ProbabilityMatrixWriter(const bool use_stdout_in, const string filename, const double precision_in) : use_stdout(use_stdout_in), precision(precision_in), precision_digits(ceil(-1 * log10(precision))) {

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

    assert(cluster_probs_1.probabilities().size() == cluster_probs_2.probabilities().size());

    if (abs(cluster_probs_1.noiseProbability() - cluster_probs_2.noiseProbability()) < precision) {

        for (size_t i = 0; i < cluster_probs_1.probabilities().size(); ++i) {

            if (abs(cluster_probs_1.probabilities().at(i) - cluster_probs_2.probabilities().at(i)) >= precision) {

                return false;
            }
        }

        return true;
    } 

    return false;
}

void ProbabilityMatrixWriter::writeCollapsedProbabilities(const vector<pair<double, vector<uint32_t> > > & collpased_probs, const bool write_zero) {

    for (auto & prob: collpased_probs) {

        if (write_zero || !doubleCompare(prob.first, 0)) {

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
}

void ProbabilityMatrixWriter::writeReadPathProbabilityCluster(const vector<pair<ReadPathProbabilities, uint32_t> > & cluster_probs, const vector<Path> & cluster_paths) {

    assert(!cluster_paths.empty());

    *writer_stream << "#" << endl;
    *writer_stream << setprecision(3);
    *writer_stream << cluster_paths.front().name << "," << cluster_paths.front().length << "," << cluster_paths.front().effective_length;

    for (size_t i = 1; i < cluster_paths.size(); ++i) {

        *writer_stream << " " << cluster_paths.at(i).name << "," << cluster_paths.at(i).length << "," << cluster_paths.at(i).effective_length;
    }

    *writer_stream << endl;

    if (!cluster_probs.empty()) {

        *writer_stream << setprecision(precision_digits);

        assert(cluster_probs.front().first.probabilities().size() == cluster_paths.size());

        uint32_t read_count = cluster_probs.front().second;
        uint32_t prev_unique_probs_idx = 0;

        for (size_t i = 1; i < cluster_probs.size(); ++i) {

            if (collapseReadPathProbabilities(cluster_probs.at(prev_unique_probs_idx).first, cluster_probs.at(i).first)) {

                read_count += cluster_probs.at(i).second;
            
            } else {

                *writer_stream << read_count << " " << cluster_probs.at(prev_unique_probs_idx).first.noiseProbability() << " ";
                writeCollapsedProbabilities(cluster_probs.at(prev_unique_probs_idx).first.collapsedProbabilities(precision), false);
                *writer_stream << endl;

                read_count = cluster_probs.at(i).second;
                prev_unique_probs_idx = i;
            }
        }

        *writer_stream << read_count << " " << cluster_probs.at(prev_unique_probs_idx).first.noiseProbability() << " ";
        writeCollapsedProbabilities(cluster_probs.at(prev_unique_probs_idx).first.collapsedProbabilities(precision), false);
        *writer_stream << endl;
    }
}

