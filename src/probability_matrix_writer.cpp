
#include "probability_matrix_writer.hpp"

#include <iomanip>


ProbabilityMatrixWriter::ProbabilityMatrixWriter(const bool use_stdout_in, const string filename) : use_stdout(use_stdout_in) {

    streambuf * writer_buffer;

    if (use_stdout) {

        writer_buffer = cout.rdbuf();
    
    } else {

        writer_file.open(filename);
        assert(writer_file.is_open());

        writer_buffer = writer_file.rdbuf();
    }

    writer_stream = new ostream(writer_buffer);

    *writer_stream << fixed << setprecision(ceil(-1 * log10(probability_precision)));
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

void ProbabilityMatrixWriter::writeReadPathProbabilities(const vector<pair<ReadPathProbabilities, uint32_t> > & read_path_probs, const vector<string> & path_names) {

    if (!read_path_probs.empty()) {

    	assert(read_path_probs.front().first.read_path_probs.size() == path_names.size());

        *writer_stream << "#\nx Noise";

        for (auto & name: path_names) {

            *writer_stream << " " << name;
        }

        *writer_stream << endl;

        uint32_t read_count = read_path_probs.front().second;
        uint32_t prev_unique_probs_idx = 0;

        for (size_t j = 1; j < read_path_probs.size(); ++j) {

    		assert(read_path_probs.at(j).first.read_path_probs.size() == path_names.size());

            if (read_path_probs.at(prev_unique_probs_idx).first == read_path_probs.at(j).first) {

                read_count += read_path_probs.at(j).second;
            
            } else {

                *writer_stream << read_count << " " << read_path_probs.at(prev_unique_probs_idx).first << endl;
                read_count = read_path_probs.at(j).second;
                prev_unique_probs_idx = j;
            }
        }

        *writer_stream << read_count << " " << read_path_probs.at(prev_unique_probs_idx).first << endl;
    }
}

