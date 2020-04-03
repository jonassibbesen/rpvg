
#include "path_abundance_writer.hpp"

#include <iomanip>


PathAbundanceWriter::PathAbundanceWriter(const bool use_stdout_in, const string filename, const double precision) : use_stdout(use_stdout_in), precision_digits(ceil(-1 * log10(precision))) {

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
    *writer_stream << "Name\tGroupID\tLength\tEffectiveLength\tConfidence\tGroupExpression\tReadCount\tTPM" << endl;
}

PathAbundanceWriter::~PathAbundanceWriter() {

	delete writer_stream;

    if (!use_stdout) {

        writer_file.close();
    }
}

void PathAbundanceWriter::writeThreadedPathClusterAbundances(const vector<vector<PathAbundances> > & threaded_path_cluster_abundances) {

    double transcript_count_sum = 0;

    for (auto & path_cluster_abundances: threaded_path_cluster_abundances) {

        for (auto & path_abundances: path_cluster_abundances) {

            for (size_t i = 0; i < path_abundances.names.size(); ++i) {

                assert(path_abundances.names.size() == path_abundances.lengths.size());
                assert(path_abundances.names.size() == path_abundances.effective_lengths.size());
                assert(path_abundances.names.size() == path_abundances.abundances.confidence.cols());
                assert(path_abundances.names.size() == path_abundances.abundances.expression.cols());

                if (path_abundances.effective_lengths.at(i) > 0) {

                    transcript_count_sum += path_abundances.abundances.expression(0, i) * path_abundances.abundances.read_count / path_abundances.effective_lengths.at(i);
                }
            }
        }
    }

    uint32_t group_id = 0;

    for (auto & path_cluster_abundances: threaded_path_cluster_abundances) {

        for (auto & path_abundances: path_cluster_abundances) {

            ++group_id;

            for (size_t i = 0; i < path_abundances.names.size(); ++i) {

                double transcript_count = 0;

                if (path_abundances.effective_lengths.at(i) > 0) {

                    transcript_count = path_abundances.abundances.expression(0, i) * path_abundances.abundances.read_count / path_abundances.effective_lengths.at(i);
                }

                *writer_stream << path_abundances.names.at(i);
                *writer_stream << "\t" << group_id;
                *writer_stream << "\t" << setprecision(3) << path_abundances.lengths.at(i);
                *writer_stream << "\t" << path_abundances.effective_lengths.at(i);
                *writer_stream << "\t" << setprecision(precision_digits) << path_abundances.abundances.confidence(0, i);
                *writer_stream << "\t" << path_abundances.abundances.expression(0, i);
                *writer_stream << "\t" << path_abundances.abundances.expression(0, i) * path_abundances.abundances.read_count;
                *writer_stream << "\t" << transcript_count / transcript_count_sum * pow(10, 6);
                *writer_stream << endl;
            }
        }
    }
}

