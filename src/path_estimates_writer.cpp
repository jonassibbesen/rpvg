
#include "path_estimates_writer.hpp"

#include <iomanip>


PathEstimatesWriter::PathEstimatesWriter(const bool use_stdout_in, const string filename) : use_stdout(use_stdout_in) {

    streambuf * writer_buffer;

    if (use_stdout) {

        writer_buffer = cout.rdbuf();
    
    } else {

        writer_file.open(filename);
        assert(writer_file.is_open());

        writer_buffer = writer_file.rdbuf();
    }

    writer_stream = new ostream(writer_buffer);
}

PathEstimatesWriter::~PathEstimatesWriter() {

	delete writer_stream;

    if (!use_stdout) {

        writer_file.close();
    }
}

void PathEstimatesWriter::writePathClusterPosteriors(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const uint32_t ploidy) {

    for (uint32_t i = 0; i < ploidy; ++i) {

        *writer_stream << "Name" << i + 1 << "\t";
    }

    *writer_stream << "ClusterID\tPosterior\tClusterReadCount" << endl;

    uint32_t cluster_id = 0;

    for (auto & path_cluster_estimates_thread: threaded_path_cluster_estimates) {

        for (auto & path_cluster_estimates: path_cluster_estimates_thread) {

            ++cluster_id;

            assert(path_cluster_estimates.path_groups.size() == path_cluster_estimates.posteriors.cols());

            for (size_t i = 0; i < path_cluster_estimates.path_groups.size(); ++i) {

                assert(path_cluster_estimates.path_groups.at(i).size() == ploidy);

                for (auto & path_idx: path_cluster_estimates.path_groups.at(i)) {

                    *writer_stream << path_cluster_estimates.paths.at(path_idx).name << "\t";
                }

                *writer_stream << cluster_id;
                *writer_stream << "\t" << path_cluster_estimates.posteriors(0, i);
                *writer_stream << "\t" << path_cluster_estimates.read_count;
                *writer_stream << endl;
            }
        }
    }
}

void PathEstimatesWriter::writePathClusterAbundances(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates) {

    *writer_stream << "Name\tClusterID\tLength\tEffectiveLength\tHaplotypePosterior\tClusterRelativeExpression\tReadCount\tTPM" << endl;

    double transcript_count_sum = 0;

    for (auto & path_cluster_estimates_thread: threaded_path_cluster_estimates) {

        for (auto & path_cluster_estimates: path_cluster_estimates_thread) {

            assert(path_cluster_estimates.paths.size() == path_cluster_estimates.posteriors.cols());
            assert(path_cluster_estimates.paths.size() == path_cluster_estimates.abundances.cols());
            assert(path_cluster_estimates.path_groups.empty());

            for (size_t i = 0; i < path_cluster_estimates.paths.size(); ++i) {

                if (path_cluster_estimates.paths.at(i).effective_length > 0) {

                    transcript_count_sum += path_cluster_estimates.abundances(0, i) * path_cluster_estimates.read_count / path_cluster_estimates.paths.at(i).effective_length;
                }
            }
        }
    }

    uint32_t cluster_id = 0;

    for (auto & path_cluster_estimates_thread: threaded_path_cluster_estimates) {

        for (auto & path_cluster_estimates: path_cluster_estimates_thread) {

            ++cluster_id;

            for (size_t i = 0; i < path_cluster_estimates.paths.size(); ++i) {

                double transcript_count = 0;

                if (path_cluster_estimates.paths.at(i).effective_length > 0) {

                    transcript_count = path_cluster_estimates.abundances(0, i) * path_cluster_estimates.read_count / path_cluster_estimates.paths.at(i).effective_length;
                }

                *writer_stream << path_cluster_estimates.paths.at(i).name;
                *writer_stream << "\t" << cluster_id;
                *writer_stream << "\t" << path_cluster_estimates.paths.at(i).length;
                *writer_stream << "\t" << path_cluster_estimates.paths.at(i).effective_length;
                *writer_stream << "\t" << path_cluster_estimates.posteriors(0, i);
                *writer_stream << "\t" << path_cluster_estimates.abundances(0, i);
                *writer_stream << "\t" << path_cluster_estimates.abundances(0, i) * path_cluster_estimates.read_count;
                *writer_stream << "\t" << transcript_count / transcript_count_sum * pow(10, 6);
                *writer_stream << endl;
            }
        }
    }
}

void PathEstimatesWriter::writeBestPathGroupsPosteriors(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const PathsIndex & paths_index, const bool is_call_traversals) {

    for (auto & path_cluster_estimates_thread: threaded_path_cluster_estimates) {

        for (auto & path_cluster_estimates: path_cluster_estimates_thread) {

            assert(path_cluster_estimates.path_groups.size() == path_cluster_estimates.posteriors.cols());

            double cur_max_posterior = 0;
            uint32_t cur_max_posterior_idx = 0;

            for (size_t i = 0; i < path_cluster_estimates.path_groups.size(); ++i) {

                if (path_cluster_estimates.posteriors(0, i) > cur_max_posterior) {

                    cur_max_posterior = path_cluster_estimates.posteriors(0, i);
                    cur_max_posterior_idx = i;
                }
            }

            if (cur_max_posterior > 0) {

                for (size_t i = 0; i < path_cluster_estimates.path_groups.at(cur_max_posterior_idx).size(); ++i) {

                    auto path_idx = path_cluster_estimates.path_groups.at(cur_max_posterior_idx).at(i);
                    auto gaf_path = gbwtPathToGAFPath(path_cluster_estimates.paths.at(path_idx), paths_index, is_call_traversals);

                    stringstream path_name;
                    path_name << path_cluster_estimates.paths.at(path_idx).name;
                    path_name << "_" << i + 1;

                    writeGAFline(path_name.str(), gaf_path, cur_max_posterior, path_cluster_estimates.read_count);
                }
            }
        }
    } 
}

void PathEstimatesWriter::writeMarginalPathPosteriors(const vector<vector<PathClusterEstimates> > & threaded_path_cluster_estimates, const PathsIndex & paths_index, const bool is_call_traversals) {

    for (auto & path_cluster_estimates_thread: threaded_path_cluster_estimates) {

        for (auto & path_cluster_estimates: path_cluster_estimates_thread) {

            assert(path_cluster_estimates.path_groups.size() == path_cluster_estimates.posteriors.cols());

            unordered_map<uint32_t, double> marginal_posteriors;

            for (size_t i = 0; i < path_cluster_estimates.path_groups.size(); ++i) {

                unordered_set<uint32_t> unique_path_idxs;

                for (auto & path_idx: path_cluster_estimates.path_groups.at(i)) {

                    if (unique_path_idxs.emplace(path_idx).second) {

                        auto marginal_posteriors_it = marginal_posteriors.emplace(path_idx, 0);
                        marginal_posteriors_it.first->second += path_cluster_estimates.posteriors(0, i);
                    }
                }
            }

            for (auto & path_posterior: marginal_posteriors) {

                auto gaf_path = gbwtPathToGAFPath(path_cluster_estimates.paths.at(path_posterior.first), paths_index, is_call_traversals);
                writeGAFline(path_cluster_estimates.paths.at(path_posterior.first).name, gaf_path, path_posterior.second, path_cluster_estimates.read_count);
            }
        }
    } 
}

pair<string, uint32_t> PathEstimatesWriter::gbwtPathToGAFPath(const PathInfo & path, const PathsIndex & paths_index, const bool is_call_traversals) {

    pair<string, uint32_t> gaf_path("", 0);

    auto gbwt_path_id = path.id;

    if (paths_index.index().bidirectional()) {

        gbwt_path_id = gbwt::Path::encode(path.id, false);
    }

    pair<uint32_t, uint32_t> call_snarl_boundaries(numeric_limits<uint32_t>::max(), numeric_limits<uint32_t>::max());

    if (is_call_traversals) {

        auto call_path_name_split = splitString(path.name, '_');
        assert(call_path_name_split.size() == 4);

        call_snarl_boundaries.first = stoi(call_path_name_split.at(1));
        call_snarl_boundaries.second = stoi(call_path_name_split.at(2));
    }

    auto gbwt_path = paths_index.index().extract(gbwt_path_id);
    assert(!gbwt_path.empty());

    stringstream node_sequence;

    for (auto & gbwt_node: gbwt_path) {

        auto node_id = gbwt::Node::id(gbwt_node);

        if (is_call_traversals) {

            if (node_id == call_snarl_boundaries.first) {

                call_snarl_boundaries.first = numeric_limits<uint32_t>::max();
            }

            if (call_snarl_boundaries.first != numeric_limits<uint32_t>::max()) {

                continue;
            }

            if (call_snarl_boundaries.second == numeric_limits<uint32_t>::max()) {

                continue;
            }

            if (node_id == call_snarl_boundaries.second) {

                call_snarl_boundaries.second = numeric_limits<uint32_t>::max();
            }
        }

        node_sequence << ((gbwt::Node::is_reverse(gbwt_node)) ? "<" : ">");
        node_sequence << node_id;

        gaf_path.second += paths_index.nodeLength(gbwt::Node::id(gbwt_node));
    }

    gaf_path.first = node_sequence.str();

    return gaf_path;
}

void PathEstimatesWriter::writeGAFline(const string & path_name, const pair<string, uint32_t> & gaf_path, const double path_posterior, const uint32_t read_count) {

    *writer_stream << path_name;
    *writer_stream << "\t" << gaf_path.second;
    *writer_stream << "\t" << "0";
    *writer_stream << "\t" << gaf_path.second;
    *writer_stream << "\t" << "+";
    *writer_stream << "\t" << gaf_path.first;
    *writer_stream << "\t" << gaf_path.second;
    *writer_stream << "\t" << "0";
    *writer_stream << "\t" << gaf_path.second;
    *writer_stream << "\t" << gaf_path.second;
    *writer_stream << "\t" << gaf_path.second;
    *writer_stream << "\t" << ((doubleCompare(path_posterior, 1)) ? 60 : prob_to_phred(1 - path_posterior));
    *writer_stream << "\t" << "cs:Z::" << gaf_path.second;
    *writer_stream << "\t" << "pp:f:" << path_posterior;
    *writer_stream << "\t" << "rc:i:" << read_count;
    *writer_stream << endl;
}


