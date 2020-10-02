
#include "fragment_length_dist.hpp"

#include <sstream>
#include <string>

#include "vg/io/protobuf_iterator.hpp"
#include "utils.hpp"

static const uint32_t frag_length_buffer_size = 1000;
static const uint32_t max_length_sd_multiplicity = 10;

FragmentLengthDist::FragmentLengthDist() : mean_(0), sd_(1) {

    assert(isValid());
    setMaxLength();
}

FragmentLengthDist::FragmentLengthDist(const double mean_in, const double sd_in) : mean_(mean_in), sd_(sd_in) {

    assert(isValid());

    setMaxLength();
    setLogProbBuffer(frag_length_buffer_size);
}

FragmentLengthDist::FragmentLengthDist(istream * alignments_istream, const bool is_multipath) {

    assert(alignments_istream->good());

    if (is_multipath) {

        for (vg::io::ProtobufIterator<vg::MultipathAlignment> alignment_it(*alignments_istream); alignment_it.has_current(); ++alignment_it) {

            if (parseMultipathAlignment(*alignment_it)) {

                break;
            }
        }

    } else {

        for (vg::io::ProtobufIterator<vg::Alignment> alignment_it(*alignments_istream); alignment_it.has_current(); ++alignment_it) {

            if (parseAlignment(*alignment_it)) {

                break;
            }        
        }
    }

    assert(isValid());

    setMaxLength();
    setLogProbBuffer(frag_length_buffer_size);
}

FragmentLengthDist::FragmentLengthDist(const spp::sparse_hash_map<vector<AlignmentPath>, uint32_t> & align_paths_index, const uint32_t num_threads) {

    vector<uint32_t> frag_length_count_buffer;

    #pragma omp parallel num_threads(num_threads)
    {
        vector<uint32_t> threaded_frag_length_count_buffer;

        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < num_threads; ++i) {

            uint32_t cur_align_paths_index_pos = 0;

            for (auto & align_paths: align_paths_index) {

                if (cur_align_paths_index_pos % num_threads == i) { 

                    assert(!align_paths.first.empty());

                    uint32_t cur_frag_length = align_paths.first.front().seq_length;
                    bool cur_frag_length_is_constant = true;

                    for (size_t j = 1; j < align_paths.first.size(); ++j) {

                        if (align_paths.first.at(j).seq_length != cur_frag_length) {

                            cur_frag_length_is_constant = false;
                            break;
                        }
                    }

                    if (cur_frag_length_is_constant) {

                        if (threaded_frag_length_count_buffer.size() <= cur_frag_length) {
                            
                            threaded_frag_length_count_buffer.resize(cur_frag_length + 1, 0);
                        }

                        threaded_frag_length_count_buffer.at(cur_frag_length) += align_paths.second;
                    }   
                }

                ++cur_align_paths_index_pos;
            }
        }

        #pragma omp critical
        {
            if (frag_length_count_buffer.size() < threaded_frag_length_count_buffer.size()) {
                
                frag_length_count_buffer.resize(threaded_frag_length_count_buffer.size(), 0);
            }

            for (size_t i = 0; i < threaded_frag_length_count_buffer.size(); ++i) {

                frag_length_count_buffer.at(i) += threaded_frag_length_count_buffer.at(i);
            }
        }
    }

    uint32_t total_count = 0;
    uint64_t sum_count = 0; 

    for (size_t i = 0; i < frag_length_count_buffer.size(); ++i) {

        total_count += frag_length_count_buffer.at(i);
        sum_count += (i * frag_length_count_buffer.at(i));
    }

    mean_ = sum_count / static_cast<double>(total_count);

    if (total_count > 1) {

        double sum_var = 0; 

        for (size_t i = 0; i < frag_length_count_buffer.size(); ++i) {

            sum_var += (pow(static_cast<double>(i) - mean_, 2) * frag_length_count_buffer.at(i));
        }    

        sd_ = sqrt(sum_var / static_cast<double>(total_count - 1));

        if (total_count < 1000) {

            cerr << "WARNING: Only " << total_count << " unambiguous read pairs available to re-estimate fragment length distribution parameters from alignment paths. Consider setting --frag-mean and --frag-sd instead." << endl;
        }
    
        assert(isValid());

        setMaxLength();
        setLogProbBuffer(frag_length_count_buffer.size());

    } else {

        sd_ = 0;
    }
}

bool FragmentLengthDist::parseAlignment(const vg::Alignment & alignment) {

    if (alignment.fragment_length_distribution().size() > 0 && alignment.fragment_length_distribution().substr(0,1) != "0") {

        stringstream frag_length_ss = stringstream(alignment.fragment_length_distribution());
        string element;

        getline(frag_length_ss, element, ':');
        assert(stod(element) > 0);

        getline(frag_length_ss, element, ':');
        mean_ = stod(element);

        getline(frag_length_ss, element, ':');;
        sd_ = stod(element);

        return true;     
    }

    return false;
}

bool FragmentLengthDist::parseMultipathAlignment(const vg::MultipathAlignment & alignment) {

    if (alignment.has_annotation() && alignment.annotation().fields().count("fragment_length_distribution")) {

        stringstream frag_length_ss = stringstream(alignment.annotation().fields().at("fragment_length_distribution").string_value());
        string element;

        getline(frag_length_ss, element, ' ');
        assert(element == "-I");

        getline(frag_length_ss, element, ' ');
        mean_ = stod(element);

        getline(frag_length_ss, element, ' ');
        assert(element == "-D");

        getline(frag_length_ss, element);
        sd_ = stod(element);

        return true;     
    }

    return false;
}

double FragmentLengthDist::mean() const {

    return mean_;
}

double FragmentLengthDist::sd() const {

    return sd_;
}

bool FragmentLengthDist::isValid() const {

    return (mean_ >= 0 && sd_ > 0);
}

uint32_t FragmentLengthDist::maxLength() const {

    assert(max_length_ > 0);
    return max_length_;
}

double FragmentLengthDist::logProb(const uint32_t value) const {

    if (value < log_prob_buffer.size()) {

        return log_prob_buffer.at(value);
    
    } else {

        return log_normal_pdf<double>(value, mean_, sd_);
    }
}

void FragmentLengthDist::setMaxLength() {

    assert(isValid());

    max_length_ = ceil(mean_ + sd_ * max_length_sd_multiplicity);
    assert(max_length_ > 0);
}

void FragmentLengthDist::setLogProbBuffer(const uint32_t size) {

    assert(isValid());

    log_prob_buffer = vector<double>(size);

    for (size_t i = 0; i < size; ++i) {

        log_prob_buffer.at(i) = log_normal_pdf<double>(i, mean_, sd_);
    }
}


