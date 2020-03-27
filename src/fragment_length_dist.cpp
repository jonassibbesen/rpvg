
#include "fragment_length_dist.hpp"

#include <sstream>
#include <string>

#include "vg/io/protobuf_iterator.hpp"
#include "utils.hpp"


static const uint32_t max_length_sd_multiplicity = 10;

FragmentLengthDist::FragmentLengthDist() : mean_(0), sd_(1) {}

FragmentLengthDist::FragmentLengthDist(const double mean_in, const double sd_in): mean_(mean_in), sd_(sd_in) {

    assert(isValid());
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

    assert(isValid());
    return ceil(mean_ + sd_ * max_length_sd_multiplicity);
}

double FragmentLengthDist::logProb(const uint32_t value) const {

    assert(isValid());
    return log_normal_pdf<double>(value, mean_, sd_);
}

