
#include "fragment_length_dist.hpp"

#include <sstream>
#include <string>
#include <cmath>

#include "vg/io/protobuf_iterator.hpp"
#include "utils.hpp"

//#define debug_skew_normal_fit


FragmentLengthDist::FragmentLengthDist() : loc_(0), scale_(0), shape_(0), max_length_(0) {

    assert(!isValid());
}

FragmentLengthDist::FragmentLengthDist(const double mean_in, const double sd_in, const uint32_t sd_max_multi) : FragmentLengthDist(mean_in, sd_in, 0.0, sd_max_multi) {}

FragmentLengthDist::FragmentLengthDist(const double loc_in, const double scale_in, const double shape_in, const uint32_t sd_max_multi) : loc_(loc_in), scale_(scale_in), shape_(shape_in) {

    assert(isValid());

    setMaxLength(sd_max_multi);
    setLogProbBuffer(max_length_);
}

FragmentLengthDist::FragmentLengthDist(istream * alignments_istream, const bool is_multipath, const uint32_t sd_max_multi) {

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

    setMaxLength(sd_max_multi);
    setLogProbBuffer(max_length_);
}

FragmentLengthDist::FragmentLengthDist(const vector<uint32_t> & frag_length_counts, const bool skew_normal) {

    assert(!frag_length_counts.empty());
    assert(frag_length_counts.front() == 0);

    uint32_t sample_size = 0;
    uint64_t frag_length_sum = 0;

    for (size_t i = 0; i < frag_length_counts.size(); ++i) {
        
        sample_size += frag_length_counts.at(i);
        frag_length_sum += (i * frag_length_counts.at(i));
    }

    if (sample_size < 2) {

        loc_ = frag_length_sum;
        scale_ = 0.0;
        shape_ = 0.0;
 
        assert(!isValid());
   
    } else {

        if (sample_size < 1000) {

            cerr << "WARNING: Only " << sample_size << " unambiguous read pairs available to re-estimate fragment length distribution parameters from alignment paths." << endl;
        }

        if (!skew_normal) {
            
            loc_ = frag_length_sum / static_cast<double>(sample_size);
                        
            double sum_var = 0;
            
            for (size_t i = 0; i < frag_length_counts.size(); ++i) {
                
                sum_var += (pow(static_cast<double>(i) - loc_, 2) * frag_length_counts.at(i));
            }
            
            scale_ = sqrt(sum_var / static_cast<double>(sample_size - 1));
            shape_ = 0.0;

        } else {
            
            // compute the cumulants
            double k0 = sample_size;
            double k1 = frag_length_sum;
            double k2 = 0.0;
            double k3 = 0.0;

            for (size_t i = 0; i < frag_length_counts.size(); ++i) {
                double term = frag_length_counts.at(i) * i * i;
                k2 += term;
                term *= i;
                k3 += term;
            }
            
    #ifdef debug_skew_normal_fit
            cerr << "non-central moments:" << endl;
            cerr << "\tk0 " << k0 << endl;
            cerr << "\tk1 " << k1 << endl;
            cerr << "\tk2 " << k2 << endl;
            cerr << "\tk3 " << k3 << endl;
    #endif
            
            // convert to un-normalized central moments
            double m1 = k1 / k0;
            double m2 = k2 / k0 - m1 * m1;
            double m3 = k3 / k0 - 3.0 * m1 * m2 - m1 * m1 * m1;
            
    #ifdef debug_skew_normal_fit
            cerr << "central moments:" << endl;
            cerr << "\tm1 " << m1 << endl;
            cerr << "\tm2 " << m2 << endl;
            cerr << "\tm3 " << m3 << endl;
    #endif
            
            // convert to normalized central moments
            double mean = m1;
            double sd = sqrt(m2);
            double skew = m3 / (sd * sd * sd);
            
    #ifdef debug_skew_normal_fit
            cerr << "normalized central moments:" << endl;
            cerr << "\tmean " << mean << endl;
            cerr << "\tsd " << sd << endl;
            cerr << "\tskew " << skew << endl;
    #endif
            
            // use method of moments to derive starting point for MLE
            double alpha = 0.0;
            double sigma = 0.0;
            double mu = 0.0;
            if (skew != 0.0 && k0 > 2.0) {
                // make sure we don't shoot past the theoretical maximum skew
                // with the sample skew (which has a non-zero probability)
                double gam = pow(min<double>(abs(skew), 0.9952717464311565), 2.0 / 3.0);
                double abs_delta = sqrt((Utils::pi / 2.0) * (gam / (gam + pow((4.0 - Utils::pi) / 2.0, 2.0 / 3.0))));
                double abs_alpha = abs_delta / sqrt(1.0 - abs_delta * abs_delta);
                alpha = skew < 0.0 ? -abs_alpha : abs_alpha;
            }
            double delta = alpha / sqrt(1.0 + alpha * alpha);
            if (sd != 0.0 && k0 > 1.0) {
                sigma = sd / sqrt(1.0 - 2.0 * delta * delta / Utils::pi);
            }
            mu = mean - sigma * delta * sqrt(2.0 / Utils::pi);
            
    #ifdef debug_skew_normal_fit
            cerr << "MOM estimates:" << endl;
            cerr << "\tmu " << mu << endl;
            cerr << "\tsigma " << sigma << endl;
            cerr << "\talpha " << alpha << endl;
    #endif
            
            // alternatingly optimize alpha and mu until convergence (sigma can
            // be computed analytically)
            auto log_likelihood = [&](double mu, double sigma, double alpha) {
                double ll = 0.0;
                for (size_t i = 0; i < frag_length_counts.size(); ++i) {
                    if (frag_length_counts.at(i) == 0) {
                        continue;
                    }
                    ll += frag_length_counts.at(i) * Utils::log_skew_normal_pdf(double(i), mu, sigma, alpha);
                }
                return ll;
            };
            function<double(double)> alpha_log_likelihood = [&](double a) {
                return log_likelihood(mu, sigma, a);
            };
            function<double(double)> mu_log_likelihood = [&](double m) {
                return log_likelihood(m, sigma, alpha);
            };
            
            double tol = 1e-4;
            double prev_mu = mu + 2.0 * tol;
            double prev_alpha = alpha + 2.0 * tol;
            int max_iters = 100;
            int iter_num = 0;
            while (iter_num < max_iters &&
                   (abs(prev_mu - mu) >= tol || abs(prev_alpha - alpha) >= tol)) {
                
                ++iter_num;
                prev_mu = mu;
                prev_alpha = alpha;
                double ll = alpha_log_likelihood(alpha);
                // find an interval that contains a local maximum
                double left_radius = 1.0, right_radius = 1.0;
                while (alpha_log_likelihood(alpha - left_radius) >= ll) {
                    left_radius *= 2.0;
                }
                while (alpha_log_likelihood(alpha + right_radius) >= ll) {
                    right_radius *= 2.0;
                }
                // maximize likelihood for alpha
                alpha = Utils::golden_section_search<double>(alpha_log_likelihood,
                                                             alpha - left_radius, alpha + right_radius, tol / 4.0);
                
                // repeat for mu
                ll = mu_log_likelihood(mu);
                left_radius = 1.0;
                right_radius = 1.0;
                while (mu_log_likelihood(alpha - left_radius) >= ll) {
                    left_radius *= 2.0;
                }
                while (mu_log_likelihood(alpha + right_radius) >= ll) {
                    right_radius *= 2.0;
                }
                mu = Utils::golden_section_search<double>(mu_log_likelihood,
                                                          mu - left_radius, mu + right_radius, tol / 4.0);
                
                // analytical solution for sigma from equation 8 of Azzalini (1985)
                double sum_sq_dev = 0.0;
                for (size_t i = 0; i < frag_length_counts.size(); ++i) {
                    double dev = i - mu;
                    sum_sq_dev += frag_length_counts.at(i) * dev * dev;
                }
                sigma = sqrt(sum_sq_dev / k0);
                
    #ifdef debug_skew_normal_fit
                cerr << "iter " << iter_num << " estimates:" << endl;
                cerr << "\tmu " << mu << endl;
                cerr << "\tsigma " << sigma << endl;
                cerr << "\talpha " << alpha << endl;
    #endif
            }

            loc_ = mu;
            scale_ = sigma;
            shape_ = alpha;
        }

        assert(isValid());

        max_length_ = frag_length_counts.size();
        setLogProbBuffer(frag_length_counts.size());
    }
}

bool FragmentLengthDist::parseAlignment(const vg::Alignment & alignment) {

    if (alignment.fragment_length_distribution().size() > 0 && alignment.fragment_length_distribution().substr(0,1) != "0") {

        stringstream frag_length_ss = stringstream(alignment.fragment_length_distribution());
        string element;

        getline(frag_length_ss, element, ':');
        assert(stod(element) > 0);

        getline(frag_length_ss, element, ':');
        loc_ = stod(element);

        getline(frag_length_ss, element, ':');;
        scale_ = stod(element);
        
        shape_ = 0.0;
        
        return true;     
    
    } else if (alignment.has_annotation() && alignment.annotation().fields().count("fragment_length_distribution")) {

        stringstream frag_length_ss = stringstream(alignment.annotation().fields().at("fragment_length_distribution").string_value());
        string element;

        getline(frag_length_ss, element, ' ');
        assert(element == "-I");

        getline(frag_length_ss, element, ' ');
        loc_ = stod(element);

        getline(frag_length_ss, element, ' ');
        assert(element == "-D");

        getline(frag_length_ss, element);
        scale_ = stod(element);

        shape_ = 0.0;
        
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
        loc_ = stod(element);

        getline(frag_length_ss, element, ' ');
        assert(element == "-D");

        getline(frag_length_ss, element);
        scale_ = stod(element);

        shape_ = 0.0;
        
        return true;     
    }

    return false;
}

double FragmentLengthDist::loc() const {

    return loc_;
}

double FragmentLengthDist::scale() const {

    return scale_;
}

double FragmentLengthDist::shape() const {
    
    return shape_;
}

bool FragmentLengthDist::isValid() const {

    return (loc_ >= 0 && scale_ > 0);
}

uint32_t FragmentLengthDist::maxLength() const {

    assert(max_length_ > 0);
    return max_length_;
}

double FragmentLengthDist::logProb(const uint32_t value) const {

    if (value < log_prob_buffer.size()) {
        return log_prob_buffer.at(value);
    } else if (Utils::doubleCompare(shape_, 0.0)) {
        return Utils::log_normal_pdf<double>(value, loc_, scale_);
    } else {
        return Utils::log_skew_normal_pdf<double>(value, loc_, scale_, shape_);
    }
}

void FragmentLengthDist::setMaxLength(const uint32_t sd_max_multi) {

    assert(isValid());

    double delta = shape_ / sqrt(1.0 + shape_ * shape_);
    double sd = scale_ * (1.0 - 2.0 * delta * delta / Utils::pi);
    
    max_length_ = ceil(loc_ + sd * sd_max_multi);
    assert(max_length_ > 0);
}

void FragmentLengthDist::setLogProbBuffer(const uint32_t size) {

    assert(isValid());

    log_prob_buffer = vector<double>(size + 1);

    if (Utils::doubleCompare(shape_, 0.0)) {

        for (size_t i = 0; i < size + 1; ++i) {
         
            log_prob_buffer.at(i) = Utils::log_normal_pdf<double>(i, loc_, scale_);
        }

    } else {

        for (size_t i = 0; i < size + 1; ++i) {

            log_prob_buffer.at(i) = Utils::log_skew_normal_pdf<double>(i, loc_, scale_, shape_);
        }
    }
}


