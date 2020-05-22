
/*
The following code have been copied and modified from BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)
*/

/*
DiscreteSampler.cpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


The MIT License (MIT)

Copyright (c) 2016 Jonas Andreas Sibbesen and Lasse Maretty

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#include <limits>

#include "discrete_sampler.hpp"
#include "utils.hpp"


DiscreteSampler::DiscreteSampler(const uint32_t size_guess) {

	cumulative_probs.reserve(size_guess);
}

uint32_t DiscreteSampler::size() const {

	return cumulative_probs.size();
}

void DiscreteSampler::addOutcome(const double prob) {

	assert(std::isfinite(prob));

	if (cumulative_probs.empty()) {

		cumulative_probs.emplace_back(prob);
	
	} else {

		cumulative_probs.emplace_back(prob + cumulative_probs.back());
	}
}

uint32_t DiscreteSampler::sample(mt19937 * mt_rng) const {

	assert(!cumulative_probs.empty());

	double random01_sample_scaled = generate_canonical<double,std::numeric_limits<double>::digits>(*mt_rng) * cumulative_probs.back();	
	return search(random01_sample_scaled);
}

uint32_t DiscreteSampler::search(const double random01_sample_scaled) const {

	assert(std::isfinite(random01_sample_scaled));

	if (cumulative_probs.size() > 1) {

		auto cumulative_probs_it = upper_bound(cumulative_probs.begin(), cumulative_probs.end(), random01_sample_scaled);
		assert(cumulative_probs_it != cumulative_probs.end());

		uint sample_idx = cumulative_probs_it - cumulative_probs.begin();
		assert(sample_idx >= 0);

		return sample_idx;

	} else {

		return 0;
	}
}

LogDiscreteSampler::LogDiscreteSampler(const uint32_t size_guess) : DiscreteSampler(size_guess) {}

void LogDiscreteSampler::addOutcome(const double log_prob) {

	assert(std::isfinite(log_prob));

	if (cumulative_probs.empty()) {

		cumulative_probs.emplace_back(log_prob);

	} else {

		cumulative_probs.emplace_back(add_log(log_prob, cumulative_probs.back())); 
	}
}

uint32_t LogDiscreteSampler::sample(mt19937 * mt_rng) const {

	assert(!cumulative_probs.empty());

	double random01_log_sample_scaled = log(generate_canonical<double,std::numeric_limits<double>::digits>(*mt_rng)) + cumulative_probs.back();
	return search(random01_log_sample_scaled);
}

