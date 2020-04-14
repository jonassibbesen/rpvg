
/*
The following code have been copied and modified from BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)
*/

/*
DiscreteSampler.hpp - This file is part of BayesTyper (https://github.com/bioinformatics-centre/BayesTyper)


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


#ifndef FERSKEN_SRC_DISCRETESAMPLER_HPP
#define FERSKEN_SRC_DISCRETESAMPLER_HPP

#include <random>
#include <vector>

using namespace std;


class DiscreteSampler {

	public:

		DiscreteSampler(const uint32_t size_guess);
		virtual ~DiscreteSampler() {};

		uint32_t size() const;

		virtual void addOutcome(const double prob); 
		virtual uint32_t sample(mt19937 * mt_rng) const;

	protected:

		uint32_t search(const double random01_sample_scaled) const;

		vector<double> cumulative_probs;		
};

class LogDiscreteSampler : public DiscreteSampler {

	public:
		
		LogDiscreteSampler(const uint32_t size_guess);
		~LogDiscreteSampler() {};

		void addOutcome(const double log_prob);
		uint32_t sample(mt19937 * mt_rng) const;
};


#endif
