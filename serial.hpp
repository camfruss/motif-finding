#pragma once

#include "algorithm"
#include "array"
#include "cassert"
#include "cmath"
#include "string"
#include "unordered_map"
#include "vector"

#include "gibbs_sampler.hpp"
#include "utility.hpp"
#include "data.hpp"

#include "iostream"

template <typename T>
class Serial : public GibbsSampler<T> {
    public: 
        Serial(const Data& data);

        Result find_motifs(int k, T pseudocount) override;
};

template <typename T>
Serial<T>::Serial(const Data& data) : GibbsSampler<T>(data) {}

template <typename T>
Result Serial<T>::find_motifs(int k, T pseudocount)
{
    auto [num_sequences, sequence_length] { this->m_data.size() };

    std::vector<int> positions { this->init_positions(k) }; 
    std::vector<T> pwm { this->init_pwm(positions, k, pseudocount) };

	std::vector<T> tmp_pwm {};

    auto has_converged = [this, &pwm, &tmp_pwm, &k](const int max_iters = 10'000, const int stable_consensus = 200) {
        static int iter_count {};
		static int iters_since_change {};
		
		iters_since_change = this->consensus(pwm, k) == this->consensus(tmp_pwm, k) ?
			iters_since_change + 1 :
			0;
		// std::cout << "iters_since_change: " << iters_since_change << "iter_count: " << iter_count << std::endl;
	
        return iter_count++ > max_iters; // || iters_since_change > stable_consensus;
    };

    // for (auto s : pwm) {
    //     std::cout << s << " ";
    // }
    // std::cout << std::endl;

    int withheld { 0 }; // TODO: should consider 1 iter 1 full iteration through all sequences
    this->update_counts(pwm, withheld, positions[withheld], k, pseudocount, false); 

    do {
        std::vector<T> scores { this->score(pwm, k, withheld) };

		positions[withheld] = this->sample(scores);

        int new_withheld { (withheld + 1) % num_sequences }; 

        tmp_pwm = this->update_pwm(pwm, positions, k, pseudocount, withheld, new_withheld);
		withheld = new_withheld;
		// std::cout << "num correct pos: " << num_correct(positions, k) << std::endl;
    } while (!has_converged());

    Result result {
        .positions = positions,
	    .num_correct = this->num_correct(positions, k),
	    .consensus = this->consensus(pwm, k)
    };
    return result;
}
