#pragma once

#include "algorithm"
#include "array"
#include "cassert"
#include "cmath"
#include "string"
#include "vector"

#include "serial.hpp"
#include "utility.hpp"
#include "data.hpp"

#include "iostream"

template <typename T>
class Serial {
    public: 
        Serial(const Data& data);

        [[nodiscard]] std::vector<std::size_t> find_motifs(
			std::size_t k,
			T pseudocount
		);

    private:
        const Data m_data;
        const std::array<T, 4> m_background;

		T sum_log_probs(T a, T b);
        std::array<T, 4> calculate_noise(std::size_t sample_size = 100);
        std::vector<std::size_t> init_positions(std::size_t end_buffer);
       
		void update_counts(
			std::vector<T>& pwm, 
			std::size_t seq_index,  // which sequence to consider
			int k, 					// target motif length
			T pseudocount,
			std::size_t start_pos,	// current idx where motif supposedly starts
			bool increment = true
		);

		std::vector<T> init_pwm(
			int k,
			T pseudocount,
			std::size_t withheld,
			std::vector<std::size_t>& positions
		);


		void update_pwm(
			std::vector<T>& pwm, 
			std::vector<std::size_t>& positions,
			int k,
			T pseudocount,
			std::size_t old_withheld,
			std::size_t new_withheld
		);

		std::vector<T> score(
			std::vector<T>& pwm,
			int k, 
			std::size_t withheld
		);

		std::size_t sample(std::vector<T> scores);
};

template <typename T>
Serial<T>::Serial(const Data& data) 
    : m_data { data },
  	  m_background { calculate_noise() }
{
}

template <typename T>
std::vector<std::size_t> Serial<T>::find_motifs(std::size_t k, T pseudocount)
{
    auto [num_sequences, sequence_length] { m_data.size() };

    std::vector<std::size_t> positions { init_positions(k) }; 
    std::size_t withheld { 0 }; 
    std::vector<T> pwm { init_pwm(k, pseudocount, withheld, positions) };

    // TODO: experiment with different convergence criterion...
    // scoring metric -> sum of scores of complete iteration
    // fixed number of iterations / max iters
    // when motifs remain unchanged / only a few change
    // scoring metric + max_iterations
    // can add correctness plots as well
    int max_iters { 5'000 }; // TODO: this should be a parameter
    auto has_converged = [&max_iters]() {
        static std::size_t iter_count {}; 
        return iter_count++ > max_iters;
    };

    do {
        std::vector<T> scores { score(pwm, k, withheld) };
        positions[withheld] = sample(scores);

        std::size_t new_withheld { (withheld + 1) % num_sequences }; 
        update_pwm(pwm, positions, k, pseudocount, withheld, new_withheld);
		withheld = new_withheld;
    } while (!has_converged());

    return positions;
}

template <typename T>
T Serial<T>::sum_log_probs(T a, T b)
{
	return a > b ?
		a + std::log1p(std::exp(b-a)) :
		b + std::log1p(std::exp(a-b));
}

template <typename T>
std::array<T, 4> Serial<T>::calculate_noise(std::size_t sample_size)
{
    // sample 100 positions with replacement
    auto [num_sequences, sequence_length] { m_data.size() };
    std::array<T, 4> result; 

    char tmp;
    std::size_t samples_per_seq { (sample_size + num_sequences - 1) / num_sequences };
    const std::vector<Sequence> seqs { m_data.sequences() };
    for (std::size_t i {}; i < num_sequences; ++i) {
        for (std::size_t j {}; j < samples_per_seq; ++j) {
            auto idx = utility::rand_indices(sequence_length)[0]; 
            tmp = seqs[i].m_sequence[idx];    
            ++result[utility::encode(tmp)];
        }
	}

    // normalize to get probability distribution
    std::size_t total_samples { samples_per_seq * num_sequences };
    std::transform(begin(result), end(result), begin(result), [&total_samples](T x) {
        return x / total_samples;
    });

    return result;
}

template <typename T>
std::vector<std::size_t> Serial<T>::init_positions(std::size_t end_buffer)
{
	auto [num_sequences, sequence_length] { m_data.size() };
	std::vector<std::size_t> result(num_sequences);

	std::generate(begin(result), end(result), [&]() {
		return utility::rand_indices(sequence_length, end_buffer)[0]; 
	});

    return result; 
}

template <typename T>
void Serial<T>::update_counts(
	std::vector<T>& pwm,
	std::size_t seq_index,
	int k, 
	T pseudocount,
	std::size_t start_pos,
	bool increment
) {
	T delta = (increment ? 1 : -1) * 1 / (k + 4 * pseudocount) ;

	const auto& seq { m_data.sequences()[seq_index].m_sequence };
	for (std::size_t i {}; i < k; ++i) {
		std::size_t idx { 4*i + utility::encode(seq[i+start_pos]) };
		pwm[idx] += delta;    
	}
}

template <typename T>
std::vector<T> Serial<T>::init_pwm(
	int k, 
	T pseudocount, 
	std::size_t withheld,
	std::vector<std::size_t>& positions
) {
	/*
	In PWM, each nucleotide increases weight by 1/(k+4*pseudo)
	and we start with pseduo/(k+4*pseudo)
	*/
	T normalized_default { pseudocount / (k + 4*pseudocount) };
    std::vector<T> pwm(4*k, normalized_default);

    assert(m_data.sequences().size() == positions.size());
    for (std::size_t i {}; i < positions.size(); ++i) {
		if (i != withheld) {
			update_counts(pwm, i, k, pseudocount, positions[i]);
		}
	}

    return pwm;
}

template <typename T>
void Serial<T>::update_pwm(
	std::vector<T>& pwm, 
	std::vector<std::size_t>& positions,
	int k,
	T pseudocount,
	std::size_t old_withheld,
	std::size_t new_withheld
) {
	update_counts(pwm, old_withheld, k, pseudocount, positions[old_withheld]);
	update_counts(pwm, new_withheld, k, pseudocount, positions[new_withheld], false);
}

// O(seq_len * k)
template <typename T>
std::vector<T> Serial<T>::score(
	std::vector<T>& pwm,
	int k,
	std::size_t withheld
) {
	// TODO: if very slow, add thresholding, where only sample if score > some value
    auto [num_sequences, sequence_length] { m_data.size() };
	std::vector<T> score(sequence_length-k);

	const auto& seq { m_data.sequences()[withheld].m_sequence };
	for (std::size_t i {}; i < sequence_length-k; ++i) {  // iterate over possible starting positions
		T tmp {};
		for (std::size_t j { }; j < k; ++j) {  // iterates over single kmer
			std::size_t nucleotide_encoding { utility::encode(seq[i+j]) };
			tmp += 
				std::log(pwm[4*j + nucleotide_encoding]) -
				std::log(m_background[nucleotide_encoding]); 
		}

		score[i] = tmp;
	}

	T norm_factor {
		std::accumulate(std::next(begin(score)), end(score), *begin(score), [this](T a, T b) {
			return this->sum_log_probs(a, b);
		})
	};

	std::transform(begin(score), end(score), begin(score), [&norm_factor](const T& x) {
		return std::exp(x - norm_factor);
	});
	
	return score;
}

template <typename T>
std::size_t Serial<T>::sample(std::vector<T> scores) {
    std::random_device rand_device {};
	std::mt19937 num_gen { rand_device() };

	std::discrete_distribution<std::size_t> discrete_distr(begin(scores), end(scores));
	return discrete_distr(num_gen);
}

