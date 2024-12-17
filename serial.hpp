#pragma once

#include "algorithm"
#include "array"
#include "cassert"
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

        std::array<T, 4> calculate_noise(std::size_t sample_size = 100);
        std::vector<std::size_t> init_positions(std::size_t end_buffer);
       
		void update_counts(
			std::vector<T>& pwm, 
			int k, 
			std::size_t position,
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

	std::cout << "Init positions" << std::endl;
    std::vector<std::size_t> positions { init_positions(k) }; 
    std::size_t withheld { 0 }; 
	std::cout << "Init pwm" << std::endl;
    std::vector<T> pwm { init_pwm(k, pseudocount, withheld, positions) };

    // TODO: experiment with different convergence criterion...
    // scoring metric -> sum of scores of complete iteration
    // fixed number of iterations / max iters
    // when motifs remain unchanged / only a few change
    // scoring metric + max_iterations
    // can add correctness plots as well
    int max_iters { 1'000 }; // TODO: this should be a parameter
    auto has_converged = [&max_iters]() {
        static std::size_t iter_count {}; 
        return iter_count++ > max_iters;
    };

	// auto normalize = std::transform(begin(pwm), end(pwm), begin(pwm), [&k](T x){
	// 	return x / (k + 4*pseudocount);
	// });

    do {
		std::cout << "Beginning scoring" << std::endl;
        std::vector<T> scores { score(pwm, k, withheld) };
		std::cout << "Beginning sampling" << std::endl;
        positions[withheld] = sample(scores);

        std::size_t new_withheld { (withheld + 1) % num_sequences }; 
		std::cout << "Beginning pwm update" << std::endl;
        update_pwm(pwm, positions, k, withheld, new_withheld);
		withheld = new_withheld;
    } while (!has_converged());

    return positions;
}

template <typename T>
std::array<T, 4> Serial<T>::calculate_noise(std::size_t sample_size)
{
    // sample 100 positions with replacement
    auto [num_sequences, sequence_length] { m_data.size() };
    std::array<T, 4> result; 

    char tmp;
    std::size_t samples_per_seq { (sample_size + num_sequences - 1) / num_sequences };
    const std::vector<std::string> seqs { m_data.sequences() };
    for (std::size_t i {}; i < num_sequences; ++i) {
        for (std::size_t j {}; j < samples_per_seq; ++j) {
            auto idx = utility::rand_indices(sequence_length)[0]; 
            tmp = seqs[i][idx];    
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
	int k, 
	std::size_t position,
	bool increment
) {
	auto seq { m_data.sequences()[position] };
	for (std::size_t i { position }; i < position + k; ++i) {
		std::size_t idx { 4*i + utility::encode(seq[i]) };
		increment ? ++pwm[idx] : --pwm[idx];    
	}
}

template <typename T>
std::vector<T> Serial<T>::init_pwm(
	int k, 
	T pseudocount, 
	std::size_t withheld,
	std::vector<std::size_t>& positions
) {
    std::vector<T> pwm(4*k, pseudocount);

    assert(m_data.sequences().size() == positions.size());
    for (std::size_t i_pos {}; i_pos < positions.size(); ++i_pos) {
		if (i_pos != withheld) {
			update_counts(pwm, k, pseudocount, i_pos);
		}
	}

    return pwm;
}

template <typename T>
void Serial<T>::update_pwm(
	std::vector<T>& pwm, 
	std::vector<std::size_t>& positions,
	int k,
	std::size_t old_withheld,
	std::size_t new_withheld
) {
	update_counts(pwm, k, positions[old_withheld]);
	update_counts(pwm, k, positions[new_withheld], false);
}

template <typename T>
std::vector<T> Serial<T>::score(
	std::vector<T>& pwm,
	int k,
	std::size_t withheld
) {
	// TODO: if very slow, add thresholding, where only sample if score > some value
	// TODO: need to work with log probabilities
    auto [num_sequences, sequence_length] { m_data.size() };
	std::vector<T> score(sequence_length-k);

	// O(seq_len * k)
	T tmp {};
	auto seq { m_data.sequences()[withheld] };
	for (std::size_t i {}; i < sequence_length - k; ++i) {  // iterate over possible starting positions
		for (std::size_t j { }; j < k; ++j) {  // iterates over single kmer
			std::size_t nucleotide_encoding { utility::encode(seq[i+j]) };
			tmp *= pwm[4*j + nucleotide_encoding] / m_background[nucleotide_encoding]; 
		}
		score[i] = tmp;
	}
	
	return score;
}

template <typename T>
std::size_t Serial<T>::sample(std::vector<T> scores) {
    std::random_device rand_device {};
	std::mt19937 num_gen { rand_device() };

	std::discrete_distribution<std::size_t> discrete_distr(begin(scores), end(scores));
	return discrete_distr(num_gen);
}

