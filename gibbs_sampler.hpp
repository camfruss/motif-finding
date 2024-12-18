#pragma once

#include "cassert"
#include "vector"

#include "data.hpp"

// enum class Criteria 
// {
// 	MAX_REPEAT,
// 	IS_STABLE
// };

struct Result
{
	std::vector<int> positions;
	int num_correct;
	std::string consensus;
};

template <typename T>
class GibbsSampler {
    public: 
        GibbsSampler(const Data& data);
		virtual ~GibbsSampler() = default;

        [[nodiscard]] virtual Result find_motifs(int k, T pseudocount) = 0;

    protected:
        const Data m_data;

		/* Returns the number of correctly estimated motif starting positions 
		 * Note: overlap is considered "correct"
		 */
		int num_correct(std::vector<int>& positions, int k);
		
		/* Calculates the consensus motif based on a current PWM */
		std::string consensus(std::vector<T>& pwm, int k);
        
		/* Initializes random motif starting positions for each sequence 
		 * in m_data 
		 */
		std::vector<int> init_positions(int end_buffer);

		/* Initializes a PWM from ALL sequences. 
		 * positions : motif starting positions for each sequence
		 */
		std::vector<T> init_pwm(std::vector<int>& positions, int k, 
			T pseudocount);

		std::vector<T> update_pwm(std::vector<T>& pwm, 
			std::vector<int>& positions, int k, T pseudocount, int old_withheld,
			int new_withheld);

		/* Helper function to update the PWM
		 * seq_index : sequence to use to update the PWM
		 * start_pos : current idx where the motif is estimated to start
		 * increment : true if adding to the PWM, false if removing
		 */
		void update_counts(std::vector<T>& pwm, int seq_index, int start_pos, 
			int k, T pseudocount, bool increment = true);

		/* Scores each k-mer in the withheld sequence using the PWM 
		 * Returns a probability distribution
		 */
		std::vector<T> score(std::vector<T>& pwm, int k, int withheld);

		/* Samples index space covered by the prob distribution in scores */
		int sample(std::vector<T> scores);

	private:
        const std::array<T, 4> m_background;

		/* Calculates the sum of 2 log probabilities */
		T sum_log_probs(T a, T b);

		/* Calculates the background taxon distribution in m_data */
		std::array<T, 4> calculate_noise(int sample_size = 100);
};

template <typename T>
GibbsSampler<T>::GibbsSampler(const Data& data) 
	: m_data { data },
	  m_background { calculate_noise() }
{
}

template <typename T>
T GibbsSampler<T>::sum_log_probs(T a, T b)
{
	return a > b ?
		a + std::log1p(std::exp(b-a)) :
		b + std::log1p(std::exp(a-b));
}

template <typename T>
std::array<T, 4> GibbsSampler<T>::calculate_noise(int sample_size)
{
    // sample 100 positions with replacement
    auto [num_sequences, sequence_length] { m_data.size() };
    std::array<T, 4> result; 

    char tmp;
    int samples_per_seq { (sample_size + num_sequences - 1) / num_sequences };
    const std::vector<Sequence> seqs { m_data.sequences() };
    for (int i {}; i < num_sequences; ++i) {
        for (int j {}; j < samples_per_seq; ++j) {
            auto idx = utility::rand_indices(sequence_length)[0]; 
            tmp = seqs[i].m_sequence[idx];    
            ++result[utility::encode(tmp)];
        }
	}

    // normalize to get probability distribution
    int total_samples { samples_per_seq * num_sequences };
    std::transform(begin(result), end(result), begin(result), [&total_samples](T x) {
        return x / total_samples;
    });

    return result;
}

template <typename T>
int GibbsSampler<T>::num_correct(std::vector<int>& positions, int k)
{
	// what if have multiple motifs and looking for 1
	// score based on all
	// return max score
	std::unordered_map<int, int> results {};
	std::vector<Sequence> seqs { m_data.sequences() };
	for (int i {}; i < positions.size(); ++i) {
		for (int j {}; j < seqs[i].m_motifs.size(); ++j) {
			auto& motif { seqs[i].m_motifs[j] };
			if (std::abs(static_cast<int>(positions[i] - motif.m_startingIndex)) < k) {
				++results[motif.m_motifId];
				break;
			}
		}
	}

	auto elm_it = std::max_element(begin(results), end(results), 
		[](const std::pair<int, int>& a, const std::pair<int, int>& b) {
			return a.second < b.second;
		});
	int result { results.empty() ? 0 : elm_it->second };
	return result;
}

template <typename T>
std::string GibbsSampler<T>::consensus(std::vector<T>& pwm, int k)
{
	std::string result {};

	for (int i {}; i < k; ++i) {
		auto index {
			std::distance(
				begin(pwm)+4*i,
				std::max_element(begin(pwm)+4*i, begin(pwm)+4*(i+1))
			)
		};
		assert(index < 4); // decode requirement
		result.push_back(utility::decode(index));
	}

	return result;
}

template <typename T>
std::vector<int> GibbsSampler<T>::init_positions(int width)
{
	auto [num_sequences, sequence_length] { m_data.size() };
	std::vector<int> result(num_sequences);

	std::generate(begin(result), end(result), [&]() {
		return utility::rand_indices(sequence_length, width)[0]; 
	});

    return result; 
}

template <typename T>
std::vector<T> GibbsSampler<T>::init_pwm(std::vector<int>& positions, int k, T pseudocount) 
{
	/*
	In PWM, each nucleotide increases weight by 1/(k+4*pseudo)
	and we start with pseduo/(k+4*pseudo)
	*/
	T normalized_default { pseudocount / (k + 4*pseudocount) };
    std::vector<T> pwm(4*k, normalized_default);

    assert(m_data.sequences().size() == positions.size());
    for (int i {}; i < positions.size(); ++i) {
		update_counts(pwm, i, positions[i], k, pseudocount);
	}

    return pwm;
}

template <typename T>
std::vector<T> GibbsSampler<T>::update_pwm(std::vector<T>& pwm, 
	std::vector<int>& positions, int k, T pseudocount, int old_withheld,
	int new_withheld) 
{
	std::vector<T> old_pwm { pwm };

	update_counts(pwm, old_withheld, positions[old_withheld], k, pseudocount);
	update_counts(pwm, new_withheld, positions[new_withheld], k, pseudocount, false);

	return old_pwm;
}

template <typename T>
void GibbsSampler<T>::update_counts(std::vector<T>& pwm, int seq_index, int start_pos, 
	int k, T pseudocount, bool increment) 
{
	T delta = (increment ? 1 : -1) * 1 / (k + 4 * pseudocount) ;

	const auto& seq { m_data.sequences()[seq_index].m_sequence };
	for (int i {}; i < k; ++i) {
		int idx { 4*i + utility::encode(seq[i+start_pos]) };
		pwm[idx] += delta;    
	}
}

// O(seq_len * k)
template <typename T>
std::vector<T> GibbsSampler<T>::score(std::vector<T>& pwm, int k, int withheld) 
{
	// TODO: if very slow, add thresholding, where only sample if score > some value
    auto [num_sequences, sequence_length] { m_data.size() };
	std::vector<T> score(sequence_length-k);

	const auto& seq { m_data.sequences()[withheld].m_sequence };
	for (int i {}; i < sequence_length-k; ++i) {  // iterate over possible starting positions
		T tmp {};
		for (int j { }; j < k; ++j) {  // iterates over single kmer
			int nucleotide_encoding { utility::encode(seq[i+j]) };
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
int GibbsSampler<T>::sample(std::vector<T> scores) 
{
    std::random_device rand_device {};
	std::mt19937 num_gen { rand_device() };

	std::discrete_distribution<int> discrete_distr(begin(scores), end(scores));
	return discrete_distr(num_gen);
}
