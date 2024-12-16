/*
calculate background noise
choose random starting positions 
take out 1 sequence 
create position weight matrix 
score every kmer in selected sequence
sample a new starting position weighted by scores
*/
#include "algorithm"
#include "array"
#include "cassert"
#include "string"
#include "vector"

#include "serial.hpp"
#include "data.hpp"
#include "utility.hpp"

template <typename T>
Serial<T>::Serial(const Data& data) 
    : m_background { calculate_noise() }, 
      m_data { data }
{
}

template <typename T>
std::vector<std::size_t> Serial<T>::find_motifs(std::size_t k)
{
    auto [num_sequences, sequence_length] { m_data.size() };

    std::vector<std::size_t> positions { init_positions(k) }; 
    auto pwm { init_pwm(k, positions) };

    // TODO: experiment with different convergence criterion...
    // scoring metric -> sum of scores of complete iteration
    // fixed number of iterations / max iters
    // when motifs remain unchanged / only a few change
    // scoring metric + max_iterations
    // can add correctness plots as well
    int max_iters { 1'000 }; // TODO: this should be a parameter
    auto has_converged = [&]() {
        static std::size_t iter_count {}; 
        return iter_count++ < max_iters;
    };

    std::size_t withheld_index { 0 }; 
    std::vector<std::size_t> scores(sequence_length);  // technically seq_len - k
    std::size_t sampled_index { 0 };
    while (!has_converged()) {
        update_pwm(pwm, withheld_index);
        scores = score(pwm, withheld_index);  // TODO
        sampled_index = sample(scores);  // TODO
        positions[withheld_index] = sampled_index; 
        withheld_index = (withheld_index + 1) % num_sequences; 
    }

    return positions;
}

template <typename T>
std::array<T, 4> Serial<T>::calculate_noise(int sample_size)
{
    // sample 100 positions with replacement
    auto [num_sequences, sequence_length] { m_data.size() };
    std::array<T, 4> result; 

    char tmp;
    int samples_per_seq { (sample_size + num_sequences - 1) / num_sequences };
    const std::vector<std::string> seqs { m_data.sequences() };
    for (std::size_t i {}; i < num_sequences; ++i) {
        for (std::size_t j {}; j < samples_per_seq; ++j) {
            auto idx = utility::rand_indices(sequence_length)[0]; 
            tmp = seqs[idx];    
            ++result[utility::encode(tmp)];
        }
    }

    // normalize to get probability distribution
    int total_samples { samples_per_seq * num_sequences };
    auto normalize = std::transform(begin(result), end(result), begin(result), [&total_samples](float x) {
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
        return utility::rand_indices(sequence_length, end_buffer);
    });
    
    return result; 
}

template <typename T>
std::vector<T> Serial<T>::init_pwm(int k, std::vector<std::size_t>& positions)
{
    std::vector<T> result(4*k, 0);

    assert(m_data.sequences().size() == positions.size());
    for (std::size_t i {}; i < positions.size(); ++i) {
        auto starting_pos { positions[i] };
        auto seq { m_data.sequences()[i] };
        for (std::size_t j { starting_pos }; j < starting_pos + k; ++j) {
            std::size_t idx { 4*i + utility::encode(seq[j]) };
            ++result[idx];    
        }
    }

    return result;
}

template <typename T>
void Serial<T>::update_pwm()
{
}


