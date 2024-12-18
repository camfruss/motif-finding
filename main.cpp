#include "iostream"
#include "string"
#include "vector"

#include "data.hpp"
#include "gibbs_sampler.hpp"
#include "serial.hpp"

int main() {    
    std::vector<int> motif_lengths { 12 };       
	int num_sequences { 10 };
	int sequence_length { 500 };

    Data data { motif_lengths, num_sequences, sequence_length };
	std::cout << data << std::endl;

	Serial<float> serial { data };
	Result result { serial.find_motifs(motif_lengths[0], 0.1) };
	std::vector<int> ending_positions { result.positions };

	auto str { std::accumulate(
		std::next(begin(ending_positions)), 
		end(ending_positions),
		std::to_string(*begin(ending_positions)),
		[](const std::string& a, int b) {
			return a + " " + std::to_string(b);
		}
	)};
	std::cout << "num correct: " << result.num_correct << "\n";
	std::cout << str << std::endl;

    return 0;
}
