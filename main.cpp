#include "iostream"
#include "string"
#include "vector"

#include "data.hpp"
#include "serial.hpp"

int main() {    
    std::vector<std::size_t> motif_lengths { 12 };       
	std::size_t num_sequences { 50 };
	std::size_t sequence_length { 1000 };

    Data data { motif_lengths, num_sequences, sequence_length };
	std::cout << data << std::endl;

	Serial<float> serial { data };
	std::vector<std::size_t> ending_positions { serial.find_motifs(motif_lengths[0], 1) };

	auto str { std::accumulate(
		std::next(begin(ending_positions)), 
		end(ending_positions),
		std::to_string(*begin(ending_positions)),
		[](const std::string& a, std::size_t b) {
			return a + " " + std::to_string(b);
		}
	)};
	std::cout << str << std::endl;

    return 0;
}
