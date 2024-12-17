#include "cstdio"
#include "vector"

#include "data.hpp"
#include "serial.hpp"

int main() {    
    std::vector<std::size_t> motif_lengths { 10 };       
	std::size_t num_sequences { 10 };
	std::size_t sequence_length { 100 };

    Data data { motif_lengths, num_sequences, sequence_length };
	Serial<float> serial { data };
	
	std::cout << "Beginning motif search" << std::endl;
	std::vector<std::size_t> starting_positions { serial.find_motifs(10, 1) };
	

	
	std::cout << data << std::endl;

    return 0;
}
