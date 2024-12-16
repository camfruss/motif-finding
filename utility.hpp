#pragma once

#include "random"
#include "unordered_map"
#include "vector"

namespace utility {
    namespace {
        std::random_device rand_device {};
        std::mt19937 num_gen { rand_device() };
        std::discrete_distribution<> discrete_distr { 25, 25, 25, 25 };

        std::unordered_map<std::size_t, char> nucleotide_map =  { {0, 'A'}, {1, 'C'}, {2, 'T'}, {3, 'G'} };
        std::unordered_map<char, std::size_t> encoding_map =  { {'A', 0}, {'C', 1}, {'T', 2}, {'G', 3} };
    }
   
    // c must be in {A, C, T, G} 
    inline std::size_t encode(char c)
    {
        return encoding_map[c];
    }

    // i must be in {0, 1, 2, 3}
    inline char decode(std::size_t i)
    {
        return nucleotide_map[i];
    }
   
    inline std::size_t rand_nucleotide()
    {
        return nucleotide_map[discrete_distr(num_gen)];
    }

    /*
    * max: selected random index from tthe uniform distribution [0, max)
    * count: number of indices to return
    * end_buffer: pseudo-length of index, such that if count > 1, they are 
    * guaranteed to be end_buffer indices away from each over 
    */ 
    std::vector<std::size_t> rand_indices(
        std::size_t max, 
        std::size_t end_buffer = 1, 
        std::size_t count = 1
    );
}

