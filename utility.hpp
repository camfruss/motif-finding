#pragma once

#include "random"
#include "unordered_map"
#include "vector"

namespace utility {
    namespace {
        std::random_device rand_device {};
		std::mt19937 num_gen { rand_device() };
        std::discrete_distribution<> discrete_distr { 25, 25, 25, 25 };

        std::unordered_map<int, char> nucleotide_map =  { {0, 'A'}, {1, 'C'}, {2, 'T'}, {3, 'G'} };
        std::unordered_map<char, int> encoding_map =  { {'A', 0}, {'C', 1}, {'T', 2}, {'G', 3} };
    }
   
	/* c must be in {A, C, T, G} */
    inline int encode(char c)
    {
        return encoding_map[c];
    }

    /* i must be in {0, 1, 2, 3} */
    inline char decode(int i)
    {
        return nucleotide_map[i];
    }
   
    /* Returns a random char in {A, C, T, G}*/
    inline char rand_nucleotide()
    {
        return nucleotide_map[discrete_distr(num_gen)];
    }

    /* Returns count random indices within the range [0, max-width]
     * max: selected random index from the uniform distribution [0, max)
     * count: number of indices to return
     * width: pseudo-length of index, such that if count > 1, they are 
     * guaranteed to be width indices away from each over to prevent overlap
     */ 
    std::vector<int> rand_indices(int max, int width = 1, int count = 1);
}
