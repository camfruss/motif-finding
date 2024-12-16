#pragma once

#include "algorithm"
#include "iostream"
#include "random"
#include "set"
#include "string"
#include "unordered_map"
#include "utility"
#include "vector"

#include "utility.hpp"

class Data
{
    public: 
        Data(
            const std::vector<std::size_t>& motif_lengths, 
            std::size_t num_sequences = 10,
            std::size_t sequence_length = 1'000
        );
        const std::vector<std::string>& sequences() const;
        const std::pair<std::size_t, std::size_t> size() const;
        
        friend std::ostream& operator<<(std::ostream& os, const Data& obj);

    private:
        const std::size_t m_numSequences; 
        const std::size_t m_sequenceLength;
        const std::vector<std::size_t> m_motifLengths;
        const std::vector<std::string> m_motifs;
        std::vector<std::string> m_sequences;    

        // Returns N motifs with lengths corresponding to motif_lengths
        std::vector<std::string> generate_motifs(); 

        // Generates a sequence with a set of motifs with lengths specified in m_motifLengths
        // returns: a string consisting of the alphabet {A, C, T, G}
        std::string generate_sequence();
};
