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

struct Motif
{
	std::string m_motif;
	std::string m_baseMotif;
	std::size_t m_startingIndex;
	std::size_t m_motifId;
};

struct Sequence
{
	std::string m_sequence;
	std::vector<Motif> m_motifs;
};

class Data
{
    public: 
        Data(
            const std::vector<std::size_t>& motif_lengths, 
            std::size_t num_sequences = 10,
            std::size_t sequence_length = 1'000
        );
        const std::vector<Sequence>& sequences() const;
        const std::pair<std::size_t, std::size_t> size() const;
        
        friend std::ostream& operator<<(std::ostream& os, const Data& obj);

    private:
        const std::size_t m_numSequences; 
        const std::size_t m_sequenceLength;
        const std::vector<std::size_t> m_motifLengths;

		/* simulated consensus motifs before random obfuscation */
        const std::vector<std::string> m_motifs;

        std::vector<Sequence> m_sequences;    

        /* Returns N motifs with lengths corresponding to motif_lengths */
        std::vector<std::string> generate_motifs(); 

        /* Generates a Sequence with a set of motifs with lengths specified in m_motifLengths */
        Sequence generate_sequence();
};

