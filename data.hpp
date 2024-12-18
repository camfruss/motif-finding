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
    /* Possibly obfuscated motif due to simulated read errors */
	std::string m_motif;

    /* True embedded motif before obfuscation */
	std::string m_baseMotif;

    /* Location where the motif starts */
	int m_startingIndex;

    /* Unique id of base motif */
	int m_motifId;
};

struct Sequence
{
    /* Full ACTG representation */
	std::string m_sequence;

    /* Motifs within the sequence */
	std::vector<Motif> m_motifs;
};

class Data
{
    public: 
        /* Initializes a sequence dataset with embedded motifs
         * motif_lengths : vector containing the length of motifs to embed
         */
        Data(
            const std::vector<int>& motif_lengths, 
            int num_sequences = 10,
            int sequence_length = 1'000
        );

        /* Returns all created Sequences */
        const std::vector<Sequence>& sequences() const;

        /* Returns (num_sequences, sequence_length) */
        const std::pair<int, int> size() const;
        
        /* Allows for pretty printing Data */
        friend std::ostream& operator<<(std::ostream& os, const Data& obj);

    private:
        const int m_numSequences; 
        const int m_sequenceLength;
        const std::vector<int> m_motifLengths;

		/* simulated consensus motifs before random obfuscation */
        const std::vector<std::string> m_motifs;

        std::vector<Sequence> m_sequences;    

        /* Returns N motifs with lengths corresponding to motif_lengths */
        std::vector<std::string> generate_motifs(); 

        /* Generates a Sequence with a set of motifs with lengths specified in m_motifLengths */
        Sequence generate_sequence();
};

