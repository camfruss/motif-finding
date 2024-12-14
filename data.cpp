#include "random"
#include "string"
#include "unordered_map"
#include "vector"

class Data
{
    public: 
        Data(
            const std::vector<std::size_t>& motif_lengths, 
            std::size_t num_sequences = 10,
            std::size_t sequence_length = 1'000
        ) : m_motifLengths { motif_lengths }, 
            m_numSequences { num_sequences }, 
            m_sequenceLength { sequence_length }, 
            m_sequences(num_sequences, generate_sequence()) 
        {
        }

    private:
        const std::vector<std::size_t> m_motifLengths;
        std::size_t m_numSequences; 
        std::size_t m_sequenceLength;
        std::vector<std::string> m_sequences;    

        // Generates a sequence with a set of motifs with lengths specified in m_motifLengths
        // returns: a string consisting of the alphabet {A, C, T, G}
        std::string generate_sequence()
        {
            std::random_device rd {};
            std::mt19937 gen { rd() };
            std::discrete_distribution<> dd {25, 25, 25, 25};

            std::unordered_map<int, char> map =  { {0, 'A'}, {1, 'C'}, {2, 'T'}, {3, 'G'} };

            std::string result {};
            result.reserve(m_sequenceLength);
            for (std::size_t i {}; i < m_sequenceLength; ++i) {
                result.push_back(map[dd(gen)]);
            }
            return result; 
        }

};

