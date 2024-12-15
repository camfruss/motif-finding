#include "algorithm"
#include "iostream"
#include "random"
#include "set"
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
        ) : m_numSequences { num_sequences }, 
            m_sequenceLength { sequence_length }, 
            m_motifs { generate_motifs(motif_lengths) }
        {
            m_sequences.resize(m_numSequences);
            std::generate(begin(m_sequences), end(m_sequences), [this]() {
                return this->generate_sequence();
            });
        }

        const std::vector<std::string> sequences() 
        {
            return m_sequences;
        }

        friend std::ostream& operator<<(std::ostream& os, const Data& obj)
        {
            os << "Motifs:\n";
            for (const auto& motif : obj.m_motifs) {
                os << motif << "\n";
            }    

            std::size_t count {};
            for (const auto& seq : obj.m_sequences) {
                os << "\n> sequence " << count+1 << "\n";  // std::format supported in g++ 13.1
                for (std::size_t i {}; i < seq.length(); i+=80) {
                    os << seq.substr(i, 80) << "\n";
                }
                ++count;
            }
            return os;
        }

    private:
        std::random_device ms_randDevice {};
        std::mt19937 ms_numGen { ms_randDevice() };
        std::discrete_distribution<> ms_discreteDistr { 25, 25, 25, 25 };
        std::unordered_map<int, char> ms_nucleotideMap =  { {0, 'A'}, {1, 'C'}, {2, 'T'}, {3, 'G'} };

        const std::size_t m_numSequences; 
        const std::size_t m_sequenceLength;
        const std::vector<std::string> m_motifs;
        std::vector<std::string> m_sequences;    

        inline char generate_nucleotide() 
        {
            return ms_nucleotideMap[ms_discreteDistr(ms_numGen)];
        }

        std::vector<std::string> generate_motifs(const std::vector<std::size_t>& motif_lengths) 
        {
            std::vector<std::string> result {};
            for (const auto& size : motif_lengths) {
                std::string tmp {};
                for (std::size_t i {}; i < size; ++i) {
                    tmp.push_back(generate_nucleotide());
                }
                result.push_back(tmp);
            }
            return result;
        }

        // Generates a sequence with a set of motifs with lengths specified in m_motifLengths
        // returns: a string consisting of the alphabet {A, C, T, G}
        std::string generate_sequence() 
        {
            std::string result {};
            result.reserve(m_sequenceLength);
            
            // create initial values
            for (std::size_t i {}; i < m_sequenceLength; ++i) {
                result.push_back(generate_nucleotide());
            }

            // insert motifs
            // pick m_motifs.size() random spots such that at least N away from each other
            std::uniform_int_distribution<std::size_t> dist(0, m_sequenceLength);
            std::set<std::size_t> positions {};
            for (const auto& motif : m_motifs) {
            
                bool valid { true };
                std::vector<std::size_t> tmp {};
                std::size_t pos;
                do {  // WARNING: may cause an infinite loop if too many motifs 
                    pos = dist(ms_numGen);
                    for (std::size_t i { pos }; i < pos + motif.size(); ++i) {
                        if (positions.count(i) || i == m_sequenceLength) { 
                            valid = false; 
                            tmp.clear();
                        }
                        tmp.push_back(i);
                    } 
                } while (!valid);
                result.replace(pos, motif.size(), motif);
                positions.insert(begin(tmp), end(tmp));                    

            }

            return result; 
        }
};
