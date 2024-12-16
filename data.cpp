#include "algorithm"
#include "iostream"
#include "random"
#include "set"
#include "string"
#include "unordered_map"
#include "utility"
#include "vector"

#include "data.hpp"
#include "utility.hpp"

Data::Data(
    const std::vector<std::size_t>& motif_lengths, 
    std::size_t num_sequences,
    std::size_t sequence_length
) : m_numSequences { num_sequences },
    m_sequenceLength { sequence_length }, 
    m_motifLengths { motif_lengths }, 
    m_motifs { generate_motifs() }
{
    m_sequences.resize(m_numSequences);
    std::generate(begin(m_sequences), end(m_sequences), [this]() {
        return this->generate_sequence();
    });
}

const std::vector<std::string>& Data::sequences() const
{
    return m_sequences;
}

const std::pair<std::size_t, std::size_t> Data::size() const
{
    return { m_numSequences, m_sequenceLength };
}

std::ostream& operator<<(std::ostream& os, const Data& obj)
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

std::vector<std::string> Data::generate_motifs() 
{
    std::vector<std::string> result {};
    for (const auto& size : m_motifLengths) {
        std::string tmp {};
        for (std::size_t i {}; i < size; ++i) {
            tmp.push_back(utility::rand_nucleotide());
        }
        result.push_back(tmp);
    }
    return result;
}

std::string Data::generate_sequence() 
{
    std::string result {};
    result.reserve(m_sequenceLength);
    
    // create initial values
    for (std::size_t i {}; i < m_sequenceLength; ++i) {
        result.push_back(utility::rand_nucleotide());
    }

    // insert motifs
    std::size_t end_buffer { *std::max_element(begin(m_motifLengths), end(m_motifLengths)) };
    auto indices { utility::rand_indices(m_sequenceLength, end_buffer, m_motifs.size()) }; 
    for (std::size_t i {}; i < m_motifs.size(); ++i) {
        result.replace(indices[i], m_motifLengths[i], m_motifs[i]);
    }

    return result; 
}

