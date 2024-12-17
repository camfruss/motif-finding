#include "algorithm"
#include "format"
#include "iostream"
#include "random"
#include "ranges"
#include "set"
#include "string"
#include "unordered_map"
#include "utility"
#include "vector"

#include "data.hpp"
#include "utility.hpp"

#include "iostream"

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

const std::vector<Sequence>& Data::sequences() const
{
    return m_sequences;
}

const std::pair<std::size_t, std::size_t> Data::size() const
{
    return { m_numSequences, m_sequenceLength };
}

std::ostream& operator<<(std::ostream& os, const Data& obj)
{
    os << "CONSENSUS MOTIFS:\n";
	const auto& motifs { obj.m_motifs };
    for (std::size_t i {}; i < motifs.size(); ++i) {
        os << std::format("{:02} > {}\n", i+1, motifs[i]);
    }    

    std::size_t count {};
    for (const auto& seq : obj.m_sequences) {
	
		// print out indices where motifs were inserted for each sequence
		auto indices = seq.m_motifs | std::views::transform([](const Motif& m) {
        	return std::to_string(m.m_startingIndex);
		});
		std::string motif_indices {
			std::accumulate(
				std::next(begin(indices)), end(indices), *begin(indices),
				[](const std::string& a, const std::string& b) {
		        	return a + ", " + b;
				})
		};
		os << std::format("> sequence {} | motif indices: {}\n", count+1, motif_indices); // std::format supported in g++ 13.1

		for (std::size_t i {}; i < seq.m_sequence.length(); i+=80) {
            os << seq.m_sequence.substr(i, 100) << "\n";
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

Sequence Data::generate_sequence() 
{
   	std::vector<Motif> motifs {};

    // create initial valuesi
	std::string sequence {};
	sequence.reserve(m_sequenceLength);
    for (std::size_t i {}; i < m_sequenceLength; ++i) {
        sequence.push_back(utility::rand_nucleotide());
    }

    // insert motifs
    std::size_t end_buffer { *std::max_element(begin(m_motifLengths), end(m_motifLengths)) };
    auto indices { utility::rand_indices(m_sequenceLength, end_buffer, m_motifs.size()) }; 
    for (std::size_t i {}; i < m_motifs.size(); ++i) {
        sequence.replace(indices[i], m_motifLengths[i], m_motifs[i]);
		motifs.push_back({
			m_motifs[i],
			m_motifs[i], // TODO: obfuscate w/ some SMALL probability
			indices[i],
			i
		});
    }

	Sequence result { sequence, motifs };
    return result; 
}

