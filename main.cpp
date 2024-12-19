#include "data.hpp"
#include "gibbs_sampler.hpp"
#include "iostream"
#include "serial.hpp"
#include "string"
#include "vector"

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <num_motifs> <motif_lengths> <num_sequences> "
                     "<sequence_length>\n";
        return 1;
    }

    int num_m = std::stoi(argv[1]);
    int m_len = std::stoi(argv[2]);
    int num_s = std::stoi(argv[3]);
    int s_len = std::stoi(argv[4]);
    std::vector<int> motif_lengths(num_m, m_len);
    int num_sequences{num_s};
    int sequence_length{s_len};

    Data data{motif_lengths, num_sequences, sequence_length};
    std::cout << data << std::endl;

    Serial<float> serial{data};
    Result result{serial.find_motifs(motif_lengths[0], 0.1)};
    std::vector<int> ending_positions{result.positions};

    auto str{std::accumulate(std::next(begin(ending_positions)),
                             end(ending_positions),
                             std::to_string(*begin(ending_positions)),
                             [](const std::string& a, int b) {
                                 return a + " " + std::to_string(b);
                             })};
    std::cout << "num correct: " << result.num_correct << "\n";
    std::cout << str << std::endl;

    return 0;
}
