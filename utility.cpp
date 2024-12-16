#include "random"
#include "set"
#include "unordered_map"
#include "vector"

#include "utility.hpp"

std::vector<std::size_t> utility::rand_indices(
    std::size_t max, 
    std::size_t end_buffer, 
    std::size_t count
) {
    std::vector<std::size_t> result {};
    std::uniform_int_distribution<std::size_t> dist(0, max);
    std::set<std::size_t> positions {};

    for (std::size_t i {}; i < count; ++i) {
        bool valid { true };
        std::vector<std::size_t> tmp {};
        std::size_t pos;

        do {  // WARNING: may cause an infinite loop if too many motifs 
            pos = dist(num_gen);
            for (std::size_t i { pos }; i < pos + end_buffer; ++i) {
                if (positions.count(i) || i == max) { 
                    valid = false; 
                    tmp.clear();
                }
                tmp.push_back(i);
            } 
        } while (!valid);
        
        result.push_back(pos);
        positions.insert(begin(tmp), end(tmp));                    
    }
    
    return result;
}

