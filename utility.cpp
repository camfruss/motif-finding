#include "random"
#include "set"
#include "unordered_map"
#include "vector"

#include "utility.hpp"

#include "iostream"

std::vector<int> utility::rand_indices(int max, int width, int count) 
{
    std::random_device rand_device {};
    std::mt19937 num_gen { rand_device() };
    
    std::vector<int> result {};
    std::uniform_int_distribution<int> dist(0, max-width);
    std::set<int> positions {};

    for (int i {}; i < count; ++i) {
        bool valid;
        std::vector<int> tmp {};
        int pos;

        do {  // WARNING: may cause an infinite loop if too many motifs 


            valid = true;
            pos = dist(num_gen);
            for (int i { pos }; i < pos + width; ++i) {
                if (positions.count(i) || i == max) { 
                    valid = false; 
                    tmp.clear();
					break;
                } else {
					tmp.push_back(i);
				}
            } 
        } while (!valid);
        
        result.push_back(pos);
        positions.insert(begin(tmp), end(tmp));                    
    }
    
    return result;
}

