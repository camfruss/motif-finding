#include "cstdio"
#include "vector"

#include "data.hpp"

int main() {    
    std::vector<std::size_t> motif_lengths { 10 };       

    Data data { motif_lengths };
    std::cout << data << std::endl;

    return 0;
}
