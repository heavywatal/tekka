#include "individual.hpp"

#include <iostream>

int main() {
    std::cout << "sizeof(Individual): " << sizeof(pbf::Individual) << "\n";
    pbf::Individual x(false);
    std::cout << x << std::endl;
    return 0;
}
