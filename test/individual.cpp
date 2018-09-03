#include "individual.hpp"

#include <iostream>

int main() {
    std::cout << "sizeof(Individual): " << sizeof(pbt::Individual) << "\n";
    pbt::Individual::set_default_values();
    pbt::Individual x;
    std::cout << x << std::endl;
    return 0;
}
