#include "individual.hpp"

#include <iostream>

int main() {
    std::cout << "sizeof(Individual): " << sizeof(pbt::Individual) << "\n";
    pbt::Individual x;
    std::cout << x << std::endl;
    return 0;
}
