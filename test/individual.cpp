#include "individual.hpp"

#include <iostream>

int main() {
    pbt::Individual::set_default_values();
    pbt::Individual x;
    std::cout << x << std::endl;
    return 0;
}
