#include "population.hpp"

#include <iostream>
#include <random>

int main() {
    pbt::Population pop(1000u, std::random_device{}());
    pop.run(10u, 0.01);
    return 0;
}
