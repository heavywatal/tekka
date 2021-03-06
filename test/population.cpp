#include "population.hpp"

#include <iostream>
#include <random>

int main() {
    pbf::Population pop(1000u, std::random_device{}());
    pop.run(10u);
    return 0;
}
