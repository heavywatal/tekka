#include "population.hpp"

#include <random>

int main() {
    pbf::Population pop(200u, std::random_device{}());
    pop.run();
    return 0;
}
