#include "population.hpp"
#include "individual.hpp"

#include <iostream>

int main() {
    pbt::Individual::set_default_values();
    pbt::Population pop(42u);
    pop.run(8u, 8u);
    return 0;
}