#include "population.hpp"
#include "individual.hpp"

#include <iostream>

int main() {
    pbt::Individual::set_default_values();
    pbt::Population pop(1000u);
    pop.run(10u, 0.01);
    return 0;
}