#include "population.hpp"

#include <iostream>

int main() {
    pbt::Population pop(1000u);
    pop.run(10u, 0.01);
    return 0;
}
