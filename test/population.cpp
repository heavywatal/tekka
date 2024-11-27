#include "population.hpp"
#include "parameters.hpp"

int main() {
    pbf::Parameters params;
    pbf::Population pop(params);
    pop.run();
    return 0;
}
