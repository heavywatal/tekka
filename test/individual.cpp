#include "individual.hpp"

#include <iostream>
#include <type_traits>

static_assert(std::is_nothrow_default_constructible_v<pbf::Individual>);
static_assert(!std::is_copy_constructible_v<pbf::Individual>);
static_assert(!std::is_move_constructible_v<pbf::Individual>);
static_assert(std::is_nothrow_destructible_v<pbf::Individual>);

int main() {
    std::cout << "sizeof(Individual): " << sizeof(pbf::Individual) << "\n";
    pbf::Individual x{};
    std::cout << x << std::endl;
    return 0;
}
