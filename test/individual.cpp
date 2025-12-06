#include "individual.hpp"

#include <fmt/base.h>

#include <type_traits>

static_assert(std::is_nothrow_default_constructible_v<pbf::Individual>);
static_assert(!std::is_copy_constructible_v<pbf::Individual>);
static_assert(!std::is_move_constructible_v<pbf::Individual>);
static_assert(std::is_nothrow_destructible_v<pbf::Individual>);

int main() {
    fmt::println("sizeof(Individual): {}", sizeof(pbf::Individual));
    pbf::Individual x{};
    fmt::println("{}", x.str());
    return 0;
}
