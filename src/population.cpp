// -*- mode: c++; coding: utf-8 -*-
/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>
#include <sfmt.hpp>

#include <iostream>

namespace pbt {

std::ostream& Population::write(std::ostream& ost) const {HERE;
    for (const auto& x: individuals_) {
        ost << x;
    }
    return ost;
}

//! shortcut Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

void Population::test() {HERE;
    Population x;
    std::cout << x << std::endl;
}

} // namespace pbt
