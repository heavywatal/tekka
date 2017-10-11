// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>
#include <sfmt.hpp>

#include <iostream>

namespace pbt {

std::ostream& Individual::write(std::ostream& ost) const {HERE;
    return ost << id_ << ":" << is_male_;
}

//! shortcut Individual::write(ost)
std::ostream& operator<<(std::ostream& ost, const Individual& x) {
    return x.write(ost);
}

void Individual::test() {HERE;
    Individual x;
    std::cout << x << std::endl;
}

} // namespace pbt
