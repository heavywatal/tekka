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

Population::Population(const size_t initial_size) {HERE;
    const size_t half = initial_size / 2UL;
    males_.resize(half);
    females_.resize(half);
}

void Population::reproduce() {HERE;
    ++time_;
    std::vector<Individual> boys;
    std::vector<Individual> girls;
    constexpr size_t clutch_size = 2;
    std::poisson_distribution<unsigned int> poisson_eggs(clutch_size);
    for (const auto& mother: females_) {
        const Individual& father = *wtl::choice(males_.begin(), males_.end(), wtl::sfmt());
        const unsigned int num_eggs = poisson_eggs(wtl::sfmt());
        for (unsigned int i=0; i<num_eggs; ++i) {
            if (wtl::sfmt().canonical() < 0.5) {
                boys.emplace_back(father, mother, time_);
            } else {
                girls.emplace_back(father, mother, time_);
            }
        }
    }
    std::copy(boys.begin(), boys.end(), std::back_inserter(males_));
    std::copy(girls.begin(), girls.end(), std::back_inserter(females_));
}

std::ostream& Population::write(std::ostream& ost) const {HERE;
    auto write_impl = [&ost](const decltype(males_)& v) {
        for (const auto& x: v) {
            ost << x << "\n";
        }
    };
    write_impl(males_);
    write_impl(females_);
    return ost;
}

//! shortcut Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

void Population::test() {HERE;
    Population x(6);
    std::cout << x << std::endl;
    x.reproduce();
    std::cout << x << std::endl;
}

} // namespace pbt
