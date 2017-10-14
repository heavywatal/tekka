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
        if (!mother.is_matured(time_)) continue;
        const Individual& father = *wtl::choice(males_.begin(), males_.end(), wtl::sfmt());
        // TODO: weighted sampling
        if (!father.is_matured(time_)) continue;
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

void Population::survive() {HERE;
   std::vector<Individual> survivors;
   survivors.reserve(males_.size());
   std::copy_if(males_.begin(), males_.end(),
                std::back_inserter(survivors), [this](const Individual& x) {
                    return x.has_survived(time_);
                });
   males_.swap(survivors);
   survivors.clear();
   std::copy_if(females_.begin(), females_.end(),
                std::back_inserter(survivors), [this](const Individual& x) {
                    return x.has_survived(time_);
                });
   females_.swap(survivors);
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
    x.survive();
    std::cout << x << std::endl;
}

} // namespace pbt
