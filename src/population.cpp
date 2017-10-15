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
   auto impl = [t=time_](const decltype(males_)& v) {
       decltype(males_) survivors;
       survivors.reserve(v.size());
       std::copy_if(v.begin(), v.end(), std::back_inserter(survivors),
                    [&t](const Individual& x) {return x.has_survived(t);});
       return survivors;
   };
   males_ = impl(males_);
   females_ = impl(females_);
}

std::ostream& Population::write(std::ostream& ost) const {HERE;
    auto impl = [&ost](const decltype(males_)& v) {
        for (const auto& x: v) {
            ost << x << "\n";
        }
    };
    impl(males_);
    impl(females_);
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
