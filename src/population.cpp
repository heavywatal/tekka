/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>

#include <iostream>

namespace pbt {

Population::Population(const size_t initial_size) {HERE;
    const size_t half = initial_size / 2UL;
    males_.resize(half);
    females_.resize(half);
}

void Population::run(const uint_fast32_t years) {HERE;
    urbg_t urbg(std::random_device{}());
    for (year_ = 3u; year_ < years; ++year_) {
        reproduce(urbg);
        survive(urbg);
        survive(urbg);
        survive(urbg);
        survive(urbg);
        migrate(urbg);
        std::cerr << year_ << ": " << males_.size() + females_.size() << std::endl;
    }
}

void Population::reproduce(urbg_t& urbg) {
    std::vector<Individual> boys;
    std::vector<Individual> girls;
    std::vector<std::vector<const Individual*>> adult_males(Individual::num_locations());
    for (const auto& x: males_) {
        adult_males[x.location()].push_back(&x);
    }
    for (const auto& mother: females_) {
        if (!mother.is_in_breeding_place()) continue;
        const auto& potential_fathers = adult_males[mother.location()];
        if (potential_fathers.empty()) continue;
        const unsigned int num_eggs = mother.clutch_size(urbg);
        const Individual* father = *wtl::choice(potential_fathers.begin(), potential_fathers.end(), urbg);
        // TODO: multiple fathers
        for (unsigned int i=0; i<num_eggs; ++i) {
            if (urbg.canonical() < 0.5) {
                boys.emplace_back(*father, mother, year_);
            } else {
                girls.emplace_back(*father, mother, year_);
            }
        }
    }
    std::copy(boys.begin(), boys.end(), std::back_inserter(males_));
    std::copy(girls.begin(), girls.end(), std::back_inserter(females_));
}

void Population::survive(urbg_t& urbg) {
    auto impl = [&urbg, year = year_](const decltype(males_)& v) {
        decltype(males_) survivors;
        survivors.reserve(v.size());
        std::copy_if(v.begin(), v.end(), std::back_inserter(survivors),
                     [&urbg, &year](const Individual& x) {return x.has_survived(year, urbg);});
        return survivors;
    };
    males_ = impl(males_);
    females_ = impl(females_);
}

void Population::migrate(urbg_t& urbg) {
    for (auto& x: males_) {x.migrate(year_, urbg);}
    for (auto& x: females_) {x.migrate(year_, urbg);}
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
    Population x(12);
    std::cout << x << std::endl;
    x.run(15);
    std::cout << x << std::endl;
}

} // namespace pbt
