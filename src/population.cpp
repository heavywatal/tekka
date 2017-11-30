/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/random.hpp>

#include <iostream>

namespace pbt {

Population::Population(const size_t initial_size) {HERE;
    engine_.seed(wtl::random_device_64{}());
    const size_t half = initial_size / 2UL;
    males_.resize(half);
    females_.resize(initial_size - half);
}

void Population::run(const uint_fast32_t years) {HERE;
    urbg_t urbg(std::random_device{}());
    for (year_ = 3u; year_ < years; ++year_) {
        reproduce();
        survive();
        survive();
        survive();
        survive();
        migrate();
        std::cerr << year_ << ": " << males_.size() + females_.size() << std::endl;
    }
}

void Population::reproduce() {
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
        const uint_fast32_t num_juveniles = mother.recruitment(year_, engine_);
        const Individual* father = *wtl::choice(potential_fathers.begin(), potential_fathers.end(), engine_);
        // TODO: multiple fathers
        for (uint_fast32_t i=0; i<num_juveniles; ++i) {
            if (wtl::generate_canonical(engine_) < 0.5) {
                boys.emplace_back(*father, mother, year_);
            } else {
                girls.emplace_back(*father, mother, year_);
            }
        }
    }
    std::copy(boys.begin(), boys.end(), std::back_inserter(males_));
    std::copy(girls.begin(), girls.end(), std::back_inserter(females_));
}

void Population::survive() {
    auto impl = [this](const decltype(males_)& v) {
        decltype(males_) survivors;
        survivors.reserve(v.size());
        std::copy_if(v.begin(), v.end(), std::back_inserter(survivors),
                     [&](const Individual& x) {return x.has_survived(year_, engine_);});
        return survivors;
    };
    males_ = impl(males_);
    females_ = impl(females_);
}

void Population::migrate() {
    for (auto& x: males_) {x.migrate(year_, engine_);}
    for (auto& x: females_) {x.migrate(year_, engine_);}
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
