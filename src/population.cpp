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
    const size_t rest = initial_size - half;
    males_.reserve(half);
    females_.reserve(rest);
    auto founder = std::make_shared<Individual>();
    for (size_t i=0; i<half; ++i) {males_.emplace_back(new Individual(founder, founder, 0));}
    for (size_t i=0; i<rest; ++i) {females_.emplace_back(new Individual(founder, founder, 0));}
}

Population::~Population() {} // to allow forward declaration of Individual

void Population::run(const uint_fast32_t simulating_duration,
                     const size_t sample_size_,
                     const uint_fast32_t recording_duration) {HERE;
    auto recording_start = simulating_duration - recording_duration;
    write_sample_header(std::cout);
    for (year_ = 4u; year_ < simulating_duration; ++year_) {
        DCERR(year_ << ": " << sizes() << std::endl);
        reproduce();
        survive(0u);
        survive(1u);
        survive(2u);
        survive(3u);
        migrate();
        if (year_ >= recording_start) {
            sample(sample_size_, std::cout);
        }
    }
    DCERR(year_ << ": " << sizes() << std::endl);
}

void Population::reproduce() {
    std::vector<std::shared_ptr<Individual>> boys;
    std::vector<std::shared_ptr<Individual>> girls;
    std::vector<std::vector<std::shared_ptr<Individual>>> males_located(Individual::num_locations());
    for (const auto& p: males_) {
        males_located[p->location()].emplace_back(p);
    }
    for (const auto& mother: females_) {
        if (!mother->is_in_breeding_place()) continue;
        const auto& potential_fathers = males_located[mother->location()];
        if (potential_fathers.empty()) continue;
        const uint_fast32_t num_juveniles = mother->recruitment(year_, engine_);
        const std::shared_ptr<Individual> father = *wtl::choice(potential_fathers.begin(), potential_fathers.end(), engine_);
        // TODO: multiple fathers
        for (uint_fast32_t i=0; i<num_juveniles; ++i) {
            if (wtl::generate_canonical(engine_) < 0.5) {
                boys.emplace_back(new Individual(father, mother, year_));
            } else {
                girls.emplace_back(new Individual(father, mother, year_));
            }
        }
    }
    std::copy(boys.begin(), boys.end(), std::back_inserter(males_));
    std::copy(girls.begin(), girls.end(), std::back_inserter(females_));
}

void Population::survive(const uint_fast32_t quarter) {
    auto impl = [quarter,this](const decltype(males_)& v) {
        decltype(males_) survivors;
        survivors.reserve(v.size());
        std::copy_if(v.begin(), v.end(), std::back_inserter(survivors),
                     [&](const auto& p) {return p->has_survived(year_, quarter, engine_);});
        return survivors;
    };
    males_ = impl(males_);
    females_ = impl(females_);
}

void Population::migrate() {
    for (auto& p: males_) {p->migrate(year_, engine_);}
    for (auto& p: females_) {p->migrate(year_, engine_);}
}

void Population::sample(const size_t n, std::ostream& ost) {
    auto impl = [this,&ost](const decltype(males_)& v, const size_t sample_size) {
        decltype(males_) survivors;
        survivors.reserve(v.size() - sample_size);
        auto indices = wtl::sample(v.size(), sample_size, engine_);
        for (size_t i=0; i<v.size(); ++i) {
            if (indices.find(i) == indices.end()) {
                survivors.emplace_back(std::move(v[i]));
            } else {
                write_sample(*v[i], ost);
            }
        }
        return survivors;
    };
    males_ = impl(males_, n / 2ul);
    females_ = impl(females_, n - n / 2ul);
}

std::ostream& Population::write_sample_header(std::ostream& ost) {
    return wtl::join(Individual::names(), ost, "\t") << "\tcapture_year\n";
}

std::ostream& Population::write_sample(const Individual& x, std::ostream& ost) const {
    return ost << x << "\t" << year_ << "\n";
}

std::vector<size_t> Population::sizes() const {
    std::vector<size_t> counter(Individual::num_locations(), 0u);
    for (const auto& p: males_) {++counter[p->location()];}
    for (const auto& p: females_) {++counter[p->location()];}
    return counter;
}

std::ostream& Population::write(std::ostream& ost) const {
    for (const auto& p: males_) {write_sample(*p, ost);}
    for (const auto& p: females_) {write_sample(*p, ost);}
    return ost;
}

//! shortcut Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

} // namespace pbt
