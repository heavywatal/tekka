/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/random.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/exception.hpp>

namespace pbf {

Population::Population(const size_t initial_size, std::random_device::result_type seed)
: engine_(std::make_unique<URBG>(seed)) {
    individuals_.reserve(initial_size);
    const size_t half = initial_size / 2UL;
    for (size_t i=0; i<initial_size; ++i) {
        individuals_.emplace_back(std::make_shared<Individual>(i < half));
    }
}

Population::~Population() = default;

void Population::run(const int_fast32_t simulating_duration,
                     const std::vector<size_t>& sample_size_adult,
                     const std::vector<size_t>& sample_size_juvenile,
                     const int_fast32_t recording_duration) {
    auto recording_start = simulating_duration - recording_duration;
    append_demography(0);
    for (year_ = 0; year_ < simulating_duration; ++year_) {
        reproduce();
        survive(0);
        survive(1);
        survive(2);
        survive(3);
        merge_juveniles();
        if (year_ >= recording_start) {
            sample(sample_size_adult, sample_size_juvenile);
        }
        migrate();
    }
    append_demography(0);
}

void Population::reproduce() {
    double female_biomass = 0.0;
    double male_biomass = 0.0;
    std::vector<std::vector<std::shared_ptr<Individual>>> males_located(Individual::num_breeding_places());
    for (const auto& p: individuals_) {
        if (p->is_in_breeding_place()) {
            if (p->is_male()) {
                male_biomass += p->weight(year_);
                males_located[p->location()].emplace_back(p);
            } else {
                female_biomass += p->weight(year_);
            }
        }
    }
    std::vector<std::discrete_distribution<unsigned>> mate_distrs;
    mate_distrs.reserve(males_located.size());
    for (size_t i=0; i<males_located.size(); ++i) {
        std::vector<double> fitnesses;
        fitnesses.reserve(males_located[i].size());
        for (const auto& male: males_located[i]) {
            fitnesses.push_back(male->weight(year_));
        }
        mate_distrs.emplace_back(fitnesses.begin(), fitnesses.end());
    }
    const double biomass = female_biomass + male_biomass;
    const double popsize = biomass / Individual::weight_for_age().back();
    const double density_effect = std::max(0.0, 1.0 - popsize / Individual::param().CARRYING_CAPACITY);
    const double exp_recruitment = density_effect * Individual::param().RECRUITMENT_COEF * female_biomass;
    const size_t n = individuals_.size();
    juveniles_.reserve(static_cast<size_t>(exp_recruitment * 1.1));
    juveniles_demography_.assign(4u, std::vector<uint_fast32_t>(Individual::num_breeding_places()));
    for (size_t i=0; i<n; ++i) {
        const auto& mother = individuals_[i];
        if (mother->is_male() || !mother->is_in_breeding_place()) continue;
        const auto& potential_fathers = males_located[mother->location()];
        if (potential_fathers.empty()) continue;
        auto& mate_distr = mate_distrs[mother->location()];
        uint_fast32_t num_juveniles = mother->recruitment(year_, density_effect, *engine_);
        for (unsigned season: {0u, 1u, 2u, 3u}) {
            std::binomial_distribution<uint_fast32_t> binom(num_juveniles, Individual::survival_rate()[season]);
            num_juveniles = binom(*engine_);
            juveniles_demography_[season][mother->location()] += num_juveniles;
        }
        for (uint_fast32_t i=0; i<num_juveniles; ++i) {
            const auto& father = potential_fathers[mate_distr(*engine_)];
            const bool is_male = (wtl::generate_canonical(*engine_) < 0.5);
            juveniles_.emplace_back(std::make_shared<Individual>(father, mother, year_, is_male));
        }
    }
}

void Population::survive(const int_fast32_t season) {
    size_t n = individuals_.size();
    for (size_t i=0; i<n; ++i) {
        auto& p = individuals_[i];
        while (!p->has_survived(year_, season, *engine_)) {
            if (i < n) {
                p = std::move(individuals_.back());
            }
            individuals_.pop_back();
            --n;
        }
    }
    append_demography(season);
}

void Population::merge_juveniles() {
    individuals_.reserve(individuals_.size() + juveniles_.size());
    for (auto& p: juveniles_) {
        individuals_.emplace_back(std::move(p));
    }
    juveniles_.clear();
}

void Population::migrate() {
    for (auto& p: individuals_) {p->migrate(year_, *engine_);}
}

void Population::sample(std::vector<size_t> sample_size_adult,
                        std::vector<size_t> sample_size_juvenile) {
    size_t total_sampled = 0u;
    total_sampled += std::accumulate(sample_size_juvenile.begin(), sample_size_juvenile.end(), 0u);
    total_sampled += std::accumulate(sample_size_adult.begin(), sample_size_adult.end(), 0u);
    sample_size_juvenile.resize(Individual::num_locations());
    sample_size_adult.resize(Individual::num_locations());
    std::vector<std::shared_ptr<Individual>>& sampled = year_samples_[year_];
    sampled.reserve(total_sampled);
    std::shuffle(individuals_.begin(), individuals_.end(), *engine_);
    size_t n = individuals_.size();
    for (size_t i=0; i<n; ++i) {
        auto& p = individuals_[i];
        auto& sample_size = (p->birth_year() == year_) ? sample_size_juvenile : sample_size_adult;
        auto& to_be_sampled = sample_size[p->location()];
        if (to_be_sampled) {
            sampled.emplace_back(p);
            --to_be_sampled;
            if (i < n) {
                p = std::move(individuals_.back());
            }
            individuals_.pop_back();
            --n;
            --i;
        }
    }
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    if (year_samples_.empty()) return ost;
    wtl::join(Individual::names(), ost, "\t") << "\tcapture_year\n";
    std::map<const Individual*, size_t> ids;
    ids.emplace(nullptr, 0u);
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            ids.emplace(p.get(), ids.size());
        }
    }
    for (const auto& ys: year_samples_) {
        for (const auto& p: ys.second) {
            p->trace_back(ost, &ids, ys.first);
        }
    }
    return ost;
}

std::vector<std::map<int_fast32_t, size_t>> Population::count(const int_fast32_t season) const {
    std::vector<std::map<int_fast32_t, size_t>> counter(Individual::num_locations());
    if (!juveniles_demography_.empty()) {
        const auto& jd_season = juveniles_demography_.at(season);
        for (unsigned loc=0; loc<jd_season.size(); ++loc) {
            counter[loc][0u] += jd_season[loc];
        }
    }
    for (const auto& p: individuals_) {
        ++counter[p->location()][year_ - p->birth_year()];
    }
    return counter;
}

void Population::append_demography(const int_fast32_t season) {
    demography_.emplace_hint(
      demography_.end(),
      std::pair<int_fast32_t, int_fast32_t>(year_, season),
      count(season)
    );
}

std::ostream& Population::write_demography(std::ostream& ost) const {
    ost << "year\tseason\tlocation\tage\tcount\n";
    for (const auto& time_structure: demography_) {
        const auto time = time_structure.first;
        for (size_t loc=0; loc<Individual::num_locations(); ++loc) {
            const auto structure = time_structure.second[loc];
            for (const auto& age_count: structure) {
                ost << time.first << "\t" << time.second << "\t"
                    << loc << "\t"
                    << age_count.first << "\t"
                    << age_count.second << "\n";
            }
        }
    }
    return ost;
}

std::ostream& Population::write(std::ostream& ost) const {
    for (const auto& p: individuals_) {ost << *p << "\n";}
    return ost;
}

//! shortcut Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

} // namespace pbf
