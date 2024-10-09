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
: subpopulations_(4u), juveniles_subpops_(2u),
  engine_(std::make_unique<URBG>(seed)) {
    subpopulations_[0u].reserve(initial_size);
    const size_t half = initial_size / 2UL;
    for (size_t i=0; i<initial_size; ++i) {
        subpopulations_[0u].emplace_back(std::make_shared<Individual>(i < half));
    }
}

void Population::run(const int_fast32_t simulating_duration,
                     const std::vector<size_t>& sample_size_adult,
                     const std::vector<size_t>& sample_size_juvenile,
                     const int_fast32_t recording_duration) {
    if (!Individual::is_ready(simulating_duration)) {
        Individual::set_dependent_static(simulating_duration);
    }
    loc_year_samples_.resize(std::min(num_subpops(),
                                      std::max(sample_size_adult.size(),
                                               sample_size_juvenile.size())));
    auto recording_start = simulating_duration - recording_duration;
    append_demography(3);
    for (year_ = 1; year_ <= simulating_duration; ++year_) {
        reproduce();
        if (year_ == 1) subpopulations_[0u].clear();
        append_demography(0);
        survive();
        if (year_ > recording_start) {
            sample(&subpopulations_, sample_size_adult);
            sample(&juveniles_subpops_, sample_size_juvenile);
        }
        migrate();
        append_demography(3);
    }
}

void Population::reproduce() {
    const auto num_breeding_places = static_cast<uint_fast32_t>(juveniles_subpops_.size());
    juveniles_demography_.assign(4u, std::vector<uint_fast32_t>(num_breeding_places));
    size_t popsize = 0;
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        popsize += subpopulations_[loc].size();
    }
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        reproduce(loc, popsize);
    }
}

void Population::reproduce(const uint_fast32_t location, const size_t popsize) {
    const double K = Individual::param().CARRYING_CAPACITY;
    const double density_effect = 1.0 - popsize / K;
    if (density_effect <= 0.0) return;
    const auto& adults = subpopulations_[location];
    auto& juveniles = juveniles_subpops_[location];
    double female_biomass = 0.0;
    std::vector<uint_fast32_t> male_indices;
    std::vector<double> fitnesses;
    const size_t males_capa = (6u * adults.size() / 10u);
    male_indices.reserve(males_capa);
    fitnesses.reserve(males_capa);
    for (uint_fast32_t i=0u; i<adults.size(); ++i) {
        const auto& p = adults[i];
        if (p->is_male()) {
            male_indices.push_back(i);
            fitnesses.push_back(p->weight(year_));
        } else {
            female_biomass += p->weight(year_);
        }
    }
    if (male_indices.size() == 0u) return;
    std::discrete_distribution<uint_fast32_t> mate_distr(fitnesses.begin(), fitnesses.end());
    const double exp_recruitment = density_effect * Individual::param().RECRUITMENT_COEF * female_biomass;
    juveniles.reserve(static_cast<size_t>(exp_recruitment * 1.1));
    const double d0 = Individual::death_rate(0, year_);
    for (const auto& mother: adults) {
        if (mother->is_male()) continue;
        uint_fast32_t num_juveniles = mother->recruitment(year_, density_effect, *engine_);
        juveniles_demography_[0u][location] += num_juveniles;
        num_juveniles -= std::binomial_distribution<uint_fast32_t>(num_juveniles, d0)(*engine_);
        juveniles_demography_[3u][location] += num_juveniles;
        const auto num_boys = std::binomial_distribution<uint_fast32_t>(num_juveniles, 0.5)(*engine_);
        for (uint_fast32_t j=0; j<num_juveniles; ++j) {
            const auto& father = adults[male_indices[mate_distr(*engine_)]];
            juveniles.emplace_back(std::make_shared<Individual>(father, mother, year_, j < num_boys));
        }
    }
}

void Population::survive() {
    for (auto& individuals: subpopulations_) {
        size_t n = individuals.size();
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            if (wtl::generate_canonical(*engine_) < p->death_rate(year_)) {
                p = std::move(individuals.back());
                individuals.pop_back();
                --n;
                --i;
            }
        }
    }
}

void Population::migrate() {
    std::vector<size_t> subpop_sizes(num_subpops());
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        subpop_sizes[loc] = subpopulations_[loc].size();
    }
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        auto& individuals = subpopulations_[loc];
        size_t n = subpop_sizes[loc];
        size_t num_immigrants = individuals.size() - n;
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            auto new_loc = p->migrate(loc, year_, *engine_);
            if (new_loc == loc) continue;
            subpopulations_[new_loc].emplace_back(std::move(p));
            p = std::move(individuals.back());
            individuals.pop_back();
            if (num_immigrants == 0u) {
                --n;
                --i;
            } else {
                --num_immigrants;
            }
        }
    }
    for (uint_fast32_t loc=0; loc<juveniles_subpops_.size(); ++loc) {
        auto& juveniles = juveniles_subpops_[loc];
        for (auto& p: juveniles) {
            auto new_loc = p->migrate(loc, year_, *engine_);
            subpopulations_[new_loc].emplace_back(std::move(p));
        }
        juveniles.clear();
    }
}

void Population::sample(std::vector<std::vector<std::shared_ptr<Individual>>>* subpops,
                        const std::vector<size_t>& sample_sizes) {
    const auto max_loc = std::min(subpops->size(), sample_sizes.size());
    for (uint_fast32_t loc=0u; loc<max_loc; ++loc) {
        auto& individuals = subpops->at(loc);
        std::shuffle(individuals.begin(), individuals.end(), *engine_);
        const auto n = std::min(individuals.size(), sample_sizes[loc]);
        std::vector<std::shared_ptr<Individual>>& sampled = loc_year_samples_[loc][year_];
        sampled.reserve(sampled.size() + n);
        for (size_t i=0; i<n; ++i) {
            sampled.emplace_back(std::move(individuals.back()));
            individuals.pop_back();
        }
    }
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    if (loc_year_samples_.empty() || loc_year_samples_[0u].empty()) return ost;
    wtl::join(Individual::names(), ost, "\t") << "\tlocation\tcapture_year\n";
    // unordered_set suffices to remove duplicates, but address is not reproducible.
    std::unordered_map<const Individual*, uint_fast32_t> ids;
    ids.emplace(nullptr, 0u);
    for (uint_fast32_t loc=0u; loc<loc_year_samples_.size(); ++loc) {
        for (const auto& ys: loc_year_samples_.at(loc)) {
            for (const auto& p: ys.second) {
                ids.emplace(p.get(), static_cast<uint_fast32_t>(ids.size()));
            }
        }
    }
    for (uint_fast32_t loc=0u; loc<loc_year_samples_.size(); ++loc) {
        for (const auto& ys: loc_year_samples_.at(loc)) {
            for (const auto& p: ys.second) {
                p->trace_back(ost, &ids, loc, ys.first);
            }
        }
    }
    return ost;
}

std::vector<std::vector<uint_fast32_t>> Population::count(const int_fast32_t season) const {
    std::vector<std::vector<uint_fast32_t>> counter(num_subpops(), std::vector<uint_fast32_t>(Individual::MAX_AGE));
    if (!juveniles_demography_.empty()) {
        const auto& jd_season = juveniles_demography_.at(season);
        for (uint_fast32_t loc=0; loc<jd_season.size(); ++loc) {
            counter[loc][0u] += jd_season[loc];
        }
    }
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        auto& counter_loc = counter[loc];
        for (const auto& p: subpopulations_[loc]) {
            ++counter_loc[p->age(year_)];
        }
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
        for (uint_fast32_t loc=0; loc<num_subpops(); ++loc) {
            const auto structure = time_structure.second[loc];
            for (uint_fast32_t age=0u; age<structure.size(); ++age) {
                if (structure[age] == 0u) continue;
                ost << time.first << "\t" << time.second << "\t"
                    << loc << "\t"
                    << age << "\t"
                    << structure[age] << "\n";
            }
        }
    }
    return ost;
}

std::ostream& Population::write(std::ostream& ost) const {
    for (const auto& individuals: juveniles_subpops_) {
        for (const auto& p: individuals) {ost << *p << "\n";}
    }
    for (const auto& individuals: subpopulations_) {
        for (const auto& p: individuals) {ost << *p << "\n";}
    }
    return ost;
}

//! Shortcut for Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

} // namespace pbf
