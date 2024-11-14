/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/random.hpp>
#include <wtl/debug.hpp>
#include <pcglite/pcglite.hpp>

#include <algorithm>

namespace pbf {

Population::Population(const size_t initial_size, std::random_device::result_type seed,
  const double carrying_capacity,
  const double recruitment_coef,
  const double negative_binom_k)
: subpopulations_(4u), juveniles_subpops_(2u),
  carrying_capacity_(carrying_capacity),
  recruitment_coef_(recruitment_coef),
  k_nbinom_(negative_binom_k),
  engine_(std::make_unique<URBG>(seed)) {
    subpopulations_[0u].reserve(initial_size);
    const size_t half = initial_size / 2UL;
    for (size_t i=0; i<initial_size; ++i) {
        subpopulations_[0u].emplace_back(std::make_shared<Individual>(i < half));
    }
}

Population::~Population() = default;

void Population::run(const int_fast32_t simulating_duration,
                     const std::vector<size_t>& sample_size_adult,
                     const std::vector<size_t>& sample_size_juvenile,
                     const int_fast32_t recording_duration) {
    if (!Individual::is_ready(simulating_duration)) {
        Individual::set_dependent_static(simulating_duration);
    }
    const auto num_sample_locs = std::max(sample_size_adult.size(), sample_size_juvenile.size());
    loc_year_samples_.resize(std::min(num_subpops(), num_sample_locs));
    auto recording_start = simulating_duration - recording_duration;
    append_demography(3);
    for (year_ = 1; year_ <= simulating_duration; ++year_) {
        reproduce();
        if (year_ == 1) subpopulations_[0u].clear();
        for (int_fast32_t q = 0; q < 4; ++q) {
          append_demography(q);
          survive(q);
        }
        if (year_ > recording_start) {
            sample(subpopulations_, sample_size_adult);
            sample(juveniles_subpops_, sample_size_juvenile);
        }
        migrate();
    }
}

namespace {

template <class T> inline wtl::negative_binomial_distribution<T>
nbinom_distribution(double k, double mu) {
    const double prob = k / (mu + k);
    return wtl::negative_binomial_distribution<T>(k, prob);
}

template <class T> inline T
rnbinom(double k, double mu, URBG& engine) {
    if (k > 0.0) {
        return nbinom_distribution<T>(k, mu)(engine);
    } else {
        return std::poisson_distribution<T>(mu)(engine);
    }
}

} // anonymous namespace

void Population::reproduce() {
    const auto num_breeding_places = juveniles_subpops_.size();
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
    const auto N = static_cast<double>(popsize);
    if (N > carrying_capacity_) return;
    const double rec_rate = recruitment_coef_ * (1.0 - N / carrying_capacity_);
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
    std::vector<double> d(4u);
    double s0 = 1.1;
    for (int_fast32_t q = 0; q < 4; ++q) {
        d[q] = Individual::DEATH_RATE(q, year_);
        s0 *= (1.0 - d[q]);
    }
    juveniles.reserve(static_cast<size_t>(s0 * rec_rate * female_biomass));
    for (const auto& mother: adults) {
        if (mother->is_male()) continue;
        const double rec_mean = rec_rate * mother->weight(year_);
        uint_fast32_t rec = rnbinom<uint_fast32_t>(k_nbinom_, rec_mean, *engine_);
        for (int_fast32_t q = 0; q < 4; ++q) {
            juveniles_demography_[q][location] += rec;
            rec -= std::binomial_distribution<uint_fast32_t>(rec, d[q])(*engine_);
        }
        const auto num_boys = std::binomial_distribution<uint_fast32_t>(rec, 0.5)(*engine_);
        for (uint_fast32_t j=0; j<rec; ++j) {
            const auto& father = adults[male_indices[mate_distr(*engine_)]];
            juveniles.emplace_back(std::make_shared<Individual>(father, mother, year_, j < num_boys));
        }
    }
}

void Population::survive(const int_fast32_t season) {
    for (auto& individuals: subpopulations_) {
        size_t n = individuals.size();
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            if (wtl::generate_canonical(*engine_) < p->death_rate(year_, season)) {
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
        for (size_t i=0; i<n;) {
            auto& p = individuals[i];
            uint_fast32_t new_loc = p->destination(loc, year_, *engine_);
            if (new_loc == loc) {++i; continue;}
            subpopulations_[new_loc].emplace_back(std::move(p));
            --n;
            p = std::move(individuals[n]);
            individuals[n] = std::move(individuals.back());
            individuals.pop_back();
        }
    }
    for (uint_fast32_t loc=0; loc<juveniles_subpops_.size(); ++loc) {
        auto& juveniles = juveniles_subpops_[loc];
        for (auto& p: juveniles) {
            auto new_loc = p->destination(loc, year_, *engine_);
            subpopulations_[new_loc].emplace_back(std::move(p));
        }
        juveniles.clear();
    }
}

void Population::sample(std::vector<std::vector<std::shared_ptr<Individual>>>& subpops,
                        const std::vector<size_t>& sample_sizes) {
    const auto max_loc = std::min(subpops.size(), sample_sizes.size());
    for (uint_fast32_t loc=0u; loc<max_loc; ++loc) {
        auto& individuals = subpops.at(loc);
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
    Individual::write_trace_back_header(ost);
    // unordered_set suffices to remove duplicates, but address is not reproducible.
    std::unordered_map<const Individual*, uint_fast32_t> ids;
    ids.emplace(nullptr, 0u);
    for (uint_fast32_t loc=0u; loc<loc_year_samples_.size(); ++loc) {
        for (const auto& [_year, samples]: loc_year_samples_.at(loc)) {
            for (const auto& p: samples) {
                ids.emplace(p.get(), static_cast<uint_fast32_t>(ids.size()));
            }
        }
    }
    for (uint_fast32_t loc=0u; loc<loc_year_samples_.size(); ++loc) {
        for (const auto& [year, samples]: loc_year_samples_.at(loc)) {
            for (const auto& p: samples) {
                p->trace_back(ost, ids, loc, year);
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
    for (const auto& [time, values]: demography_) {
        const auto [year, season] = time;
        for (uint_fast32_t loc=0; loc<num_subpops(); ++loc) {
            const auto& structure = values[loc];
            for (uint_fast32_t age=0u; age<structure.size(); ++age) {
                if (structure[age] == 0u) continue;
                ost << year << "\t" << season << "\t"
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
