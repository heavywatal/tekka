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
: subpopulations_(4u), juveniles_subpops_(2u), loc_year_samples_(4u),
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
    auto recording_start = simulating_duration - recording_duration;
    append_demography(3);
    for (year_ = 1; year_ <= simulating_duration; ++year_) {
        reproduce();
        if (year_ == 1) subpopulations_[0u].clear();
        append_demography(0);
        survive();
        merge_juveniles();
        if (year_ > recording_start) {
            sample(sample_size_adult, sample_size_juvenile);
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
    const double density_effect = std::max(0.0, 1.0 - popsize / Individual::param().CARRYING_CAPACITY);
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        reproduce(loc, density_effect);
    }
}

void Population::reproduce(const uint_fast32_t location, const double density_effect) {
    const auto& adults = subpopulations_[location];
    auto& juveniles = juveniles_subpops_[location];
    const size_t n = adults.size();
    double female_biomass = 0.0;
    std::vector<uint_fast32_t> male_indices;
    std::vector<double> fitnesses;
    male_indices.reserve(static_cast<size_t>(adults.size() * 0.6));
    fitnesses.reserve(static_cast<size_t>(adults.size() * 0.6));
    for (uint_fast32_t i=0u; i<n; ++i) {
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
    const double survival_rate = 1.0 - Individual::death_rate()[0u];
    for (size_t i=0; i<n; ++i) {
        const auto& mother = adults[i];
        if (mother->is_male()) continue;
        uint_fast32_t num_juveniles = mother->recruitment(year_, density_effect, *engine_);
        juveniles_demography_[0u][location] += num_juveniles;
        std::binomial_distribution<uint_fast32_t> binom(num_juveniles, survival_rate);
        num_juveniles = binom(*engine_);
        juveniles_demography_[3u][location] += num_juveniles;
        for (uint_fast32_t j=0; j<num_juveniles; ++j) {
            const auto& father = adults[male_indices[mate_distr(*engine_)]];
            const bool is_male = (wtl::generate_canonical(*engine_) < 0.5);
            juveniles.emplace_back(std::make_shared<Individual>(father, mother, year_, is_male));
        }
    }
}

void Population::survive() {
    for (auto& individuals: subpopulations_) {
        size_t n = individuals.size();
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            if (p->is_dead(year_, *engine_)) {
                p = std::move(individuals.back());
                individuals.pop_back();
                --n;
                --i;
            }
        }
    }
}

void Population::merge_juveniles() {
    for (uint_fast32_t loc=0; loc<juveniles_subpops_.size(); ++loc) {
        auto& individuals = subpopulations_[loc];
        auto& juveniles = juveniles_subpops_[loc];
        individuals.reserve(individuals.size() + juveniles.size());
        for (auto& p: juveniles) {
            individuals.emplace_back(std::move(p));
        }
        juveniles.clear();
    }
}

void Population::migrate() {
    std::vector<size_t> subpopsizes(num_subpops());
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        subpopsizes[loc] = subpopulations_[loc].size();
    }
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        auto& individuals = subpopulations_[loc];
        size_t n = subpopsizes[loc];
        size_t num_immigrants = individuals.size() - n;
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            auto newloc = p->migrate(loc, year_, *engine_);
            if (newloc == loc) continue;
            subpopulations_[newloc].emplace_back(std::move(p));
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
}

void Population::sample(std::vector<size_t> sample_size_adult,
                        std::vector<size_t> sample_size_juvenile) {
    size_t total_sampled = 0u;
    total_sampled += std::accumulate(sample_size_juvenile.begin(), sample_size_juvenile.end(), 0u);
    total_sampled += std::accumulate(sample_size_adult.begin(), sample_size_adult.end(), 0u);
    sample_size_juvenile.resize(num_subpops());
    sample_size_adult.resize(num_subpops());
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        std::vector<std::shared_ptr<Individual>>& sampled = loc_year_samples_[loc][year_];
        sampled.reserve(total_sampled);
        auto& individuals = subpopulations_[loc];
        std::shuffle(individuals.begin(), individuals.end(), *engine_);
        size_t n = individuals.size();
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            auto& sample_size = (p->birth_year() == year_) ? sample_size_juvenile : sample_size_adult;
            auto& to_be_sampled = sample_size[loc];
            if (to_be_sampled) {
                sampled.emplace_back(p);
                --to_be_sampled;
                p = std::move(individuals.back());
                individuals.pop_back();
                --n;
                --i;
            }
        }
    }
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    if (loc_year_samples_.empty() || loc_year_samples_[0u].empty()) return ost;
    wtl::join(Individual::names(), ost, "\t") << "\tlocation\tcapture_year\n";
    std::map<const Individual*, size_t> ids;
    ids.emplace(nullptr, 0u);
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        const auto& year_samples = loc_year_samples_[loc];
        for (const auto& ys: year_samples) {
            for (const auto& p: ys.second) {
                ids.emplace(p.get(), ids.size());
            }
        }
        for (const auto& ys: year_samples) {
            for (const auto& p: ys.second) {
                p->trace_back(ost, &ids, loc, ys.first);
            }
        }
    }
    return ost;
}

std::vector<std::vector<size_t>> Population::count(const int_fast32_t season) const {
    std::vector<std::vector<size_t>> counter(num_subpops(), std::vector<size_t>(80u));
    if (!juveniles_demography_.empty()) {
        const auto& jd_season = juveniles_demography_.at(season);
        for (uint_fast32_t loc=0; loc<jd_season.size(); ++loc) {
            counter[loc][0u] += jd_season[loc];
        }
    }
    for (uint_fast32_t loc=0u; loc<num_subpops(); ++loc) {
        auto& counter_loc = counter[loc];
        for (const auto& p: subpopulations_[loc]) {
            ++counter_loc[year_ - p->birth_year()];
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

//! shortcut Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

} // namespace pbf
