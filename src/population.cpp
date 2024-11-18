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
  const Parameters& params,
  const int_fast32_t simulating_duration,
  const double carrying_capacity,
  const double recruitment_coef,
  const double negative_binom_k)
: subpopulations_(4u),
  simulating_duration_{simulating_duration},
  carrying_capacity_{carrying_capacity},
  recruitment_coef_{recruitment_coef},
  k_nbinom_{negative_binom_k},
  engine_{std::make_unique<URBG>(seed)} {
    auto& subpop0 = subpopulations_[0u];
    subpop0[Sex::F].reserve(initial_size);
    subpop0[Sex::M].reserve(initial_size);
    const size_t half = initial_size / 2UL;
    for (size_t i=0; i<half; ++i) {
        subpop0[Sex::F].emplace_back(std::make_shared<Individual>());
    }
    for (size_t i=half; i<initial_size; ++i) {
        subpop0[Sex::M].emplace_back(std::make_shared<Individual>());
    }
    if (!Individual::is_ready(simulating_duration)) {
        Individual::set_dependent_static(params, simulating_duration);
    }
}

Population::~Population() = default;

void Population::run(const std::vector<size_t>& sample_size_adult,
                     const std::vector<size_t>& sample_size_juvenile,
                     const int_fast32_t recording_duration) {
    const auto num_sample_locs = std::max(sample_size_adult.size(), sample_size_juvenile.size());
    auto recording_start = simulating_duration_ - recording_duration;
    init_demography(simulating_duration_ + 1);
    record_demography(3);
    for (year_ = 1; year_ <= simulating_duration_; ++year_) {
        reproduce();
        if (year_ == 1) subpopulations_[0u].clear();
        for (int_fast32_t q = 0; q < 4; ++q) {
          record_demography(q);
          survive(q);
        }
        if (year_ > recording_start) {
            for (uint_fast32_t loc = 0; loc < num_sample_locs; ++loc) {
                sample(subpopulations_[loc], sample_size_adult[loc], sample_size_juvenile[loc]);
            }
        }
        migrate();
    }
}

void Population::init_demography(const int_fast32_t duration) {
    for (auto& subpop: subpopulations_) {
        subpop.demography.resize(duration);
        for (auto& seasons: subpop.demography) {
            for (auto& record: seasons) {
                record.resize(Individual::MAX_AGE);
            }
        }
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
    constexpr uint_fast32_t num_breeding_places = 2u;
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
    auto& subpop = subpopulations_[location];
    const auto& females = subpop[Sex::F];
    const auto& males = subpop[Sex::M];
    if (females.size() == 0u || males.size() == 0u) return;
    const std::vector<double> vw = weights(males);
    std::discrete_distribution<uint_fast32_t> mate_distr(vw.begin(), vw.end());
    std::vector<double> d(4u);
    for (int_fast32_t q = 0; q < 4; ++q) {
        d[q] = Individual::DEATH_RATE(q, year_);
    }
    auto& juveniles = subpop.juveniles;
    juveniles.reserve(subpop.size());
    auto& demography_year = subpop.demography[year_];
    for (const auto& mother: females) {
        const double rec_mean = rec_rate * mother->weight(year_);
        uint_fast32_t rec = rnbinom<uint_fast32_t>(k_nbinom_, rec_mean, *engine_);
        for (int_fast32_t q = 0; q < 4; ++q) {
            demography_year[q][0u] += rec;
            rec -= std::binomial_distribution<uint_fast32_t>(rec, d[q])(*engine_);
        }
        for (uint_fast32_t j=0; j<rec; ++j) {
            const auto& father = males[mate_distr(*engine_)];
            juveniles.emplace_back(std::make_shared<Individual>(father, mother, year_));
        }
    }
}

std::vector<double> Population::weights(const std::vector<ShPtrIndividual>& individuals) const {
    std::vector<double> res;
    res.reserve(individuals.size());
    for (const auto& p: individuals) {
        res.push_back(p->weight(year_));
    }
    return res;
}

void Population::survive(const int_fast32_t season) {
    for (auto& subpop: subpopulations_) {
      for (auto& individuals: subpop.adults) {
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
}

void Population::migrate() {
    for (const auto sex: {Sex::F, Sex::M}) {
        std::vector<size_t> subpop_sizes;
        subpop_sizes.reserve(subpopulations_.size());
        for (const auto& subpop: subpopulations_) {
            subpop_sizes.push_back(subpop[sex].size());
        }
        for (uint_fast32_t loc=0u; loc<subpopulations_.size(); ++loc) {
            auto& individuals = subpopulations_[loc][sex];
            size_t n = subpop_sizes[loc];
            for (size_t i=0; i<n;) {
                auto& p = individuals[i];
                uint_fast32_t dst = p->destination(loc, year_, *engine_);
                if (dst == loc) {++i; continue;}
                subpopulations_[dst][sex].emplace_back(std::move(p));
                if (--n)
                    p = std::move(individuals[n]);
                if (individuals.back())
                    individuals[n] = std::move(individuals.back());
                individuals.pop_back();
            }
        }
    }
    for (uint_fast32_t loc=0; loc<subpopulations_.size(); ++loc) {
        auto& juveniles = subpopulations_[loc].juveniles;
        for (auto& p: juveniles) {
            constexpr URBG::result_type half_max = URBG::max() >> 1u;
            const auto sex = engine_->operator()() < half_max ? Sex::F : Sex::M;
            const auto dst = p->destination(loc, year_, *engine_);
            subpopulations_[dst][sex].emplace_back(std::move(p));
        }
        juveniles.clear();
    }
}

void Population::sample(SubPopulation& subpop, size_t num_adults, size_t num_juveniles) {
    std::binomial_distribution<size_t> binom(num_adults, 0.5);
    auto num_males = binom(*engine_);
    sample(subpop[Sex::F], subpop.samples[year_], num_adults - num_males);
    sample(subpop[Sex::M], subpop.samples[year_], num_males);
    sample(subpop.juveniles, subpop.samples[year_], num_juveniles);
}

void Population::sample(std::vector<ShPtrIndividual>& src,
                        std::vector<ShPtrIndividual>& dst, size_t n) {
    if (n > src.size()) {
        std::cerr << "WARNING:Population::sample(): n > src.size() ("
                  << n << " > " << src.size() << ")\n";
        for (auto& p: src) dst.emplace_back(std::move(p));
        src.clear();
        return;
    }
    for (size_t i = 0u; i < n; ++i) {
        std::uniform_int_distribution<size_t> uniform(0u, src.size() - 1u);
        auto& p = src[uniform(*engine_)];
        dst.emplace_back(std::move(p));
        p = std::move(src.back());
        src.pop_back();
    }
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    Individual::write_trace_back_header(ost);
    // unordered_set suffices to remove duplicates, but address is not reproducible.
    std::unordered_map<const Individual*, uint_fast32_t> ids;
    ids.emplace(nullptr, 0u);
    for (uint_fast32_t loc=0u; loc<subpopulations_.size(); ++loc) {
        for (const auto& [_year, samples]: subpopulations_[loc].samples) {
            for (const auto& p: samples) {
                ids.emplace(p.get(), static_cast<uint_fast32_t>(ids.size()));
            }
        }
    }
    for (uint_fast32_t loc=0u; loc<subpopulations_.size(); ++loc) {
        for (const auto& [year, samples]: subpopulations_[loc].samples) {
            for (const auto& p: samples) {
                p->trace_back(ost, ids, loc, year);
            }
        }
    }
    return ost;
}

void Population::record_demography(const int_fast32_t season) {
    for (auto& subpop: subpopulations_) {
        auto& counter = subpop.demography[year_][season];
        for (const auto& inds: subpop.adults) {
            for (const auto& p: inds) {
                ++counter[p->age(year_)];
            }
        }
    }
}

std::ostream& Population::write_demography(std::ostream& ost) const {
    ost << "year\tseason\tlocation\tage\tcount\n";
    for (uint_fast32_t loc = 0u; loc < subpopulations_.size(); ++loc) {
      const auto& subpop = subpopulations_[loc];
      for (uint_fast32_t year = 0u; year < subpop.demography.size(); ++year) {
        const auto& demography_year = subpop.demography[year];
        for (uint_fast32_t season = 0u; season < demography_year.size(); ++season) {
            const auto& structure = demography_year[season];
            for (uint_fast32_t age=0u; age<structure.size(); ++age) {
                if (structure[age] == 0u) continue;
                ost << year << "\t" << season << "\t"
                    << loc << "\t"
                    << age << "\t"
                    << structure[age] << "\n";
            }
        }
      }
    }
    return ost;
}

std::ostream& Population::write(std::ostream& ost) const {
    for (const auto& subpop: subpopulations_) {
        for (const auto& individuals: subpop.adults) {
            for (const auto& p: individuals) {ost << *p << "\n";}
        }
        for (const auto& p: subpop.juveniles) {ost << *p << "\n";}
    }
    return ost;
}

//! Shortcut for Population::write(ost)
std::ostream& operator<<(std::ostream& ost, const Population& pop) {
    return pop.write(ost);
}

} // namespace pbf
