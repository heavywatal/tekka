/*! @file population.cpp
    @brief Implementation of Population class
*/
#include "population.hpp"
#include "individual.hpp"

#include <wtl/random.hpp>
#include <wtl/debug.hpp>
#include <pcglite/pcglite.hpp>

#include <algorithm>
#include <cmath>

namespace pbf {

namespace {

template <class T> inline
void elongate(std::vector<T>& v, size_t n) noexcept {
    for (size_t i=v.size(); i<n; ++i) {
        v.emplace_back(v.back());
    }
}

template <class T> inline
void copy_elongate(const std::vector<T>& src, std::vector<T>& dst, size_t n) noexcept {
    dst.reserve(n);
    for (const auto& x: src) {
        dst.emplace_back(x);
    }
    elongate(dst, n);
}

inline uint_fast32_t get_dest(const std::vector<double>& v) {
    uint_fast32_t idx = 0;
    uint_fast32_t num_positive = 0u;
    for (uint_fast32_t i=0u; i<v.size(); ++i) {
        if (v[i] > 0.0) {
            ++num_positive;
            idx = i;
        }
    }
    return num_positive == 1u ? idx : std::numeric_limits<uint_fast32_t>::max();
}

inline uint_fast32_t sub_sat(const uint_fast32_t x, const uint_fast32_t y) {
    return (x > y) ? x - y : uint_fast32_t{};
}

} // anonymous namespace


Population::Population(const Parameters& params)
: params_{params},
  subpopulations_(4u),
  engine_{std::make_unique<URBG>(params_.seed)} {
    if (!is_ready()) {
        propagate_params();
    }
    init_demography(params_.years + 1);
    const uint_fast32_t initial_size = params_.origin ? params_.origin :
      (params_.med_recruitment ? params_.med_recruitment : params_.carrying_capacity);
    auto& subpop0 = subpopulations_[0u];
    auto origin = std::make_shared<Individual>();
    subpop0[Sex::F].assign({origin});
    subpop0[Sex::M].assign({origin});
    reproduce_impl(subpop0, {initial_size});
    subpop0.clear();
    migrate();
}

Population::~Population() = default;

void Population::run() {
    const auto num_sample_locs = std::max(params_.sample_size_adult.size(), params_.sample_size_juvenile.size());
    auto recording_start = params_.years - params_.last;
    for (year_ = 1; year_ <= params_.years; ++year_) {
        reproduce();
        for (int_fast32_t q = 0; q < 4; ++q) {
          record_demography(q);
          survive(q);
        }
        if (year_ > recording_start) {
            for (uint_fast32_t loc = 0; loc < num_sample_locs; ++loc) {
                sample(subpopulations_[loc], params_.sample_size_adult[loc], params_.sample_size_juvenile[loc]);
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
                record.resize(MAX_AGE);
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
    if (params_.med_recruitment > 0.0) {
        reproduce_lognormal();
    } else {
        reproduce_logistic();
    }
}

void Population::reproduce_lognormal() {
    constexpr uint_fast32_t num_breeding_places = 2u;
    std::lognormal_distribution lognormal{std::log(params_.med_recruitment), params_.sigma_recruitment};
    const size_t recruitment = params_.sigma_recruitment ? lognormal(*engine_) : params_.med_recruitment;
    std::vector<double> subtotal_weights(num_breeding_places);
    std::vector<std::vector<double>> female_weights(num_breeding_places);
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        const auto& females = subpopulations_[loc][Sex::F];
        female_weights[loc] = weights(females);
        subtotal_weights[loc] = std::reduce(female_weights[loc].begin(), female_weights[loc].end());
    }
    wtl::multinomial_distribution multinomial(subtotal_weights.begin(), subtotal_weights.end());
    const auto recruitment_sub = multinomial(*engine_, recruitment);
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        auto& subpop = subpopulations_[loc];
        auto& wl = female_weights[loc];
        wtl::multinomial_distribution<uint_fast32_t> multinomial{wl.begin(), wl.end()};
        reproduce_impl(subpop, multinomial(*engine_, recruitment_sub[loc]));
    }
}

void Population::reproduce_logistic() {
    constexpr uint_fast32_t num_breeding_places = 2u;
    size_t popsize = 0;
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        popsize += subpopulations_[loc].size();
    }
    for (uint_fast32_t loc=0u; loc<num_breeding_places; ++loc) {
        auto& subpop = subpopulations_[loc];
        const auto& females = subpop[Sex::F];
        reproduce_impl(subpop, litter_sizes_logistic(females, popsize));
    }
}

void Population::reproduce_impl(SubPopulation& subpop, const std::vector<uint_fast32_t>& litter_sizes) {
    const auto& females = subpop[Sex::F];
    const auto& males = subpop[Sex::M];
    if (females.empty() || males.empty()) return;
    const std::vector<double> vw = weights(males);
    std::discrete_distribution<uint_fast32_t> mate_distr(vw.begin(), vw.end());
    std::vector<double> d(4u);
    for (int_fast32_t q = 0; q < 4; ++q) {
        d[q] = death_rate(0, q);
    }
    auto& juveniles = subpop.juveniles;
    juveniles.reserve(subpop.size());
    auto& demography_year = subpop.demography[year_];
    for (size_t i = 0u; i < litter_sizes.size(); ++i) {
        auto rec = litter_sizes[i];
        for (int_fast32_t q = 0; q < 4; ++q) {
            demography_year[q][0u] += rec;
            rec -= std::binomial_distribution<uint_fast32_t>(rec, d[q])(*engine_);
        }
        auto& mother = females[i];
        for (decltype(rec) j = 0; j < rec; ++j) {
            const auto& father = males[mate_distr(*engine_)];
            juveniles.emplace_back(std::make_shared<Individual>(father, mother, year_));
        }
    }
}

std::vector<uint_fast32_t>
Population::litter_sizes_logistic(const std::vector<ShPtrIndividual>& females, const size_t popsize) {
    const auto N = static_cast<double>(popsize);
    if (N > params_.carrying_capacity) return {};
    const double rec_rate = params_.recruitment * (1.0 - N / params_.carrying_capacity);
    std::vector<uint_fast32_t> res;
    res.reserve(females.size());
    for (const auto& mother: females) {
        const double rec_mean = rec_rate * WEIGHT_FOR_AGE_[mother->age(year_)];
        res.push_back(rnbinom<uint_fast32_t>(params_.overdispersion, rec_mean, *engine_));
    }
    return res;
}

std::vector<double> Population::weights(const std::vector<ShPtrIndividual>& individuals) const {
    std::vector<double> res;
    res.reserve(individuals.size());
    for (const auto& p: individuals) {
        res.push_back(WEIGHT_FOR_AGE_[p->age(year_)]);
    }
    return res;
}

void Population::survive(const int_fast32_t season) {
    for (auto& subpop: subpopulations_) {
      for (auto& individuals: subpop.adults) {
        size_t n = individuals.size();
        for (size_t i=0; i<n; ++i) {
            auto& p = individuals[i];
            if (wtl::generate_canonical(*engine_) < death_rate(p->age(year_), season)) {
                p = std::move(individuals.back());
                individuals.pop_back();
                --n;
                --i;
            }
        }
      }
    }
}

double Population::death_rate(const int_fast32_t age, const int_fast32_t season) const {
    const auto q_age = 4 * age + season;
    const auto m = NATURAL_MORTALITY_[q_age];
    const auto f = FISHING_MORTALITY_[q_age] * FISHING_COEF_[year_];
    return 1.0 - std::exp(-m - f);
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
                uint_fast32_t dst = destination(p->age(year_), loc);
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
            const auto dst = destination(p->age(year_), loc);
            subpopulations_[dst][sex].emplace_back(std::move(p));
        }
        juveniles.clear();
    }
}

uint_fast32_t Population::destination(int_fast32_t age, uint_fast32_t loc) {
    auto& [dest, dist] = MIGRATION_DESTINATION_[age][loc];
    return (dest == std::numeric_limits<uint_fast32_t>::max()) ? dist(*engine_) : dest;
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

void Population::propagate_params() {
    init_mortality();
    init_migration();
    init_weight();
}

void Population::init_migration() {
    using dist_t = PairDestDist::second_type;
    MIGRATION_DESTINATION_.clear();
    MIGRATION_DESTINATION_.reserve(MAX_AGE);
    for (const auto& matrix: params_.migration_matrices) {
        std::vector<PairDestDist> pairs;
        pairs.reserve(matrix.size());
        for (const auto& row: matrix) {
            pairs.emplace_back(get_dest(row), dist_t(row.begin(), row.end()));
        }
        MIGRATION_DESTINATION_.emplace_back(std::move(pairs));
    }
    elongate(MIGRATION_DESTINATION_, MAX_AGE);
}

void Population::init_mortality() {
    copy_elongate(params_.natural_mortality, NATURAL_MORTALITY_, 4u * MAX_AGE);
    copy_elongate(params_.fishing_mortality, FISHING_MORTALITY_, 4u * MAX_AGE);
    NATURAL_MORTALITY_.back() = 1e9;
    const auto fc_size = static_cast<uint_fast32_t>(params_.fishing_coef.size());
    const auto offset = sub_sat(fc_size, params_.years);
    FISHING_COEF_.assign(params_.years, 1.0);
    std::copy_backward(params_.fishing_coef.begin() + offset, params_.fishing_coef.end(),
                       FISHING_COEF_.end());
}

void Population::init_weight() {
    WEIGHT_FOR_AGE_.reserve(MAX_AGE);
    WEIGHT_FOR_AGE_.resize(params_.weight_for_age.size() / 4u);
    for (size_t year=0; year<WEIGHT_FOR_AGE_.size(); ++year) {
        WEIGHT_FOR_AGE_[year] = params_.weight_for_age[4u * year];
    }
    elongate(WEIGHT_FOR_AGE_, MAX_AGE);
}

bool Population::is_ready() const {
    return (
      FISHING_COEF_.size() >= static_cast<size_t>(params_.years)
      && NATURAL_MORTALITY_.size() >= 4u * MAX_AGE
      && FISHING_MORTALITY_.size() >= 4u * MAX_AGE
      && WEIGHT_FOR_AGE_.size() >= MAX_AGE
      && MIGRATION_DESTINATION_.size() >= MAX_AGE
    );
}

} // namespace pbf
