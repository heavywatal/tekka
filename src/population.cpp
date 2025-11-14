#include "population.hpp"
#include "individual.hpp"

#include <wtl/random.hpp>
#include <wtl/debug.hpp>
#include <pcglite/pcglite.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace pbf {

namespace {

// C++20 <iterator>
template <class T> inline
constexpr auto ssize(const T& x) {
    return static_cast<std::make_signed_t<decltype(x.size())>>(x.size());
}

template <class T> inline
constexpr T ramp(const T x) {
    return x < T{} ? T{} : x;
}

template <class T> inline
constexpr std::make_unsigned_t<T> cast_u(T x) {
    return static_cast<std::make_unsigned_t<T>>(x);
}

template <class T> inline
void elongate(std::vector<T>& v, ptrdiff_t n) noexcept {
    for (ptrdiff_t i=ssize(v); i<n; ++i) {
        v.emplace_back(v.back());
    }
}

template <class T> inline
void copy_elongate(const std::vector<T>& src, std::vector<T>& dst, ptrdiff_t n) noexcept {
    dst.reserve(cast_u(n));
    for (const auto& x: src) {
        dst.emplace_back(x);
    }
    elongate(dst, n);
}

inline int_fast32_t get_dest(const std::vector<double>& v) {
    int_fast32_t idx = 0;
    int_fast32_t num_positive = 0;
    for (int_fast32_t i=0; i<ssize(v); ++i) {
        if (v[cast_u(i)] > 0.0) {
            ++num_positive;
            idx = i;
        }
    }
    return num_positive == 1 ? idx : std::numeric_limits<int_fast32_t>::max();
}

} // anonymous namespace


Population::Population(const Parameters& params)
: params_{params},
  subpopulations_(4u),
  engine_{std::make_unique<URBG>(params_.seed)} {
    if (!is_ready()) {
        propagate_params();
    }
    params_.sample_size_adult.resize(subpopulations_.size());
    params_.sample_size_juvenile.resize(subpopulations_.size());
    init_demography(params_.years + 1);
    const int_fast32_t initial_size = params_.origin ? params_.origin :
      static_cast<int_fast32_t>(
        params_.med_recruitment > 0.0 ? params_.med_recruitment : params_.carrying_capacity
      );
    auto& subpop0 = subpopulations_[0u];
    auto origin = std::make_shared<Individual>();
    subpop0[Sex::F].assign({origin});
    subpop0[Sex::M].assign({origin});
    reproduce_impl(0, {initial_size});
    subpop0.clear();
    migrate();
}

Population::~Population() = default;

void Population::run() {
    for (year_ = 1; year_ <= params_.years; ++year_) {
        reproduce();
        for (int_fast32_t q = 0; q < 4; ++q) {
          record_demography(q);
          survive(q);
        }
        migrate();
    }
}

void Population::init_demography(const int_fast32_t duration) {
    for (auto& subpop: subpopulations_) {
        subpop.demography.resize(cast_u(duration));
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
    constexpr int_fast32_t num_breeding_places = 2;
    std::lognormal_distribution lognormal{std::log(params_.med_recruitment), params_.sigma_recruitment};
    const auto recruitment = static_cast<int_fast32_t>(
      params_.sigma_recruitment > 0.0 ? lognormal(*engine_) : params_.med_recruitment
    );
    int_fast32_t num_mothers = 0;
    std::vector<double> subtotal_weights(num_breeding_places);
    std::vector<std::vector<double>> female_weights(num_breeding_places);
    for (int_fast32_t loc=0; loc<num_breeding_places; ++loc) {
        const auto& females = subpopulations_[cast_u(loc)][Sex::F];
        num_mothers += ssize(females);
        female_weights[cast_u(loc)] = weights(females);
        subtotal_weights[cast_u(loc)] = std::reduce(female_weights[cast_u(loc)].begin(), female_weights[cast_u(loc)].end());
    }
    if (num_mothers == 0) return;
    wtl::multinomial_distribution<int_fast32_t> dist_rec(subtotal_weights.begin(), subtotal_weights.end());
    const auto recruitment_sub = dist_rec(*engine_, recruitment);
    for (int_fast32_t loc=0; loc<num_breeding_places; ++loc) {
        auto& wl = female_weights[cast_u(loc)];
        wtl::multinomial_distribution<int_fast32_t> dist_mothers{wl.begin(), wl.end()};
        reproduce_impl(loc, dist_mothers(*engine_, recruitment_sub[cast_u(loc)]));
    }
}

void Population::reproduce_logistic() {
    constexpr int_fast32_t num_breeding_places = 2;
    int_fast32_t popsize = 0;
    for (int_fast32_t loc=0; loc<num_breeding_places; ++loc) {
        popsize += ssize(subpopulations_[cast_u(loc)]);
    }
    for (int_fast32_t loc=0; loc<num_breeding_places; ++loc) {
        auto& subpop = subpopulations_[cast_u(loc)];
        reproduce_impl(loc, litter_sizes_logistic(subpop[Sex::F], popsize));
    }
}

void Population::reproduce_impl(int_fast32_t loc, std::vector<int_fast32_t> litter_sizes) {
    auto& subpop = subpopulations_[cast_u(loc)];
    const auto& females = subpop[Sex::F];
    const auto& males = subpop[Sex::M];
    if (females.empty() || males.empty()) return;
    const std::vector<double> vw = weights(males);
    std::discrete_distribution<int_fast32_t> mate_distr(vw.begin(), vw.end());
    const auto sample_size = sample_size_juvenile(loc);
    const auto dead3 = survive_juvenile(litter_sizes, subpop.demography[cast_u(year_)], sample_size);
    if (sample_size > 0) {
      auto& samples_year = subpop.samples[year_];
      std::discrete_distribution<int_fast32_t> mother_distr(dead3.begin(), dead3.end());
      for (int_fast32_t i = 0; i < sample_size; ++i) {
        const auto& mother = females[cast_u(mother_distr(*engine_))];
        const auto& father = males[cast_u(mate_distr(*engine_))];
        samples_year.emplace_back(std::make_shared<Individual>(father, mother, year_));
      }
    }
    auto& juveniles = subpop.juveniles;
    juveniles.reserve(subpop.size());
    for (int_fast32_t i = 0; i < ssize(litter_sizes); ++i) {
        auto& mother = females[cast_u(i)];
        auto rec = litter_sizes[cast_u(i)];
        for (decltype(rec) j = 0; j < rec; ++j) {
            const auto& father = males[cast_u(mate_distr(*engine_))];
            juveniles.emplace_back(std::make_shared<Individual>(father, mother, year_));
        }
    }
}

std::vector<int_fast32_t>
Population::litter_sizes_logistic(const std::vector<ShPtrIndividual>& females, const int_fast32_t popsize) {
    const auto N = static_cast<double>(popsize);
    if (N > params_.carrying_capacity) return {};
    const double rec_rate = params_.recruitment * (1.0 - N / params_.carrying_capacity);
    std::vector<int_fast32_t> res;
    res.reserve(females.size());
    for (const auto& mother: females) {
        const double rec_mean = rec_rate * WEIGHT_FOR_AGE_[cast_u(mother->age(year_))];
        res.push_back(rnbinom<int_fast32_t>(params_.overdispersion, rec_mean, *engine_));
    }
    return res;
}

std::vector<double> Population::weights(const std::vector<ShPtrIndividual>& individuals) const {
    std::vector<double> res;
    res.reserve(individuals.size());
    for (const auto& p: individuals) {
        res.push_back(WEIGHT_FOR_AGE_[cast_u(p->age(year_))]);
    }
    return res;
}

std::vector<double> Population::survive_juvenile(
  std::vector<int_fast32_t>& litter_sizes,
  std::array<std::vector<int_fast32_t>, 4>& demography_y,
  const int_fast32_t sample_size) {
    std::vector<double> d0(4u);
    for (int_fast32_t q = 0; q < 4; ++q) {
        d0[cast_u(q)] = death_rate(0, q);
    }
    std::vector<double> dead3(sample_size ? litter_sizes.size() : 0u);
    for (int_fast32_t q = 0; q < 4; ++q) {
        auto& demography_q0 = demography_y[cast_u(q)][0u];
        for (int_fast32_t i = 0; i < ssize(litter_sizes); ++i) {
            auto& n = litter_sizes[cast_u(i)];
            demography_q0 += n;
            std::binomial_distribution<int_fast32_t> distr(n, d0[cast_u(q)]);
            const auto death = distr(*engine_);
            n -= death;
            if (sample_size > 0 && q == 3) dead3[cast_u(i)] = static_cast<double>(death);
        }
    }
    return dead3;
}

int_fast32_t Population::sample_size_juvenile(int_fast32_t loc) const {
    if (year_ > (params_.years - params_.last)) {
        return params_.sample_size_juvenile[cast_u(loc)];
    }
    return 0;
}

void Population::survive(const int_fast32_t season) {
    int_fast32_t total_n{0};
    const bool is_sampling_time = season == 3 && year_ > (params_.years - params_.last);
    for (int_fast32_t loc = 0; loc < ssize(subpopulations_); ++loc) {
      auto ss_loc = is_sampling_time ? params_.sample_size_adult[cast_u(loc)] : 0;
      std::vector<ShPtrIndividual> carcasses;
      carcasses.reserve(cast_u(ss_loc));
      auto& subpop = subpopulations_[cast_u(loc)];
      for (const auto sex: {Sex::F, Sex::M}) {
        auto& individuals = subpop[sex];
        auto sub_n = ssize(individuals);
        for (decltype(sub_n) i=0; i<sub_n; ++i) {
            auto& p = individuals[cast_u(i)];
            if (wtl::generate_canonical(*engine_) < death_rate(p->age(year_), season)) {
                if (ss_loc > 0) {
                    carcasses.emplace_back(std::move(p));
                }
                p = std::move(individuals.back());
                individuals.pop_back();
                --sub_n;
                --i;
            }
        }
        total_n += sub_n;
      }
      if (ss_loc == 0) continue;
      if (ssize(carcasses) < ss_loc) {
        std::cerr << "WARNING:Population::survive():"
                  << "y" << year_ << "-l" << loc
                  << ": " << ss_loc - ssize(carcasses) << " fewer samples\n";
        ss_loc = static_cast<int_fast32_t>(ssize(carcasses));
      }
      std::sample(carcasses.begin(), carcasses.end(),
                  std::back_inserter(subpop.samples[year_]), ss_loc, *engine_);
    }
    if (total_n == 0) {
        std::ostringstream oss{"runtime_error:", std::ios_base::ate};
        oss << "extinction in " << year_ << "-q" << season;
        throw std::runtime_error(oss.str());
    }
}

double Population::death_rate(const int_fast32_t age, const int_fast32_t season) const {
    const auto q_age = 4 * age + season;
    const auto m = NATURAL_MORTALITY_[cast_u(q_age)];
    const auto f = FISHING_MORTALITY_[cast_u(q_age)] * FISHING_COEF_[cast_u(year_)];
    return 1.0 - std::exp(-m - f);
}

void Population::migrate() {
    for (const auto sex: {Sex::F, Sex::M}) {
        // sizes before migration is needed for loc > 0
        std::vector<ptrdiff_t> subpop_sizes;
        subpop_sizes.reserve(subpopulations_.size());
        for (const auto& subpop: subpopulations_) {
            subpop_sizes.push_back(ssize(subpop[sex]));
        }
        for (int_fast32_t loc=0; loc<ssize(subpopulations_); ++loc) {
            auto& individuals = subpopulations_[cast_u(loc)][sex];
            auto n = subpop_sizes[cast_u(loc)];
            for (decltype(n) i=0; i<n;) {
                auto& p = individuals[cast_u(i)];
                int_fast32_t dst = destination(p->age(year_), loc);
                if (dst == loc) {++i; continue;}
                subpopulations_[cast_u(dst)][sex].emplace_back(std::move(p));
                if (i < --n)
                    p = std::move(individuals[cast_u(n)]);
                if (individuals.back())
                    individuals[cast_u(n)] = std::move(individuals.back());
                individuals.pop_back();
            }
        }
    }
    for (int_fast32_t loc=0; loc<ssize(subpopulations_); ++loc) {
        auto& juveniles = subpopulations_[cast_u(loc)].juveniles;
        for (auto& p: juveniles) {
            constexpr URBG::result_type half_max = URBG::max() >> 1u;
            const auto sex = engine_->operator()() < half_max ? Sex::F : Sex::M;
            const auto dst = destination(0, loc);
            subpopulations_[cast_u(dst)][sex].emplace_back(std::move(p));
        }
        juveniles.clear();
    }
}

int_fast32_t Population::destination(int_fast32_t age, int_fast32_t loc) {
    auto& [dest, dist] = MIGRATION_DESTINATION_[cast_u(age)][cast_u(loc)];
    return (dest == std::numeric_limits<int_fast32_t>::max()) ? dist(*engine_) : dest;
}

std::ostream& Population::write_sample_family(std::ostream& ost) const {
    Individual::write_trace_back_header(ost);
    // unordered_set suffices to remove duplicates, but address is not reproducible.
    std::unordered_map<const Individual*, int_fast32_t> ids;
    ids.emplace(nullptr, 0);
    for (int_fast32_t loc=0; loc<ssize(subpopulations_); ++loc) {
        for (const auto& [_year, samples]: subpopulations_[cast_u(loc)].samples) {
            for (const auto& p: samples) {
                ids.emplace(p.get(), ssize(ids));
            }
        }
    }
    for (int_fast32_t loc=0; loc<ssize(subpopulations_); ++loc) {
        for (const auto& [year, samples]: subpopulations_[cast_u(loc)].samples) {
            for (const auto& p: samples) {
                p->trace_back(ost, ids, loc, year);
            }
        }
    }
    return ost;
}

void Population::record_demography(const int_fast32_t season) {
    for (auto& subpop: subpopulations_) {
        auto& counter = subpop.demography[cast_u(year_)][cast_u(season)];
        for (const auto& inds: subpop.adults) {
            for (const auto& p: inds) {
                ++counter[cast_u(p->age(year_))];
            }
        }
    }
}

std::ostream& Population::write_demography(std::ostream& ost) const {
    ost << "year\tseason\tlocation\tage\tcount\n";
    for (int_fast32_t loc = 0; loc < ssize(subpopulations_); ++loc) {
      const auto& subpop = subpopulations_[cast_u(loc)];
      for (int_fast32_t year = 0; year < ssize(subpop.demography); ++year) {
        const auto& demography_year = subpop.demography[cast_u(year)];
        for (int_fast32_t season = 0; season < ssize(demography_year); ++season) {
            const auto& structure = demography_year[cast_u(season)];
            for (int_fast32_t age=0; age<ssize(structure); ++age) {
                if (structure[cast_u(age)] == 0) continue;
                ost << year << "\t" << season << "\t"
                    << loc << "\t"
                    << age << "\t"
                    << structure[cast_u(age)] << "\n";
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
    copy_elongate(params_.natural_mortality, NATURAL_MORTALITY_, 4 * MAX_AGE);
    copy_elongate(params_.fishing_mortality, FISHING_MORTALITY_, 4 * MAX_AGE);
    NATURAL_MORTALITY_.back() = 1e9;
    const auto fc_size = ssize(params_.fishing_coef);
    const auto offset = ramp(fc_size - params_.years);
    FISHING_COEF_.assign(cast_u(params_.years), 1.0);
    std::copy_backward(params_.fishing_coef.begin() + offset, params_.fishing_coef.end(),
                       FISHING_COEF_.end());
}

void Population::init_weight() {
    WEIGHT_FOR_AGE_.reserve(MAX_AGE);
    WEIGHT_FOR_AGE_.resize(params_.weight_for_age.size() / 4u);
    for (int_fast32_t year=0; year<ssize(WEIGHT_FOR_AGE_); ++year) {
        WEIGHT_FOR_AGE_[cast_u(year)] = params_.weight_for_age[cast_u(4 * year)];
    }
    elongate(WEIGHT_FOR_AGE_, MAX_AGE);
}

bool Population::is_ready() const {
    return (
      ssize(FISHING_COEF_) >= params_.years
      && ssize(NATURAL_MORTALITY_) >= 4 * MAX_AGE
      && ssize(FISHING_MORTALITY_) >= 4 * MAX_AGE
      && ssize(WEIGHT_FOR_AGE_) >= MAX_AGE
      && ssize(MIGRATION_DESTINATION_) >= MAX_AGE
    );
}

} // namespace pbf
