/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"
#include "random.hpp"
#include "config.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <clippson/json.hpp>

#include <type_traits>

namespace pbf {

Individual::param_type Individual::PARAM_;
IndividualJson Individual::JSON_;

static_assert(std::is_nothrow_default_constructible<Individual>{}, "");
static_assert(std::is_nothrow_copy_constructible<Individual>{}, "");
static_assert(std::is_nothrow_move_constructible<Individual>{}, "");

bool Individual::has_survived(const uint_fast32_t year, const uint_fast32_t quarter, URBG& engine) const {
    uint_fast32_t qage = 4u * (year - birth_year_) + quarter;
    return (wtl::generate_canonical(engine) < JSON_.SURVIVAL_RATE.at(qage));
}

//! Translate parameter `mean` to `prob`
template <class T> inline wtl::negative_binomial_distribution<T>
nbinom_distribution(double k, double mu) {
    const double prob = k / (mu + k);
    return wtl::negative_binomial_distribution<T>(k, prob);
}

uint_fast32_t Individual::recruitment(const uint_fast32_t year, const size_t popsize, URBG& engine) const noexcept {
    const double coef = (1.0 - popsize / param().CARRYING_CAPACITY);
    if (coef < 0.0) return 0u;
    const double mean = coef * param().RECRUITMENT_COEF * weight(year);
    const double k = param().NEGATIVE_BINOM_K;
    if (k > 0.0) {
        return nbinom_distribution<uint_fast32_t>(k, mean)(engine);
    } else {
        return std::poisson_distribution<uint_fast32_t>(mean)(engine);
    }
}

void Individual::migrate(const uint_fast32_t year, URBG& engine) {
    location_ = JSON_.MIGRATION_DISTRIBUTIONS.at(year - birth_year_)[location_](engine);
}

void Individual::trace_back(std::ostream& ost, std::map<const Individual*, size_t>* ids, uint_fast32_t year) const {
    if (!ids->emplace(this, ids->size()).second) return;
    if (father_) {
        father_->trace_back(ost, ids);
        mother_->trace_back(ost, ids);
    }
    write(ost, *ids) << "\t";
    if (year > 0u) {
        ost << year;
    }
    ost << "\n";
}

std::vector<std::string> Individual::names() {
    return {"id", "father_id", "mother_id", "birth_year", "location"};
}

std::ostream& Individual::write(std::ostream& ost) const {
    return ost << this << "\t"
               << father_ << "\t"
               << mother_ << "\t"
               << birth_year_ << "\t"
               << location_;
}

std::ostream& Individual::write(std::ostream& ost, const std::map<const Individual*, size_t>& ids) const {
    return ost << ids.at(this) << "\t"
               << ids.at(father_.get()) << "\t"
               << ids.at(mother_.get()) << "\t"
               << birth_year_ << "\t"
               << location_;
}

//! Shortcut of Individual::write
std::ostream& operator<<(std::ostream& ost, const Individual& x) {
    return x.write(ost);
}

//! @cond
IndividualJson::IndividualJson() {
    std::istringstream iss(default_values);
    read(iss);
}

template <class T> inline
void elongate(std::vector<T>* v, size_t n) noexcept {
    for (size_t i=v->size(); i<n; ++i) {
        v->emplace_back(v->back());
    }
}

void IndividualJson::set_dependent_static() {
    constexpr uint_fast32_t max_age = 80u;
    constexpr uint_fast32_t max_qage = 4u * (max_age + 1u);
    MIGRATION_DISTRIBUTIONS.clear();
    MIGRATION_DISTRIBUTIONS.reserve(max_age);
    for (const auto& matrix: MIGRATION_MATRICES) {
        decltype(MIGRATION_DISTRIBUTIONS)::value_type dists;
        dists.reserve(matrix.size());
        for (const auto& row: matrix) {
            dists.emplace_back(row.begin(), row.end());
        }
        MIGRATION_DISTRIBUTIONS.emplace_back(std::move(dists));
    }
    elongate(&MIGRATION_DISTRIBUTIONS, max_age);
    SURVIVAL_RATE.reserve(max_qage);
    SURVIVAL_RATE.resize(NATURAL_MORTALITY.size());
    std::transform(NATURAL_MORTALITY.begin(), NATURAL_MORTALITY.end(),
                   FISHING_MORTALITY.begin(), SURVIVAL_RATE.begin(),
                   [](double n, double f) {
                       return std::exp(-n - f);
                   });
    elongate(&SURVIVAL_RATE, max_qage);
    elongate(&WEIGHT_FOR_AGE, max_qage);
}

void IndividualJson::read(std::istream& ist) {
    nlohmann::json obj;
    ist >> obj;
    NATURAL_MORTALITY = obj.at("natural_mortality").get<decltype(NATURAL_MORTALITY)>();
    FISHING_MORTALITY = obj.at("fishing_mortality").get<decltype(FISHING_MORTALITY)>();
    WEIGHT_FOR_AGE = obj.at("weight_for_age").get<decltype(WEIGHT_FOR_AGE)>();
    MIGRATION_MATRICES = obj.at("migration_matrices").get<decltype(MIGRATION_MATRICES)>();
    set_dependent_static();
}

void IndividualJson::write(std::ostream& ost) const {
    nlohmann::json obj;
    obj["natural_mortality"] = NATURAL_MORTALITY;
    obj["fishing_mortality"] = FISHING_MORTALITY;
    obj["weight_for_age"] = WEIGHT_FOR_AGE;
    obj["migration_matrices"] = MIGRATION_MATRICES;
    ost << obj;
}
//! @endcond

std::vector<int> Individual::rnbinom(const int n, const double k, const double mu) {
    auto dist = nbinom_distribution<int>(k, mu);
    std::vector<int> v(n);
    for (int& x: v) {
        x = dist(engine64());
    }
    return v;
}

} // namespace pbf
