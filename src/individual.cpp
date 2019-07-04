/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"
#include "config.hpp"

#include <wtl/random.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <clippson/json.hpp>

#include <type_traits>

namespace pbf {

Individual::param_type Individual::PARAM_;
IndividualJson Individual::JSON_;

static_assert(!std::is_default_constructible<Individual>{}, "");
static_assert(std::is_nothrow_copy_constructible<Individual>{}, "");
static_assert(std::is_nothrow_move_constructible<Individual>{}, "");

bool Individual::is_dead(const int_fast32_t year, URBG& engine) const {
    return (wtl::generate_canonical(engine) < JSON_.DEATH_RATE[year - birth_year_]);
}

//! Translate parameter `mean` to `prob`
template <class T> inline wtl::negative_binomial_distribution<T>
nbinom_distribution(double k, double mu) {
    const double prob = k / (mu + k);
    return wtl::negative_binomial_distribution<T>(k, prob);
}

uint_fast32_t Individual::recruitment(const int_fast32_t year, const double density_effect, URBG& engine) const noexcept {
    if (density_effect < 0.0) return 0u;
    const double mean = density_effect * param().RECRUITMENT_COEF * weight(year);
    const double k = param().NEGATIVE_BINOM_K;
    if (k > 0.0) {
        return nbinom_distribution<uint_fast32_t>(k, mean)(engine);
    } else {
        return std::poisson_distribution<uint_fast32_t>(mean)(engine);
    }
}

uint_fast32_t Individual::migrate(const uint_fast32_t loc, const int_fast32_t year, URBG& engine) {
    return JSON_.MIGRATION_DISTRIBUTIONS.at(year - birth_year_)[loc](engine);
}

void Individual::trace_back(std::ostream& ost, std::map<const Individual*, size_t>* ids, uint32_t loc, int_fast32_t year) const {
    if (!ids->emplace(this, ids->size()).second && (year == 0)) return;
    if (father_) father_->trace_back(ost, ids, loc, 0);
    if (mother_) mother_->trace_back(ost, ids, loc, 0);
    write(ost, *ids);
    if (year > 0) {
        ost << "\t" << loc << "\t" << year << "\n";
    } else {
        ost << "\t\t\n";
    }
}

std::vector<std::string> Individual::names() {
    return {"id", "father_id", "mother_id", "birth_year"};
}

std::ostream& Individual::write(std::ostream& ost) const {
    return ost << this << "\t"
               << father_ << "\t"
               << mother_ << "\t"
               << birth_year_;
}

std::ostream& Individual::write(std::ostream& ost, const std::map<const Individual*, size_t>& ids) const {
    return ost << ids.at(this) << "\t"
               << ids.at(father_.get()) << "\t"
               << ids.at(mother_.get()) << "\t"
               << birth_year_;
}

//! Shortcut of Individual::write
std::ostream& operator<<(std::ostream& ost, const Individual& x) {
    return x.write(ost);
}

std::string Individual::default_json() {
    return default_values;
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
    constexpr int_fast32_t max_age = 80;
    constexpr int_fast32_t max_qage = 4 * (max_age + 1);
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
    DEATH_RATE.reserve(max_age);
    DEATH_RATE.resize(NATURAL_MORTALITY.size() / 4u);
    for (size_t year=0; year<DEATH_RATE.size(); ++year) {
        double z = 0.0;
        for (size_t q = 4u * year, q_end = q + 4u; q<q_end; ++q) {
            z += NATURAL_MORTALITY[q];
            z += FISHING_MORTALITY[q];
        }
        DEATH_RATE[year] = 1.0 - std::exp(-z);
    }
    elongate(&DEATH_RATE, max_age);
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
    URBG engine(std::random_device{}());
    auto dist = nbinom_distribution<int>(k, mu);
    std::vector<int> v(n);
    for (int& x: v) {
        x = dist(engine);
    }
    return v;
}

} // namespace pbf
