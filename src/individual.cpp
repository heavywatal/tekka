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
    return JSON_.MIGRATION_DISTRIBUTIONS[year - birth_year_][loc](engine);
}

void Individual::trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>* ids,
                            uint_fast32_t loc, int_fast32_t year) const {
    if (!ids->emplace(this, static_cast<uint_fast32_t>(ids->size())).second && (year == 0)) return;
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

std::ostream& Individual::write(std::ostream& ost, const std::unordered_map<const Individual*, uint_fast32_t>& ids) const {
    return ost << ids.at(this) << "\t"
               << ids.at(father_.get()) << "\t"
               << ids.at(mother_.get()) << "\t"
               << birth_year_;
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

inline std::function<uint_fast32_t(URBG&)> make_dist(const std::vector<double>& v) {
    uint_fast32_t idx = 0u;
    unsigned num_positive = 0u;
    for (uint_fast32_t i=0u; i<v.size(); ++i) {
        if (v[i] > 0.0) {
            ++num_positive;
            idx = i;
        }
    }
    if (num_positive == 1u) {
        return [idx](URBG&){return idx;};
    } else {
        return std::discrete_distribution<uint_fast32_t>(v.begin(), v.end());
    }
}

void IndividualJson::set_dependent_static() {
    constexpr int_fast32_t max_age = 80;
    MIGRATION_DISTRIBUTIONS.clear();
    MIGRATION_DISTRIBUTIONS.reserve(max_age);
    for (const auto& matrix: MIGRATION_MATRICES) {
        decltype(MIGRATION_DISTRIBUTIONS)::value_type dists;
        dists.reserve(matrix.size());
        for (const auto& row: matrix) {
            dists.emplace_back(make_dist(row));
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
    DEATH_RATE.back() = 1.0;
    WEIGHT_FOR_YEAR_AGE.reserve(max_age);
    WEIGHT_FOR_YEAR_AGE.resize(WEIGHT_FOR_AGE.size() / 4u);
    for (size_t year=0; year<WEIGHT_FOR_YEAR_AGE.size(); ++year) {
        WEIGHT_FOR_YEAR_AGE[year] = WEIGHT_FOR_AGE[4u * year];
    }
    elongate(&WEIGHT_FOR_YEAR_AGE, max_age);
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

} // namespace pbf
