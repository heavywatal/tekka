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

static_assert(!std::is_default_constructible<Individual>{}, "");
static_assert(std::is_nothrow_copy_constructible<Individual>{}, "");
static_assert(std::is_nothrow_move_constructible<Individual>{}, "");

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
    return MIGRATION_DISTRIBUTIONS_[year - birth_year_][loc](engine);
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

void Individual::set_static_migration() {
    MIGRATION_DISTRIBUTIONS_.clear();
    MIGRATION_DISTRIBUTIONS_.reserve(MAX_AGE);
    for (const auto& matrix: JSON.migration_matrices) {
        decltype(MIGRATION_DISTRIBUTIONS_)::value_type distributions;
        distributions.reserve(matrix.size());
        for (const auto& row: matrix) {
            distributions.emplace_back(make_dist(row));
        }
        MIGRATION_DISTRIBUTIONS_.emplace_back(std::move(distributions));
    }
    elongate(&MIGRATION_DISTRIBUTIONS_, MAX_AGE);
}

void Individual::set_static_mortality() {
    const auto size = JSON.natural_mortality.size() / 4u;
    NATURAL_MORTALITY_.reserve(MAX_AGE);
    NATURAL_MORTALITY_.assign(size, 0.0);
    FISHING_MORTALITY_.reserve(MAX_AGE);
    FISHING_MORTALITY_.assign(size, 0.0);
    for (size_t year=0; year<size; ++year) {
        for (size_t q = 0; q < 4u; ++q) {
            const auto q_i = 4u * year + q;
            NATURAL_MORTALITY_[year] += JSON.natural_mortality[q_i];
            FISHING_MORTALITY_[year] += JSON.fishing_mortality[q_i];
        }
    }
    elongate(&NATURAL_MORTALITY_, MAX_AGE);
    elongate(&FISHING_MORTALITY_, MAX_AGE);
    NATURAL_MORTALITY_.back() = 1e9;
}

void Individual::set_static_weight() {
    WEIGHT_FOR_AGE_.reserve(MAX_AGE);
    WEIGHT_FOR_AGE_.resize(JSON.weight_for_age.size() / 4u);
    for (size_t year=0; year<WEIGHT_FOR_AGE_.size(); ++year) {
        WEIGHT_FOR_AGE_[year] = JSON.weight_for_age[4u * year];
    }
    elongate(&WEIGHT_FOR_AGE_, MAX_AGE);
}

inline uint_fast32_t sub_sat(const uint_fast32_t x, const uint_fast32_t y) {
    return (x > y) ? x - y : uint_fast32_t{};
}

void Individual::set_dependent_static(const uint_fast32_t years) {
    set_static_migration();
    set_static_mortality();
    set_static_weight();
    const auto fc_size = static_cast<uint_fast32_t>(JSON.fishing_coef.size());
    const auto offset = sub_sat(fc_size, years);
    FISHING_COEF_.assign(years, 1.0);
    std::copy_backward(JSON.fishing_coef.begin() + offset, JSON.fishing_coef.end(),
                       FISHING_COEF_.end());
}

//! @cond
IndividualJson::IndividualJson() {
    std::istringstream iss(default_values);
    read(iss);
}

void IndividualJson::read(std::istream& ist) {
    nlohmann::json obj;
    ist >> obj;
    obj.at("natural_mortality").get_to(natural_mortality);
    obj.at("fishing_mortality").get_to(fishing_mortality);
    fishing_mortality = obj.value("fishing_coef", fishing_mortality);
    obj.at("weight_for_age").get_to(weight_for_age);
    obj.at("migration_matrices").get_to(migration_matrices);
}

std::ostream& IndividualJson::write(std::ostream& ost) const {
    nlohmann::json obj;
    obj["natural_mortality"] = natural_mortality;
    obj["fishing_mortality"] = fishing_mortality;
    obj["fishing_coef"] = fishing_coef;
    obj["weight_for_age"] = weight_for_age;
    obj["migration_matrices"] = migration_matrices;
    return ost << obj;
}
//! @endcond

} // namespace pbf
