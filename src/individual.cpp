/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"
#include "config.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/random.hpp>
#include <sfmt.hpp>
#include <nlohmann/json.hpp>

#include <type_traits>

namespace pbt {

Individual::param_type Individual::PARAM_;
std::vector<double> Individual::NATURAL_MORTALITY_;
std::vector<double> Individual::FISHING_MORTALITY_;
std::vector<double> Individual::SURVIVAL_RATE_;
std::vector<double> Individual::WEIGHT_FOR_AGE_;
std::vector<std::vector<std::vector<double>>> Individual::MIGRATION_MATRICES_;

//! discrete distributions for migration
static std::vector<std::vector<std::discrete_distribution<uint_fast32_t>>> MIGRATION_DISTRIBUTIONS;

static_assert(std::is_nothrow_default_constructible<Individual>{}, "");
static_assert(std::is_nothrow_copy_constructible<Individual>{}, "");
static_assert(std::is_nothrow_move_constructible<Individual>{}, "");

void Individual::param(const param_type& p) {
    PARAM_ = p;
}

void Individual::set_default_values() {HERE;
    if (!NATURAL_MORTALITY_.empty()) return;
    std::istringstream iss(default_values);
    read_json(iss);
}

//! append the last element until v->size() reaches n
template <class T> inline
void elongate(std::vector<T>* v, size_t n) noexcept {
    for (size_t i=v->size(); i<n; ++i) {
        v->emplace_back(v->back());
    }
}

void Individual::set_dependent_static() {HERE;
    constexpr uint_fast32_t max_age = 80u;
    constexpr uint_fast32_t max_qage = 4u * (max_age + 1u);
    MIGRATION_DISTRIBUTIONS.clear();
    MIGRATION_DISTRIBUTIONS.reserve(max_age);
    for (const auto& matrix: MIGRATION_MATRICES_) {
        decltype(MIGRATION_DISTRIBUTIONS)::value_type dists;
        dists.reserve(matrix.size());
        for (const auto& row: matrix) {
            dists.emplace_back(row.begin(), row.end());
        }
        MIGRATION_DISTRIBUTIONS.emplace_back(std::move(dists));
    }
    elongate(&MIGRATION_DISTRIBUTIONS, max_age);
    SURVIVAL_RATE_.reserve(max_qage);
    SURVIVAL_RATE_.resize(NATURAL_MORTALITY_.size());
    std::transform(NATURAL_MORTALITY_.begin(), NATURAL_MORTALITY_.end(),
                   FISHING_MORTALITY_.begin(), SURVIVAL_RATE_.begin(),
                   [](double n, double f) {
                       return std::exp(-n - f);
                   });
    elongate(&SURVIVAL_RATE_, max_qage);
    elongate(&WEIGHT_FOR_AGE_, max_qage);
}

bool Individual::has_survived(const uint_fast32_t year, const uint_fast32_t quarter, URBG& engine) const {
    uint_fast32_t qage = 4u * (year - birth_year_) + quarter;
    return (wtl::generate_canonical(engine) < SURVIVAL_RATE_.at(qage));
}

uint_fast32_t Individual::recruitment(const uint_fast32_t year, URBG& engine) const noexcept {
    const double mean = param().RECRUITMENT_COEF * weight(year);
    const double prob = param().NEGATIVE_BINOM_K / (mean + param().NEGATIVE_BINOM_K);
    if (prob < 1.0) {
        return wtl::negative_binomial_distribution<uint_fast32_t>(param().NEGATIVE_BINOM_K, prob)(engine);
    } else {
        return std::poisson_distribution<uint_fast32_t>(mean)(engine);
    }
}

void Individual::migrate(const uint_fast32_t year, URBG& engine) {
    location_ = MIGRATION_DISTRIBUTIONS.at(year - birth_year_)[location_](engine);
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

//! shortcut Individual::write(ost)
std::ostream& operator<<(std::ostream& ost, const Individual& x) {
    return x.write(ost);
}

void Individual::read_json(std::istream& ist) {
    nlohmann::json obj;
    ist >> obj;
    NATURAL_MORTALITY_ = obj.at("natural_mortality").get<decltype(NATURAL_MORTALITY_)>();
    FISHING_MORTALITY_ = obj.at("fishing_mortality").get<decltype(FISHING_MORTALITY_)>();
    WEIGHT_FOR_AGE_ = obj.at("weight_for_age").get<decltype(WEIGHT_FOR_AGE_)>();
    MIGRATION_MATRICES_ = obj.at("migration_matrices").get<decltype(MIGRATION_MATRICES_)>();
    set_dependent_static();
}

void Individual::write_json(std::ostream& ost) {
    nlohmann::json obj;
    obj["natural_mortality"] = NATURAL_MORTALITY_;
    obj["fishing_mortality"] = FISHING_MORTALITY_;
    obj["weight_for_age"] = WEIGHT_FOR_AGE_;
    obj["migration_matrices"] = MIGRATION_MATRICES_;
    ost << obj;
}

} // namespace pbt
