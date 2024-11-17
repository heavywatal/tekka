/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"
#include "config.hpp"

#include <wtl/debug.hpp>
#include <clippson/json.hpp>

#include <sstream>
#include <type_traits>

namespace pbf {

static_assert(std::is_nothrow_default_constructible_v<Individual>);
static_assert(!std::is_copy_constructible_v<Individual>);
static_assert(!std::is_move_constructible_v<Individual>);
static_assert(std::is_nothrow_destructible_v<Individual>);

namespace {

template <class T> inline
void elongate(std::vector<T>& v, size_t n) noexcept {
    for (size_t i=v.size(); i<n; ++i) {
        v.emplace_back(v.back());
    }
}

inline int_fast32_t get_dest(const std::vector<double>& v) {
    uint_fast32_t idx = 0;
    uint_fast32_t num_positive = 0u;
    for (uint_fast32_t i=0u; i<v.size(); ++i) {
        if (v[i] > 0.0) {
            ++num_positive;
            idx = i;
        }
    }
    return num_positive == 1u ? idx : -1;
}

inline uint_fast32_t sub_sat(const uint_fast32_t x, const uint_fast32_t y) {
    return (x > y) ? x - y : uint_fast32_t{};
}

} // anonymous namespace

void Individual::trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>& ids,
                            uint_fast32_t loc, int_fast32_t year) const {
    if (!ids.emplace(this, static_cast<uint_fast32_t>(ids.size())).second && (year == 0)) return;
    if (father_) father_->trace_back(ost, ids, loc, 0);
    if (mother_) mother_->trace_back(ost, ids, loc, 0);
    write(ost, ids);
    if (year > 0) {
        ost << "\t" << loc << "\t" << year << "\n";
    } else {
        ost << "\t\t\n";
    }
}

std::ostream& Individual::write_trace_back_header(std::ostream& ost) {
   return write_names(ost) << "\tlocation\tcapture_year\n";
}

std::ostream& Individual::write_names(std::ostream& ost) {
    return ost << "id\t"
               << "father_id\t"
               << "mother_id\t"
               << "birth_year";
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

void Individual::set_static_migration() {
    using dist_t = PairDestDist::second_type;
    MIGRATION_DESTINATION_.clear();
    MIGRATION_DESTINATION_.reserve(MAX_AGE);
    for (const auto& matrix: JSON.migration_matrices) {
        std::vector<PairDestDist> pairs;
        pairs.reserve(matrix.size());
        for (const auto& row: matrix) {
            pairs.emplace_back(get_dest(row), dist_t(row.begin(), row.end()));
        }
        MIGRATION_DESTINATION_.emplace_back(std::move(pairs));
    }
    elongate(MIGRATION_DESTINATION_, MAX_AGE);
}

void Individual::set_static_mortality() {
    elongate(JSON.natural_mortality, 4u * MAX_AGE);
    elongate(JSON.fishing_mortality, 4u * MAX_AGE);
    JSON.natural_mortality.back() = 1e9;
}

void Individual::set_static_weight() {
    WEIGHT_FOR_AGE_.reserve(MAX_AGE);
    WEIGHT_FOR_AGE_.resize(JSON.weight_for_age.size() / 4u);
    for (size_t year=0; year<WEIGHT_FOR_AGE_.size(); ++year) {
        WEIGHT_FOR_AGE_[year] = JSON.weight_for_age[4u * year];
    }
    elongate(WEIGHT_FOR_AGE_, MAX_AGE);
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
    if (!is_ready(years)) {
        throw std::logic_error("logic_error:!Individual::is_ready()");
    }
}

IndividualJson::IndividualJson() {
    std::istringstream iss(default_values);
    read(iss);
}

void IndividualJson::read(std::istream& ist) {
    nlohmann::json obj;
    ist >> obj;
    obj.at("natural_mortality").get_to(natural_mortality);
    obj.at("fishing_mortality").get_to(fishing_mortality);
    fishing_coef = obj.value("fishing_coef", fishing_coef);
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

} // namespace pbf
