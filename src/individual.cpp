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
#include <boost/program_options.hpp>

namespace pbt {

double Individual::RECRUITMENT_COEF_ = 1.0;
std::vector<double> Individual::NATURAL_MORTALITY_;
std::vector<double> Individual::FISHING_MORTALITY_;
std::vector<double> Individual::SURVIVAL_RATE_;
std::vector<double> Individual::WEIGHT_FOR_AGE_;
std::vector<std::vector<std::vector<double>>> Individual::MIGRATION_MATRICES_;
uint_fast32_t Individual::LAST_ID_ = 0;

//! discrete distributions for migration
static std::vector<std::vector<std::discrete_distribution<uint_fast32_t>>> MIGRATION_DISTRIBUTIONS;

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option | Symbol         | Variable
    ------------------- | -------------- | -------------------------------
    `-r,--recruitment`  | -              | Individual::RECRUITMENT_COEF_
*/
boost::program_options::options_description Individual::options_desc() {
    namespace po = boost::program_options;
    po::options_description desc{"Individual"};
    desc.add_options()
        ("recruitment,r", po::value(&RECRUITMENT_COEF_)->default_value(RECRUITMENT_COEF_))
    ;
    return desc;
}

void Individual::set_default_values() {HERE;
    if (!NATURAL_MORTALITY_.empty()) return;
    std::istringstream iss(default_values);
    read_json(iss);
}

void Individual::set_dependent_static() {HERE;
    MIGRATION_DISTRIBUTIONS.clear();
    MIGRATION_DISTRIBUTIONS.reserve(MAX_AGE_);
    for (const auto& matrix: MIGRATION_MATRICES_) {
        decltype(MIGRATION_DISTRIBUTIONS)::value_type dists;
        dists.reserve(matrix.size());
        for (const auto& row: matrix) {
            dists.emplace_back(row.begin(), row.end());
        }
        MIGRATION_DISTRIBUTIONS.emplace_back(std::move(dists));
    }
    for (size_t i=MIGRATION_DISTRIBUTIONS.size(); i<MAX_AGE_; ++i) {
        MIGRATION_DISTRIBUTIONS.emplace_back(MIGRATION_DISTRIBUTIONS.back());
    }
    SURVIVAL_RATE_.resize(NATURAL_MORTALITY_.size());
    std::transform(NATURAL_MORTALITY_.begin(), NATURAL_MORTALITY_.end(),
                   FISHING_MORTALITY_.begin(), SURVIVAL_RATE_.begin(),
                   [](double n, double f) {
                       return std::exp(-n - f);
                   });
}

bool Individual::has_survived(const uint_fast32_t year, const uint_fast32_t quarter, URBG& engine) const {
    const auto age = year - birth_year_;
    return (wtl::generate_canonical(engine) < SURVIVAL_RATE_[4U * age + quarter]);
}

uint_fast32_t Individual::recruitment(const uint_fast32_t year, URBG& engine) const {
    const double mean = RECRUITMENT_COEF_ * weight(year);
    std::poisson_distribution<uint_fast32_t> poisson(mean);
    return poisson(engine);
}

void Individual::migrate(const uint_fast32_t year, URBG& engine) {
    location_ = MIGRATION_DISTRIBUTIONS[year - birth_year_][location_](engine);
}

std::vector<std::string> Individual::names() {
    return {"id", "father_id", "mother_id", "birth_year", "location"};
}

std::ostream& Individual::write(std::ostream& ost) const {
    return write_ids(ost) << "\t"
               << birth_year_ << "\t"
               << location_;
}

std::ostream& Individual::write_ids(std::ostream& ost) const {
    return ost << id_ << "\t"
               << (father_ ? father_->id() : 0u) << "\t"
               << (mother_ ? mother_->id() : 0u);
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
