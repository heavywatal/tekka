/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/random.hpp>

namespace pbt {

double Individual::RECRUITMENT_COEF_ = 0.2;
std::vector<double> Individual::NATURAL_MORTALITY_;
std::vector<double> Individual::FISHING_MORTALITY_;
std::vector<double> Individual::WEIGHT_FOR_AGE_;
std::vector<std::vector<std::vector<double>>> Individual::MIGRATION_MATRICES_;
std::vector<std::vector<std::discrete_distribution<uint_fast32_t>>> Individual::MIGRATION_DISTRIBUTIONS_;
uint_fast32_t Individual::LAST_ID_ = 0;

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

void Individual::set_dependent_static() {HERE;
    for (const auto& matrix: MIGRATION_MATRICES_) {
        decltype(MIGRATION_DISTRIBUTIONS_)::value_type dists;
        dists.reserve(matrix.size());
        for (const auto& row: matrix) {
            dists.emplace_back(row.begin(), row.end());
        }
        MIGRATION_DISTRIBUTIONS_.push_back(std::move(dists));
    }
    for (size_t i=MIGRATION_DISTRIBUTIONS_.size(); i<MAX_AGE_; ++i) {
        MIGRATION_DISTRIBUTIONS_.push_back(MIGRATION_DISTRIBUTIONS_.back());
    }
}

void Individual::from_json(const json::json& obj) {HERE;
    NATURAL_MORTALITY_ = obj.at("natural_mortality").get<decltype(NATURAL_MORTALITY_)>();
    FISHING_MORTALITY_ = obj.at("fishing_mortality").get<decltype(FISHING_MORTALITY_)>();
    WEIGHT_FOR_AGE_ = obj.at("weight_for_age").get<decltype(WEIGHT_FOR_AGE_)>();
    MIGRATION_MATRICES_ = obj.at("migration_matrices").get<decltype(MIGRATION_MATRICES_)>();
}

void Individual::to_json(json::json& obj) {HERE;
    obj["natural_mortality"] = NATURAL_MORTALITY_;
    obj["fishing_mortality"] = FISHING_MORTALITY_;
    obj["weight_for_age"] = WEIGHT_FOR_AGE_;
    obj["migration_matrices"] = MIGRATION_MATRICES_;
}

bool Individual::has_survived(const uint_fast32_t year, urbg_t& g) const {
    const auto age = year - birth_year_;
    return (wtl::generate_canonical(g) > NATURAL_MORTALITY_[age])
        && (wtl::generate_canonical(g) > FISHING_MORTALITY_[age]);
}

std::vector<std::string> Individual::names() {
    return {"id", "father_id", "mother_id", "birth_year", "location"};
}

std::ostream& Individual::write(std::ostream& ost) const {
    return ost << id_ << "\t"
               << father_id_ << "\t" << mother_id_ << "\t"
               << birth_year_ << "\t" << location_;
}

//! shortcut Individual::write(ost)
std::ostream& operator<<(std::ostream& ost, const Individual& x) {
    return x.write(ost);
}

void Individual::test() {HERE;
    Individual x;
    std::cout << x << std::endl;
}

} // namespace pbt
