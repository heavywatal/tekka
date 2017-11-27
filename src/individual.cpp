/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"
#include "individual.hpp"

#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/prandom.hpp>
#include <sfmt.hpp>

#include <iostream>

namespace pbt {

uint_fast32_t Individual::MEAN_MATING_NUMBER_ = 2;
uint_fast32_t Individual::MEAN_CLUTCH_SIZE_ = 4;
std::vector<double> Individual::NATURAL_MORTALITY_(MAX_AGE_, 0.1);
std::vector<double> Individual::FISHING_MORTALITY_(MAX_AGE_, 0.1);
uint_fast32_t Individual::LAST_ID_ = 0;

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option | Symbol         | Variable
    ------------------- | -------------- | -------------------------------
    `-m,--mating`       | -              | Individual::MEAN_MATING_NUMBER_
    `-c,--clutch`       | -              | Individual::MEAN_CLUTCH_SIZE_
*/
boost::program_options::options_description Individual::options_desc() {
    namespace po = boost::program_options;
    po::options_description desc{"Individual"};
    desc.add_options()
        ("mating,m", po::value(&MEAN_MATING_NUMBER_)->default_value(MEAN_MATING_NUMBER_))
        ("clutch,c", po::value(&MEAN_CLUTCH_SIZE_)->default_value(MEAN_CLUTCH_SIZE_))
    ;
    return desc;
}

bool Individual::has_survived(const uint_fast32_t year) const {
    const auto age = year - birth_year_;
    return (wtl::sfmt().canonical() > NATURAL_MORTALITY_[age])
        && (wtl::sfmt().canonical() > FISHING_MORTALITY_[age]);
}

std::ostream& Individual::write(std::ostream& ost) const {
    return ost << id_ << ":"
               << father_id_ << ":" << mother_id_ << ":"
               << birth_year_;
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
