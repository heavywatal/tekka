// -*- mode: c++; coding: utf-8 -*-
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

uint_fast32_t Individual::MEAN_CLUTCH_SIZE_ = 4;
uint_fast32_t Individual::LAST_ID_ = 0;
std::poisson_distribution<uint_fast32_t> Individual::POISSON_CLUTCH_SIZE_;

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option | Symbol         | Variable
    ------------------- | -------------- | -------------------------------
    `-c,--clutch`       | -              | Individual::MEAN_CLUTCH_SIZE_
*/
boost::program_options::options_description Individual::options_desc() {
    namespace po = boost::program_options;
    po::options_description desc{"Tissue"};
    desc.add_options()
        ("clutch,c", po::value(&MEAN_CLUTCH_SIZE_)->default_value(MEAN_CLUTCH_SIZE_))
    ;
    return desc;
}

void Individual::init_static_members() {HERE;
    POISSON_CLUTCH_SIZE_.param(decltype(POISSON_CLUTCH_SIZE_)::param_type(MEAN_CLUTCH_SIZE_));
}

bool Individual::has_survived(const uint_fast32_t time) const {
    const uint_fast32_t age = (time - birth_date_) / 4U;
    return (wtl::sfmt()() > age);
}

std::ostream& Individual::write(std::ostream& ost) const {
    return ost << id_ << ":"
               << father_id_ << ":" << mother_id_ << ":"
               << birth_date_;
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
