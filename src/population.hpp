/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

#include "common.hpp"

#include <iosfwd>
#include <vector>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

class Individual;

/*! @brief Population class
*/
class Population {
  public:
    //! constructor
    Population(const size_t initial_size);

    //! main iteration
    void run(const uint_fast32_t years);

    //! write
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);
    //! unit test
    static void test();
  private:
    //! give birth to children
    void reproduce();

    //! evaluate survival
    void survive();

    //! evaluate migration
    void migrate();

    //! sample
    void sample(const size_t n, std::ostream&);

    //! Individual array males
    std::vector<Individual> males_;
    //! Individual array females
    std::vector<Individual> females_;
    //! year
    uint_fast32_t year_ = 0;
    //! random bit generator
    urbg_t engine_;
};

} // namespace pbt

#endif /* PBT_POPULATION_HPP_ */
