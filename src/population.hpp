/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

#include <iosfwd>
#include <vector>

#include <sfmt.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

class Individual;

/*! @brief Population class
*/
class Population {
    using URBG = wtl::sfmt19937_64;
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
    void reproduce(URBG&);

    //! evaluate survival
    void survive(URBG&);

    //! evaluate migration
    void migrate(URBG&);

    //! Individual array males
    std::vector<Individual> males_;
    //! Individual array females
    std::vector<Individual> females_;
    //! year
    uint_fast32_t year_ = 0;
};

} // namespace pbt

#endif /* PBT_POPULATION_HPP_ */
