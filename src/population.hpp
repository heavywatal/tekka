// -*- mode: c++; coding: utf-8 -*-
/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

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

    //! give birth to children
    void reproduce();

    //! write
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);
    //! unit test
    static void test();
  private:
    //! Individual array males
    std::vector<Individual> males_;
    //! Individual array females
    std::vector<Individual> females_;
    //! by the quater-year
    uint_fast32_t time_ = 0;
};

} // namespace pbt

#endif /* PBT_POPULATION_HPP_ */
