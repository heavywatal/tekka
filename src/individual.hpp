// -*- mode: c++; coding: utf-8 -*-
/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include <iosfwd>
#include <cstdint>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

/*! @brief Individual class
*/
class Individual {
  public:
    //! default constructor
    Individual() = default;

    //! write
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);
    //! unit test
    static void test();
  private:
    //! ID
    const uint_fast64_t id_ = 0UL;
};

} // namespace pbt

#endif /* PBT_INDIVIDUAL_HPP_ */
