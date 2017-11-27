/*! @file common.hpp
    @brief Define globals
*/
#pragma once
#ifndef PBT_COMMON_HPP_
#define PBT_COMMON_HPP_

#include <sfmt.hpp>
#include <json.hpp>

namespace pbt {

using urbg_t = wtl::sfmt19937_64;

namespace json = nlohmann;

} // namespace pbt

#endif /* PBT_COMMON_HPP_ */
