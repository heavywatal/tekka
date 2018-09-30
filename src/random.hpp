#pragma once
#ifndef PBT_RANDOM_HPP
#define PBT_RANDOM_HPP

#include "random_fwd.hpp"
#ifdef SFMT_FOUND
  #include <sfmt.hpp>
#endif
#include <wtl/random.hpp>

namespace pbt {

inline URBG& engine64() {
#ifdef SFMT_FOUND
    return wtl::sfmt64();
#else
    return wtl::mt64();
#endif
}

}

#endif//PBT_RANDOM_HPP
