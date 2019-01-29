#pragma once
#ifndef PBT_RANDOM_FWD_HPP
#define PBT_RANDOM_FWD_HPP

#include <random>

#ifdef SFMT_FOUND
namespace wtl {class sfmt19937_64;}
namespace pbf {
  using URBG = wtl::sfmt19937_64;
}
#else
namespace pbf {
  using URBG = std::mt19937_64;
}
#endif

#endif//PBT_RANDOM_FWD_HPP
