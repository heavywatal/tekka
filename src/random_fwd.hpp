#pragma once
#ifndef PBF_RANDOM_FWD_HPP
#define PBF_RANDOM_FWD_HPP

#include <random>

namespace wtl {
  class sfmt19937_64;
}

namespace pbf {
  using URBG = wtl::sfmt19937_64;
}

#endif//PBF_RANDOM_FWD_HPP
