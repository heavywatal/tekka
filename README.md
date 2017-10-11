# blackthunnus

Individual-based simulator of pacific bluefin tuna

[Project page on GitHub](https://github.com/heavywatal/blackthunnus)


## Dependencies

- Unix-like OS (macOS, Linux, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/)
- [Boost C++ Libraries](http://www.boost.org/) (>= 1.64.0)
- [sfmt-class](https://github.com/heavywatal/sfmt-class)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)


## Installation

Use `BOOST_ROOT` environment variable so that CMake can find your boost library.

```sh
git clone https://github.com/heavywatal/blackthunnus.git
mkdir build-blackthunnus
cd build-blackthunnus/
# export BOOST_ROOT=${HOME}/local
cmake -DCMAKE_INSTALL_PREFIX=${HOME}/local ../blackthunnus
cmake --build .
cmake --build . --target install
```
