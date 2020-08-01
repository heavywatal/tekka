# tekka

[![Build status](https://github.com/heavywatal/tekka/workflows/build/badge.svg)](https://github.com/heavywatal/tekka/actions)

Individual-based simulator of pacific bluefin tuna.

[Project page on GitHub](https://github.com/heavywatal/tekka)


## Requirements

- Unix-like environment (macOS, Linux, WSL, MinGW on MSYS2, etc.)
- C++14 compiler (clang++ >= Apple LLVM 8.1, g++ >= 5.3)
- [CMake](https://cmake.org/) (>= 3.8.0)

The following libraries are optional or automatically installed:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [zlib](https://zlib.net)


## R interface

You can install and use this program via [R package "tekkamaki"](https://heavywatal.github.io/tekkamaki/).


## Installation of command-line version

The easiest way is to use [Homebrew](https://brew.sh/).
The following command installs tekka and all the dependencies:
```sh
brew install heavywatal/tap/tekka
```

Alternatively, you can get the source code from GitHub manually:
```sh
git clone https://github.com/heavywatal/tekka.git
cd tekka/
mkdir build
cd build/
YOUR_PREFIX=${HOME}/local  # or /usr/local
cmake -DCMAKE_INSTALL_PREFIX=$YOUR_PREFIX ..
make -j2
make install
```
