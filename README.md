# tekka

[![Build status](https://github.com/heavywatal/tekka/workflows/build/badge.svg)](https://github.com/heavywatal/tekka/actions)

Individual-based simulator of pacific bluefin tuna.

[Project page on GitHub](https://github.com/heavywatal/tekka)


## Requirements

- Unix-like environment (macOS, Linux, WSL, MinGW on MSYS2, etc.)
- C++17 compiler (clang++ >= Apple clang 11.0, g++ >= 9.1)
- [CMake](https://cmake.org/)

The following libraries are optional or automatically installed:

- [clippson](https://github.com/heavywatal/clippson)
- [cxxwtl](https://github.com/heavywatal/cxxwtl)
- [pcglite](https://github.com/heavywatal/pcglite)
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
DESTINATION=${HOME}/local
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$DESTINATION
cmake --build build -j 2
cmake --install build
PATH=${DESTINATION}/bin:$PATH
```

## Usage

```sh
tekka --help
```

All the parameters are described also in [Parameters](https://heavywatal.github.io/tekka/group__parameters.html) page.
