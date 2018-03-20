# tekkamaki

R interface to [blackthunnus](https://github.com/heavywatal/blackthunnus)

## Installation

1.  Install [blackthunnus](https://github.com/heavywatal/blackthunnus) with Homebrew/Linuxbrew or CMake.

1.  Install [devtools](https://github.com/hadley/devtools) in R:
    `install.packages('devtools')`

1.  Execute `devtools::install_github('heavywatal/blackthunnus/r')` in R.
    You may need `Sys.setenv(CMAKE_PREFIX_PATH='/prefix/to/blackthunnus')` to tell R the location of blackthunnus installation if it is not in the standard paths such as `/usr/local/`.
