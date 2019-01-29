/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/program.hpp"
#include "src/population.hpp"

#include <wtl/filesystem.hpp>
#ifdef ZLIB_FOUND
  #include <wtl/zlib.hpp>
#endif

#include <iostream>
#include <fstream>
#include <stdexcept>

//! Output results to files
void write(const pbf::Program& program) {
  #ifdef ZLIB_FOUND
    using ofstream = wtl::zlib::ofstream;
    const std::string ext = ".tsv.gz";
  #else
    using ofstream = std::ofstream;
    const std::string ext = ".tsv";
  #endif
    const auto& population = program.population();
    const auto outdir = program.outdir();
    if (!outdir.empty()) {
        wtl::ChDir cd(outdir, true);
        std::ofstream{"config.json"} << program.config();
        {
            ofstream ost{"sample_family" + ext};
            population.write_sample_family(ost);
        }
        {
            ofstream ost{"demography" + ext};
            population.write_demography(ost);
        }
        ofstream ost{"msout" + ext};
        program.write_ms(ost);
    } else {
        program.write_ms(std::cout);
    }
}

//! Just instantiate and run Program
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        pbf::Program program(arguments);
        program.run();
        write(program);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
