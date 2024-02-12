/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "program.hpp"
#include "population.hpp"

#ifdef ZLIB_FOUND
  #include <wtl/zlib.hpp>
#endif

#include <filesystem>
#include <iostream>
#include <fstream>
#include <stdexcept>

//! Output results to files
void write(const pbf::Program& program) {
    namespace fs = std::filesystem;
  #ifdef ZLIB_FOUND
    using ofstream = wtl::zlib::ofstream;
    const std::string ext = ".tsv.gz";
  #else
    using ofstream = std::ofstream;
    const std::string ext = ".tsv";
  #endif
    const auto& population = program.population();
    const auto outdir = fs::path(program.outdir());
    if (!outdir.empty()) {
        fs::create_directory(outdir);
        std::ofstream{outdir / "config.json"} << program.config();
        {
            ofstream ost{outdir / ("sample_family" + ext)};
            population.write_sample_family(ost);
        }
        {
            ofstream ost{outdir / ("demography" + ext)};
            population.write_demography(ost);
        }
    } else {
        population.write_demography(std::cout);
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
