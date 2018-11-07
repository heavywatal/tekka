/*! @file main.cpp
    @brief Only defines tiny main()
*/
#include "src/program.hpp"
#include "src/population.hpp"

#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>

#include <iostream>
#include <stdexcept>

void write(const pbt::Program& program) {
    const auto& population = program.population();
    const auto outdir = program.outdir();
    if (!outdir.empty()) {
        wtl::ChDir cd(outdir, true);
        std::ofstream{"config.json"} << program.config();
        {
            wtl::zlib::ofstream ost{"sample_family.tsv.gz"};
            population.write_sample_family(ost);
        }
        {
            wtl::zlib::ofstream ost{"demography.tsv.gz"};
            population.write_demography(ost);
        }
        wtl::zlib::ofstream ost{"msout.txt.gz"};
        program.write_ms(ost);
    } else {
        program.write_ms(std::cout);
    }
}

//! Just instantiate and run Program
int main(int argc, char* argv[]) {
    std::vector<std::string> arguments(argv + 1, argv + argc);
    try {
        pbt::Program program(arguments);
        program.run();
        write(program);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
