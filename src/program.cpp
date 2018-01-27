/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "program.hpp"
#include "population.hpp"
#include "individual.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/zfstream.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/getopt.hpp>
#include <wtl/chrono.hpp>
#include <sfmt.hpp>

#include <boost/program_options.hpp>

namespace pbt {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

const std::string OUT_DIR = wtl::strftime("thunnus_%Y%m%d_%H%M_") + std::to_string(::getpid());

//! options description for general arguments
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
    ;
    return description;
}

/*! @ingroup params

    Command line option | Symbol  | Variable
    ------------------- | ------- | -------------------------
    `-n,--popsize`      | \f$N\f$ | Program::pop_size_
    `-y,--years`        |         | Program::simulating_duration_
    `-l,--last`         |         | Program::recording_duration_
    `-s,--sample`       |         | Program::sample_size_
    `-j,--parallel`     |         | Program::concurrency_
    `-o,--outdir`       |         | Program::out_dir_
*/
po::options_description Program::options_desc() {HERE;
    po::options_description description("Program");
    description.add_options()
      ("popsize,n", po::value(&pop_size_)->default_value(pop_size_), "Initial population size")
      ("years,y", po::value(&simulating_duration_)->default_value(simulating_duration_))
      ("last,l", po::value(&recording_duration_)->default_value(recording_duration_), "Sample last _ years")
      ("sample,s", po::value(&sample_rate_)->default_value(sample_rate_))
      ("parallel,j", po::value(&concurrency_)->default_value(concurrency_))
      ("outdir-predefined,O", po::bool_switch(), OUT_DIR.c_str())
      ("default,d", po::bool_switch(), "Print default parameters in json")
      ("infile,i", po::value<std::string>(), "config file in json format")
      ("outdir,o", po::value(&out_dir_));
    // description.add(Population::options_desc());
    description.add(Individual::options_desc());
    return description;
}

[[noreturn]] void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: blackthunnus [options]\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    wtl::join(arguments, std::cerr, " ") << std::endl;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    auto description = general_desc();
    description.add(options_desc());
    po::variables_map vm;
    po::store(po::command_line_parser({arguments.begin() + 1, arguments.end()}).
              options(description).run(), vm);
    if (vm["help"].as<bool>()) {help_and_exit();}
    po::notify(vm);
    Individual::set_default_values();
    if (vm["default"].as<bool>()) {
        Individual::write_json(std::cout);
        throw wtl::ExitSuccess();
    }
    if (vm.count("infile")) {
        auto ifs = wtl::make_ifs(vm["infile"].as<std::string>());
        Individual::read_json(ifs);
    }
    config_string_ = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_string_ << std::endl;
        Individual::write_json(std::cerr);
    }
    if (vm["outdir-predefined"].as<bool>()) {
        if (vm.count("outdir")) {
            throw std::runtime_error("cannot use -o and -O at the same time");
        }
        out_dir_ = OUT_DIR;
    }
}

void Program::run() {HERE;
    try {
        main();
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    }
}

void Program::main() {HERE;
    Population pop(pop_size_);
    pop.run(simulating_duration_, sample_rate_, recording_duration_);
    if (!out_dir_.empty()) {
        wtl::ChDir cd(out_dir_, true);
        wtl::make_ofs("program_options.conf") << config_string_;
        wtl::ozfstream ost{"sample_family.tsv.gz"};
        pop.write_sample_family(ost);
    } else {
        pop.write_sample_family(std::cout);
    }
}

} // namespace pbt
