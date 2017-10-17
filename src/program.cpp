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
#include <sfmt.hpp>

#include <iostream>

namespace pbt {

namespace po = boost::program_options;

//! options description for general arguments
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
        ("help,h", po::bool_switch(), "print this help")
        ("verbose,v", po::bool_switch(), "verbose output")
        ("test", po::value<int>()->default_value(0)->implicit_value(1));
    return description;
}

/*! @ingroup params

    Command line option | Symbol  | Variable
    ------------------- | ------- | -------------------------
    `-n,--popsize`      | \f$N\f$ | Program::pop_size_
    `-g,--maxtime`      |         | Program::max_time_
    `-j,--parallel`     |         | Program::concurrency_
    `-o,--outdir`       |         | Program::out_dir_
*/
po::options_description Program::options_desc() {HERE;
    po::options_description description("Program");
    description.add_options()
      ("popsize,n", po::value(&pop_size_)->default_value(pop_size_))
      ("maxtime,t", po::value(&max_time_)->default_value(max_time_))
      ("parallel,j", po::value(&concurrency_)->default_value(concurrency_))
      ("outdir,o", po::value(&out_dir_)->default_value(out_dir_));
    // description.add(Population::options_desc());
    description.add(Individual::options_desc());
    return description;
}

void Program::help_and_exit() {HERE;
    auto description = general_desc();
    description.add(options_desc());
    // do not print positional arguments as options
    std::cout << "Usage: pbt [options]\n" << std::endl;
    description.print(std::cout);
    throw wtl::ExitSuccess();
}

//! Unit test for each class
inline void test(const int flg) {HERE;
    switch (flg) {
      case 0:
        break;
      case 1:
        Individual::test();
        Population::test();
        throw wtl::ExitSuccess();
      default:
        throw std::runtime_error("Unknown argument for --test");
    }
}

Program::Program(const std::vector<std::string>& arguments) {HERE;
    std::cout << wtl::join(arguments, " ") << std::endl;
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
    Individual::init_static_members();

    config_string_ = wtl::flags_into_string(vm);
    if (vm["verbose"].as<bool>()) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_string_ << std::endl;
    }
    test(vm["test"].as<int>());
}

void Program::run() {HERE;
    try {
        main();
    } catch (const wtl::KeyboardInterrupt& e) {
        std::cerr << e.what() << std::endl;
    }
}

void Program::main() {HERE;
    wtl::ChDir cd_outdir(out_dir_, true);
    // wtl::opfstream{"program_options.conf"} << config_string_;
    Population pop(pop_size_);
    pop.run(max_time_);
}

} // namespace pbt
