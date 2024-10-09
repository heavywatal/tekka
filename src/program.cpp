/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "program.hpp"
#include "population.hpp"
#include "individual.hpp"
#include "config.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <clippson/clippson.hpp>

namespace pbf {

//! Global variables mapper of command-line arguments
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      wtl::option(vm, {"h", "help"}, false, "Print this help"),
      wtl::option(vm, {"version"}, false, "Print version"),
      wtl::option(vm, {"v", "verbose"}, false, "Verbose output"),
      wtl::option(vm, {"default"}, false, "Print default parameters in json")
    ).doc("General:");
}

//! Program options
/*! @ingroup params

    Command line option           | Symbol
    ----------------------------- | -------
    `-O,--origin`                 | -
    `-y,--years`                  | -
    `-l,--last`                   | -
    `--sa,--sample_size_adult`    | -
    `--sj,--sample_size_juvenile` | -
    `-i,--infile`                 | -
    `-o,--outdir`                 | -
*/
inline clipp::group program_options(nlohmann::json* vm) {
    const std::string OUT_DIR = wtl::strftime("thunnus_%Y%m%d_%H%M%S");
    const int seed = static_cast<int>(std::random_device{}()); // 32-bit signed integer for R
    return (
      wtl::option(vm, {"O", "origin"}, 0.2, "Initial population size relative to K"),
      wtl::option(vm, {"y", "years"}, 100, "Duration of simulation"),
      wtl::option(vm, {"l", "last"}, 3, "Sample last _ years"),
      wtl::option(vm, {"sa", "sample_size_adult"}, std::vector<size_t>{10u, 10u}, "per location"),
      wtl::option(vm, {"sj", "sample_size_juvenile"}, std::vector<size_t>{10u, 10u}, "per location"),
      wtl::option(vm, {"i", "infile"}, std::string(""), "config file in json format"),
      wtl::option(vm, {"o", "outdir"}, OUT_DIR),
      wtl::option(vm, {"seed"}, seed)
    ).doc("Program:");
}

//! Reproduction options
/*! @ingroup params

    Command line option      | Symbol   | Variable
    ------------------------ | -------- | -------------------------------
    `-r,--recruitment`       | \f$r\f$  | Population::recruitment_coef_
    `-K,--carrying_capacity` | \f$K\f$  | Population::carrying_capacity_
    `-k,--overdispersion`    | \f$k\f$  | Population::k_nbinom_
*/
inline clipp::group reproduction_options(nlohmann::json* vm) {
    return (
      wtl::option(vm, {"r", "recruitment"}, 2.0),
      wtl::option(vm, {"K", "carrying_capacity"}, 1e3),
      wtl::option(vm, {"k", "overdispersion"}, -1.0,
        "k ∈ (0, ∞); equivalent to Poisson when k→∞ (or k<0 for convience)"
      )
    ).doc("Reproduction:");
}

Program::Program(const std::vector<std::string>& arguments)
: command_args_(arguments) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    nlohmann::json vm_local;
    auto cli = (
      general_options(&vm_local),
      program_options(&VM),
      reproduction_options(&VM)
    );
    wtl::parse(cli, arguments);
    if (vm_local.at("help")) {
        auto fmt = wtl::doc_format();
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local.at("version")) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local.at("default")) {
        Individual::JSON.write(std::cout);
        throw wtl::ExitSuccess();
    }
    const std::string infile = VM.at("infile");
    if (!infile.empty()) {
        auto ifs = wtl::make_ifs(infile);
        Individual::JSON.read(ifs);
    }
    Individual::set_dependent_static(VM.at("years"));
    config_ = VM.dump(2) + "\n";
    if (vm_local.at("verbose")) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_ << std::endl;
        Individual::JSON.write(std::cerr) << std::endl;
    }
}

void Program::run() {
    const double K = VM.at("carrying_capacity");
    const double O = VM.at("origin");
    population_ = std::make_unique<Population>(
        static_cast<size_t>(K * O),
        VM.at("seed"),
        K,
        VM.at("recruitment"),
        VM.at("overdispersion")
    );
    population_->run(
        VM.at("years"),
        VM.at("sample_size_adult"),
        VM.at("sample_size_juvenile"),
        VM.at("last")
    );
}

std::string Program::outdir() const {
    return VM.at("outdir");
}

} // namespace pbf
