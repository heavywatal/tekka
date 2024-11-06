/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "program.hpp"
#include "population.hpp"
#include "individual.hpp"
#include "config.hpp"

#include <wtl/chrono.hpp>
#include <wtl/debug.hpp>
#ifdef ZLIB_FOUND
  #include <wtl/zlib.hpp>
#endif
#include <clippson/clippson.hpp>

#include <filesystem>
#include <fstream>
#include <cstdlib>

namespace pbf {

//! Global variables mapper of command-line arguments
nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      clippson::option(vm, {"h", "help"}, false, "Print this help"),
      clippson::option(vm, {"version"}, false, "Print version"),
      clippson::option(vm, {"v", "verbose"}, false, "Verbose output"),
      clippson::option(vm, {"default"}, false, "Print default parameters in json")
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
      clippson::option(vm, {"O", "origin"}, 0.2, "Initial population size relative to K"),
      clippson::option(vm, {"y", "years"}, 100, "Duration of simulation"),
      clippson::option(vm, {"l", "last"}, 3, "Sample last _ years"),
      clippson::option(vm, {"sa", "sample_size_adult"}, std::vector<size_t>{10u, 10u}, "per location"),
      clippson::option(vm, {"sj", "sample_size_juvenile"}, std::vector<size_t>{10u, 10u}, "per location"),
      clippson::option(vm, {"i", "infile"}, std::string(""), "config file in json format"),
      clippson::option(vm, {"o", "outdir"}, OUT_DIR),
      clippson::option(vm, {"seed"}, seed)
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
      clippson::option(vm, {"r", "recruitment"}, 2.0),
      clippson::option(vm, {"K", "carrying_capacity"}, 1e3),
      clippson::option(vm, {"k", "overdispersion"}, -1.0,
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
    clippson::parse(cli, arguments);
    if (vm_local.at("help")) {
        auto fmt = clippson::doc_format();
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw exit_success();
    }
    if (vm_local.at("version")) {
        std::cout << PROJECT_VERSION << "\n";
        throw exit_success();
    }
    if (vm_local.at("default")) {
        Individual::JSON.write(std::cout);
        throw exit_success();
    }
    const std::string infile = VM.at("infile");
    if (!infile.empty()) {
        std::ifstream ifs(infile);
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

void Program::write() const {
    namespace fs = std::filesystem;
  #ifdef ZLIB_FOUND
    using ofstream = wtl::zlib::ofstream;
    const std::string ext = ".tsv.gz";
  #else
    using ofstream = std::ofstream;
    const std::string ext = ".tsv";
  #endif
    const auto outdir = fs::path(VM.at("outdir"));
    if (!outdir.empty()) {
        fs::create_directory(outdir);
        std::ofstream{outdir / "config.json"} << config_;
        {
            ofstream ost{outdir / ("sample_family" + ext)};
            population_->write_sample_family(ost);
        }
        {
            ofstream ost{outdir / ("demography" + ext)};
            population_->write_demography(ost);
        }
    } else {
        population_->write_demography(std::cout);
    }
}

} // namespace pbf
