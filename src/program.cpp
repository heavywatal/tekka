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

//! Global variables mapper of commane-line arguments
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
    `-y,--years`                  | -
    `-l,--last`                   | -
    `--sa,--sample_size_adult`    | -
    `--sj,--sample_size_juvenile` | -
    `-u,--mutation`               | -
    `-i,--infile`                 | -
    `-o,--outdir`                 | -
*/
inline clipp::group program_options(nlohmann::json* vm) {
    const std::string OUT_DIR = wtl::strftime("thunnus_%Y%m%d_%H%M%S");
    const int seed = static_cast<int>(std::random_device{}()); // 32-bit signed integer for R
    return (
      wtl::option(vm, {"y", "years"}, 40u, "Duration of simulation"),
      wtl::option(vm, {"l", "last"}, 2u, "Sample last _ years"),
      wtl::option(vm, {"sa", "sample_size_adult"}, std::vector<size_t>{10u, 10u}, "per location"),
      wtl::option(vm, {"sj", "sample_size_juvenile"}, std::vector<size_t>{10u, 10u}, "per location"),
      wtl::option(vm, {"u", "mutation"}, 0.1, "per generation per haploid"),
      wtl::option(vm, {"i", "infile"}, std::string(""), "config file in json format"),
      wtl::option(vm, {"o", "outdir"}, OUT_DIR),
      wtl::option(vm, {"seed"}, seed)
    ).doc("Program:");
}

//! Individual options
/*! @ingroup params

    Command line option      | Symbol   | Variable
    ------------------------ | -------- | -------------------------------
    `-r,--recruitment`       | \f$r\f$  | IndividualParams::RECRUITMENT_COEF
    `-K,--carrying_capacity` | \f$K\f$  | IndividualParams::CARRYING_CAPACITY
    `-k,--overdispersion`    | \f$k\f$  | IndividualParams::NEGATIVE_BINOM_K
*/
inline clipp::group individual_options(nlohmann::json* vm, IndividualParams* p) {
    return (
      wtl::option(vm, {"r", "recruitment"}, &p->RECRUITMENT_COEF),
      wtl::option(vm, {"K", "carrying_capacity"}, &p->CARRYING_CAPACITY),
      wtl::option(vm, {"k", "overdispersion"}, &p->NEGATIVE_BINOM_K,
        "k ∈ (0, ∞); equivalent to Poisson when k→∞ (or k<0 for convience)"
      )
    ).doc("Individual:");
}

Program::Program(const std::vector<std::string>& arguments)
: command_args_(arguments) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    nlohmann::json vm_local;
    IndividualParams individual_params;
    auto cli = (
      general_options(&vm_local),
      program_options(&VM),
      individual_options(&VM, &individual_params)
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
    Individual::param(individual_params);
    if (vm_local.at("default")) {
        Individual::write_json(std::cout);
        throw wtl::ExitSuccess();
    }
    const std::string infile = VM.at("infile");
    if (!infile.empty()) {
        auto ifs = wtl::make_ifs(infile);
        Individual::read_json(ifs);
    }
    config_ = VM.dump(2) + "\n";
    if (vm_local.at("verbose")) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_ << std::endl;
        Individual::write_json(std::cerr);
    }
}

Program::~Program() = default;

void Program::run() {
    population_ = std::make_unique<Population>(
        VM.at("carrying_capacity"),
        VM.at("seed")
    );
    population_->run(
        VM.at("years"),
        VM.at("sample_size_adult"),
        VM.at("sample_size_juvenile"),
        VM.at("last")
    );
}

std::string Program::sample_family() const {
    std::ostringstream oss;
    population_->write_sample_family(oss);
    return oss.str();
}

std::string Program::demography() const {
    std::ostringstream oss;
    population_->write_demography(oss);
    return oss.str();
}

std::string Program::outdir() const {
    return VM.at("outdir");
}

//! std::cout.rdbuf
std::streambuf* std_cout_rdbuf(std::streambuf* buf) {
    return std::cout.rdbuf(buf);
}
//! std::cerr.rdbuf
std::streambuf* std_cerr_rdbuf(std::streambuf* buf) {
    return std::cerr.rdbuf(buf);
}

} // namespace pbf
