/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "program.hpp"
#include "population.hpp"
#include "individual.hpp"
#include "segment.hpp"
#include "random.hpp"
#include "config.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>
#include <clippson/clippson.hpp>

namespace pbt {

nlohmann::json VM;

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {HERE;
    return (
      wtl::option(vm, {"h", "help"}, false, "print this help"),
      wtl::option(vm, {"version"}, false, "print version"),
      wtl::option(vm, {"v", "verbose"}, false, "verbose output"),
      wtl::option(vm, {"default"}, false, "Print default parameters in json")
    ).doc("General:");
}

/*! @ingroup params

    Command line option | Symbol
    ------------------- | -------
    `-n,--popsize`      | \f$N\f$
    `-y,--years`        | -
    `-l,--last`         | -
    `-s,--sample`       | -
    `-u,--mutation`     | -
*/
inline clipp::group program_options(nlohmann::json* vm) {HERE;
    const std::string OUT_DIR = wtl::strftime("thunnus_%Y%m%d_%H%M%S");
    const int seed = std::random_device{}(); // 32-bit signed integer for R
    return (
      wtl::option(vm, {"n", "popsize"}, 1000u, "Initial population size"),
      wtl::option(vm, {"y", "years"}, 40u, "Duration of simulation"),
      wtl::option(vm, {"l", "last"}, 2u, "Sample last _ years"),
      wtl::option(vm, {"s", "sample"}, 0.02, "per location"),
      wtl::option(vm, {"u", "mutation"}, 0.1, "per generation per haploid"),
      wtl::option(vm, {"i", "infile"}, std::string(""), "config file in json format"),
      wtl::option(vm, {"o", "outdir"}, OUT_DIR),
      wtl::option(vm, {"q", "quiet"}, false, "suppress output"),
      wtl::option(vm, {"t", "tree"}, false, "Output family tree"),
      wtl::option(vm, {"seed"}, seed)
    ).doc("Program:");
}

//! Program options
/*! @ingroup params
    @return Program options description

    Command line option  | Symbol         | Variable
    -------------------- | -------------- | -------------------------------
    `-r,--recruitment`   | -              | Individual::RECRUITMENT_COEF_
    `-k,--overdispersion`| -              | Individual::NEGATIVE_BINOM_K_
*/
inline clipp::group individual_options(nlohmann::json* vm, IndividualParams* p) {HERE;
    return (
      wtl::option(vm, {"r", "recruitment"}, &p->RECRUITMENT_COEF),
      wtl::option(vm, {"k", "overdispersion"}, &p->NEGATIVE_BINOM_K)
    ).doc("Individual");
}

Program::Program(const std::vector<std::string>& arguments)
: command_args_(arguments) {HERE;
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
    auto fmt = wtl::doc_format();
    if (vm_local.at("help")) {
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw wtl::ExitSuccess();
    }
    if (vm_local.at("version")) {
        std::cout << PROJECT_VERSION << "\n";
        throw wtl::ExitSuccess();
    }
    Individual::param(individual_params);
    Individual::set_default_values();
    if (vm_local.at("default")) {
        Individual::write_json(std::cout);
        throw wtl::ExitSuccess();
    }
    const std::string infile = VM.at("infile").get<std::string>();
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
    HERE;
}

Program::~Program() = default;

void Program::run() {HERE;
    const bool quiet = VM.at("quiet");
    const bool writing_tree = VM.at("tree");
    const size_t popsize = VM.at("popsize");
    const uint_fast32_t years = VM.at("years");
    const uint_fast32_t last_years = VM.at("last");
    const double sampling_rate = VM.at("sample");
    const double mutation_rate = VM.at("mutation");
    const std::string outdir = VM.at("outdir");
    const int seed = VM.at("seed");
    population_ = std::make_unique<Population>(popsize, static_cast<uint32_t>(seed));
    population_->run(years, sampling_rate, last_years);
    if (quiet) return;
    if (!outdir.empty()) {
        wtl::ChDir cd(outdir, true);
        wtl::make_ofs("config.json") << config_;
        if (writing_tree) {
            wtl::zlib::ofstream ost{"sample_family.tsv.gz"};
            population_->write_sample_family(ost);
        }
        {
            wtl::zlib::ofstream ost{"demography.tsv.gz"};
            population_->write_demography(ost);
        }
        wtl::zlib::ofstream ost{"msout.txt.gz"};
        auto tree = population_->coalesce();
        population_->write_ms(tree, mutation_rate, ost);
    } else {
        if (writing_tree) {
            population_->write_sample_family(std::cout);
        } else {
            wtl::join(command_args_, std::cout, " ") << "\n";
            std::cout << seed << "\n";
            auto tree = population_->coalesce();
            population_->write_ms(tree, mutation_rate, std::cout);
        }
    }
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

std::streambuf* std_cout_rdbuf(std::streambuf* buf) {
    return std::cout.rdbuf(buf);
}
std::streambuf* std_cerr_rdbuf(std::streambuf* buf) {
    return std::cerr.rdbuf(buf);
}

} // namespace pbt
