/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "program.hpp"
#include "population.hpp"
#include "parameters.hpp"
#include "config.hpp"

#include <wtl/chrono.hpp>
#include <wtl/debug.hpp>
#include <wtl/zlib.hpp>
#include <clippson/clippson.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <cstdlib>

namespace pbf {

//! Options description for general purpose
inline clipp::group general_options(nlohmann::json* vm) {
    return (
      clippson::option(vm, {"h", "help"}, false, "Print this help"),
      clippson::option(vm, {"version"}, false, "Print version"),
      clippson::option(vm, {"default"}, false, "Print default parameters in json"),
      clippson::option(vm, {"v", "verbose"}, false, "Verbose output"),
      clippson::option(vm, {"i", "infile"}, std::string{}, "config file in json format")
    ).doc("General:");
}

//! Global options for simulation
/*! @ingroup params

    Command line option  | Variable
    -------------------- | ------------------
    `--seed`             | Parameters::seed
    `-i,--infile`        | -
    `-o,--outdir`        | Parameters::outdir
    `-y,--years`         | Parameters::years
    `-O,--origin`        | Parameters::origin
*/
inline clipp::group program_options(Parameters& params) {
    return (
      clippson::option({"seed"}, &params.seed),
      clippson::option({"o", "outdir"}, &params.outdir),
      clippson::option({"y", "years"}, &params.years, "Duration of simulation"),
      clippson::option({"O", "origin"}, &params.origin,
        "Initial number of juveniles at location 0 in season 0 in year 0; falls back to -R or -K.")
    ).doc("Initialization:");
}

//! Options for Population::reproduce
/*! @ingroup params

    Command line option      | Symbol   | Variable
    ------------------------ | -------- | -------------------------------
    `-R,--med_recruitment`   | \f$\operatorname{med}(R)\f$  | Parameters::med_recruitment
    `-S,--sigma_recruitment` | \f$\sigma_{R}\f$  | Parameters::sigma_recruitment
    `-K,--carrying_capacity` | \f$K\f$  | Parameters::carrying_capacity
    `-k,--overdispersion`    | \f$k\f$  | Parameters::overdispersion
    `-r,--recruitment`       | \f$r\f$  | Parameters::recruitment
*/
inline clipp::group reproduction_options(Parameters& params) {
    return (
      clippson::option({"R", "med_recruitment"}, &params.med_recruitment),
      clippson::option({"S", "sigma_recruitment"}, &params.sigma_recruitment),
      clippson::option({"K", "carrying_capacity"}, &params.carrying_capacity),
      clippson::option({"k", "overdispersion"}, &params.overdispersion,
        "k ∈ (0, ∞); equivalent to Poisson when k→∞ (or k<0 for convience)"
      ),
      clippson::option({"r", "recruitment"}, &params.recruitment)
    ).doc("Reproduction:");
}

//! Options for Population::sample
/*! @ingroup params

    Command line option           | Variable
    ----------------------------- | --------------------------------
    `-l,--last`                   | Parameters::last
    `--sa,--sample_size_adult`    | Parameters::sample_size_adult
    `--sj,--sample_size_juvenile` | Parameters::sample_size_juvenile
*/
inline clipp::group sampling_options(Parameters& params) {
    return (
      clippson::option({"l", "last"}, &params.last, "Sample last _ years"),
      clippson::option({"sa", "sample_size_adult"}, &params.sample_size_adult, "per location"),
      clippson::option({"sj", "sample_size_juvenile"}, &params.sample_size_juvenile, "per location")
    ).doc("Sampling:");
}

//! @cond

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(Parameters,
  natural_mortality, fishing_mortality, fishing_coef, weight_for_age, migration_matrices,
  origin, years, last, sample_size_adult, sample_size_juvenile, outdir, seed, carrying_capacity, recruitment, overdispersion
);

//! @endcond

Program::~Program() = default;

Program::Program(const std::vector<std::string>& arguments) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    nlohmann::json vm_local;
    try {
        clippson::parse(general_options(&vm_local), arguments);
    } catch (std::invalid_argument&) {}
    if (vm_local.at("version")) {
        std::cout << PROJECT_VERSION << "\n";
        throw exit_success();
    }
    Parameters params;
    if (vm_local.at("default")) {
        nlohmann::json obj{params};
        std::cout << obj << "\n";
        throw exit_success();
    }
    if (const std::string infile = vm_local.at("infile"); !infile.empty()) {
        std::ifstream ifs{infile};
        nlohmann::json obj;
        ifs >> obj;
        obj.get_to(params);
    }
    auto cli = (
      general_options(&vm_local),
      program_options(params),
      reproduction_options(params),
      sampling_options(params)
    );
    clippson::parse(cli, arguments);
    if (vm_local.at("help")) {
        auto fmt = clippson::doc_format();
        std::cout << "Usage: " << PROJECT_NAME << " [options]\n\n";
        std::cout << clipp::documentation(cli, fmt) << "\n";
        throw exit_success();
    }
    outdir_ = params.outdir;
    config_ = nlohmann::json{params}.dump() + "\n";
    if (vm_local.at("verbose")) {
        std::cerr << wtl::iso8601datetime() << std::endl;
        std::cerr << config_ << std::endl;
    }
    population_ = std::make_unique<Population>(params);
}

void Program::run() {
    population_->run();
}

void Program::write() const {
    namespace fs = std::filesystem;
    const auto outdir = fs::path(outdir_);
    if (!outdir.empty()) {
        fs::create_directory(outdir);
        std::ofstream{outdir / "config.json"} << config_;
        {
            wtl::zlib::ofstream ost{outdir / "sample_family.tsv.gz"};
            population_->write_sample_family(ost);
        }
        {
            wtl::zlib::ofstream ost{outdir / "demography.tsv.gz"};
            population_->write_demography(ost);
        }
    } else {
        population_->write_demography(std::cout);
    }
}

Parameters::Parameters():
  outdir{wtl::strftime("thunnus_%Y%m%d_%H%M%S")},
  seed{static_cast<int32_t>(std::random_device{}())} // 32-bit signed integer for R
  {
    std::istringstream iss(default_values);
    read(iss);
}

void Parameters::read(std::istream& ist) {
    nlohmann::json obj;
    ist >> obj;
    obj.at("natural_mortality").get_to(natural_mortality);
    obj.at("fishing_mortality").get_to(fishing_mortality);
    fishing_coef = obj.value("fishing_coef", fishing_coef);
    obj.at("weight_for_age").get_to(weight_for_age);
    obj.at("migration_matrices").get_to(migration_matrices);
}

std::ostream& Parameters::write(std::ostream& ost) const {
    nlohmann::json obj;
    obj["natural_mortality"] = natural_mortality;
    obj["fishing_mortality"] = fishing_mortality;
    obj["fishing_coef"] = fishing_coef;
    obj["weight_for_age"] = weight_for_age;
    obj["migration_matrices"] = migration_matrices;
    return ost << obj;
}

} // namespace pbf
