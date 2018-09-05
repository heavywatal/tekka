/*! @file program.cpp
    @brief Implementation of Program class
    @defgroup params Parameters
*/
#include "program.hpp"
#include "population.hpp"
#include "individual.hpp"
#include "segment.hpp"

#include <wtl/exception.hpp>
#include <wtl/debug.hpp>
#include <wtl/iostr.hpp>
#include <wtl/chrono.hpp>
#include <wtl/zlib.hpp>
#include <wtl/filesystem.hpp>
#include <wtl/getopt.hpp>
#include <sfmt.hpp>

#include <boost/program_options.hpp>

namespace pbt {

namespace fs = boost::filesystem;
namespace po = boost::program_options;

//! options description for general arguments
inline po::options_description general_desc() {HERE;
    po::options_description description("General");
    description.add_options()
      ("help,h", po::bool_switch(), "print this help")
      ("verbose,v", po::bool_switch(), "verbose output")
      ("quiet,q", po::bool_switch(), "suppress output")
    ;
    return description;
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
po::options_description Program::options_desc() {HERE;
    const std::string OUT_DIR = wtl::strftime("thunnus_%Y%m%d_%H%M%S");
    po::options_description description("Program");
    description.add_options()
      ("popsize,n", po::value<size_t>()->default_value(1000u), "Initial population size")
      ("years,y", po::value<uint_fast32_t>()->default_value(40u), "Duration of simulation")
      ("last,l", po::value<uint_fast32_t>()->default_value(2u), "Sample last _ years")
      ("sample,s", po::value<double>()->default_value(0.02), "per location")
      ("mutation,u", po::value<double>()->default_value(0.1), "per generation per haploid")
      ("default,d", po::bool_switch(), "Print default parameters in json")
      ("infile,i", po::value<std::string>(), "config file in json format")
      ("outdir,o", po::value<std::string>()->default_value("")->implicit_value(OUT_DIR))
      ("tree,t", po::bool_switch(), "Output family tree")
      ("seed", po::value<std::random_device::result_type>()->default_value(std::random_device{}()))
    ;
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

Program::Program(const std::vector<std::string>& arguments)
: vars_(std::make_unique<po::variables_map>()),
  command_args_(arguments) {HERE;
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);
    std::cout.precision(15);
    std::cerr.precision(6);

    auto description = general_desc();
    description.add(options_desc());
    // po::variables_map vm;
    auto& vm = *vars_;
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
}

Program::~Program() = default;

void Program::run() {HERE;
    auto& vm = *vars_;
    const auto quiet = vm["quiet"].as<bool>();
    const auto writing_tree = vm["tree"].as<bool>();
    const auto popsize = vm["popsize"].as<size_t>();
    const auto years = vm["years"].as<uint_fast32_t>();
    const auto last_years = vm["last"].as<uint_fast32_t>();
    const auto sampling_rate = vm["sample"].as<double>();
    const auto mutation_rate = vm["mutation"].as<double>();
    const auto outdir = vm["outdir"].as<std::string>();
    const auto seed = vm["seed"].as<std::random_device::result_type>();
    population_ = std::make_unique<Population>(popsize, seed);
    population_->run(years, sampling_rate, last_years);
    if (quiet) return;
    if (!outdir.empty()) {
        wtl::ChDir cd(outdir, true);
        wtl::make_ofs("program_options.conf") << config_string_;
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

} // namespace pbt
