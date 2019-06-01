/*! @file program.hpp
    @brief Interface of Program class
*/
#pragma once
#ifndef PBT_PROGRAM_HPP_
#define PBT_PROGRAM_HPP_

#include <vector>
#include <string>
#include <memory>

namespace pbf {

class Population;

/*! @brief Program class
*/
class Program {
  public:
    //! parse command arguments
    Program(const std::vector<std::string>& args);
    //! destructor
    ~Program();
    //! top level function that should be called once from global main
    void run();

    //! @name Getter for main()
    //@{
    //! Get #population_
    const Population& population() const noexcept {return *population_;}
    //! Get #config_
    const std::string& config() const noexcept {return config_;}
    //! Get VM["outdir"]
    std::string outdir() const;
    //@}

    //! @name Output for Rcpp
    //@{
    //! Using Population.write_sample_family
    std::string sample_family() const;
    //! Using Population.write_demography
    std::string demography() const;
    //@}

  private:
    //! command line arguments
    std::vector<std::string> command_args_;
    //! writen to "config.json"
    std::string config_;
    //! Population instance
    std::unique_ptr<Population> population_;
};

//! @name Workaround for R/Rcpp
//@{
std::streambuf* std_cout_rdbuf(std::streambuf*);
std::streambuf* std_cerr_rdbuf(std::streambuf*);
//@}

} // namespace pbf

#endif /* PBT_PROGRAM_HPP_ */
