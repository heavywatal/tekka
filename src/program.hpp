#pragma once
#ifndef PBT_PROGRAM_HPP_
#define PBT_PROGRAM_HPP_

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

namespace pbf {

class Population;

//! @cond
class exit_success: public std::logic_error {
  public:
    exit_success() noexcept: std::logic_error("") {}
};
//! @endcond

/*! @brief Handle command line arguments and a Population instance
*/
class Program {
  public:
    //! Initialize with command-line arguments.
    Program(const std::vector<std::string>& args);
    ~Program();
    //! Top level function that should be called once from global main.
    void run();
    //! Output results to files
    void write() const;

  private:
    //! Written to "config.json"
    std::string config_{};
    //! Output directory
    std::string outdir_{};
    //! Population instance
    std::unique_ptr<Population> population_{nullptr};
};

} // namespace pbf

#endif /* PBT_PROGRAM_HPP_ */
