/*! @file program.hpp
    @brief Interface of Program class
*/
#pragma once
#ifndef PBT_PROGRAM_HPP_
#define PBT_PROGRAM_HPP_

#include <vector>
#include <string>

#include <boost/program_options.hpp>

namespace pbt {

/*! @brief Program class
*/
class Program {
  public:
    //! parse command arguments
    Program(const std::vector<std::string>& args);
    //! parse command arguments
    Program(int argc, char* argv[])
    : Program(std::vector<std::string>(argv, argv + argc)) {}

    //! top level function that should be called once from global main
    void run();

  private:
    //! called from run()
    void main();
    //! options description for Program class
    boost::program_options::options_description options_desc();
    //! print help message and exit
    [[noreturn]] void help_and_exit();

    //! initial population size
    size_t pop_size_ = 1000u;
    //! number of samples per year
    size_t sample_size_ = 10u;
    //! maximum years to simulate
    uint_fast32_t simulating_duration_ = 40u;
    //! last years to record samples
    uint_fast32_t recording_duration_ = 1u;
    //! number of threads
    unsigned int concurrency_ = 1u;
    //! `-w`
    bool is_writing_ = false;
    //! name of output directory
    std::string out_dir_ = "";
    //! writen to "program_options.conf"
    std::string config_string_;
};

} // namespace pbt

#endif /* PBT_PROGRAM_HPP_ */
