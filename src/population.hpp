/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

#include <cstdint>
#include <iosfwd>
#include <vector>
#include <map>
#include <memory>

namespace wtl {class sfmt19937_64;}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

using URBG = wtl::sfmt19937_64;
class Individual;

/*! @brief Population class
*/
class Population {
  public:
    //! constructor
    Population(const size_t initial_size);
    //! destructor
    ~Population();

    //! main iteration
    void run(const uint_fast32_t simulating_duration,
             const double sample_rate,
             const uint_fast32_t recording_duration=1u);

    //! count individuals for each location
    std::vector<size_t> sizes() const;
    //! write sampled individuals
    std::ostream& write_sample(std::ostream&) const;
    //! Construct and write tree from samples
    std::ostream& write_sample_family(std::ostream& ost) const;
    //! write
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

  private:
    //! give birth to children
    void reproduce();

    //! evaluate survival
    void survive(const uint_fast32_t quarter);

    //! evaluate migration
    void migrate();

    //! sample individuals
    void sample(const double rate);

    //! write column names for write_sample()
    std::ostream& write_sample_header(std::ostream&) const;

    //! Individual array males
    std::vector<std::shared_ptr<Individual>> males_;
    //! Individual array females
    std::vector<std::shared_ptr<Individual>> females_;
    //! samples: capture_year => individuals
    std::map<uint_fast32_t, std::vector<std::shared_ptr<Individual>>> year_samples_;
    //! year
    uint_fast32_t year_ = 0;
    //! random bit generator
    std::unique_ptr<URBG> engine_;
};

} // namespace pbt

#endif /* PBT_POPULATION_HPP_ */
