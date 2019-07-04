/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

#include "random_fwd.hpp"

#include <cstdint>
#include <iosfwd>
#include <vector>
#include <list>
#include <map>
#include <memory>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbf {

class Individual;

/*! @brief Population class
*/
class Population {
  public:
    //! constructor
    Population(const size_t initial_size, uint_fast32_t seed);
    //! destructor
    ~Population();

    //! main iteration
    void run(const int_fast32_t simulating_duration,
             const std::vector<size_t>& sample_size_adult={1u, 1u},
             const std::vector<size_t>& sample_size_juvenile={1u,1u},
             const int_fast32_t recording_duration=1);

    //! Construct and write tree from samples
    std::ostream& write_sample_family(std::ostream& ost) const;
    //! write #demography_
    std::ostream& write_demography(std::ostream&) const;
    //! write
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

  private:
    //! give birth to children
    void reproduce();

    //! give birth to children
    void reproduce(uint_fast32_t location, double density_effect);

    //! evaluate survival
    void survive();

    //! evaluate migration
    void migrate();

    //! sample individuals
    void sample(std::vector<std::vector<std::shared_ptr<Individual>>>* subpops,
                const std::vector<size_t>& sample_sizes);

    //! append current state to #demography_
    void append_demography(int_fast32_t season);

    //! Count individuals for each location and age
    std::vector<std::vector<uint_fast32_t>> count(int_fast32_t season) const;

    size_t num_subpops() const noexcept {return subpopulations_.size();}

    //! Individual array for each subpopulation
    std::vector<std::vector<std::shared_ptr<Individual>>> subpopulations_;
    //! first-year individuals
    std::vector<std::vector<std::shared_ptr<Individual>>> juveniles_subpops_;
    //! Counts of juveniles; [[number for each location] for each season]
    std::vector<std::vector<uint_fast32_t>> juveniles_demography_;
    //! samples: capture_year => individuals
    std::vector<std::map<int_fast32_t, std::vector<std::shared_ptr<Individual>>>> loc_year_samples_;
    //! (year, season) => [[count for each age] for each location]
    std::map<std::pair<int_fast32_t, int_fast32_t>, std::vector<std::vector<uint_fast32_t>>> demography_;
    //! year
    int_fast32_t year_ = 0;
    //! random bit generator
    std::unique_ptr<URBG> engine_;
};

} // namespace pbf

#endif /* PBT_POPULATION_HPP_ */
