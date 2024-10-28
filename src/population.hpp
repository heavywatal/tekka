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
#include <map>
#include <memory>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbf {

class Individual;

/*! @brief Population class
*/
class Population {
  public:
    //! Constructor
    Population(const size_t initial_size, std::random_device::result_type seed,
               const double carrying_capacity = 1e3,
               const double recruitment_coef = 2.0,
               const double negative_binom_k = -1.0);
    ~Population();

    //! Main iteration
    void run(const int_fast32_t simulating_duration,
             const std::vector<size_t>& sample_size_adult={1u, 1u},
             const std::vector<size_t>& sample_size_juvenile={1u,1u},
             const int_fast32_t recording_duration=1);

    //! Construct and write tree from #loc_year_samples_.
    //! Serial IDs are assigned to sampled individuals.
    //! IDs larger than sample size are assigned to unsampled ancestors.
    std::ostream& write_sample_family(std::ostream& ost) const;
    //! Write #demography_.
    std::ostream& write_demography(std::ostream&) const;
    //! Write #subpopulations_ and #juveniles_subpops_ for debugging.
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

  private:
    //! give birth to children
    void reproduce();

    //! Append new individuals to #juveniles_subpops_ at `location`.
    //! The expected number of children for each adult is proportional to its weight.
    //! All the females are evaluated for recruitment,
    //! whereas males are stochastically chosen.
    void reproduce(uint_fast32_t location, size_t popsize);

    //! Evaluate survival.
    void survive();

    //! Evaluate migration. Individuals in #juveniles_subpops_ moves to #subpopulations_.
    void migrate();

    //! Sample individuals.
    void sample(std::vector<std::vector<std::shared_ptr<Individual>>>* subpops,
                const std::vector<size_t>& sample_sizes);

    //! Append current state to #demography_
    void append_demography(int_fast32_t season);

    //! Count individuals for each location and age
    std::vector<std::vector<uint_fast32_t>> count(int_fast32_t season) const;

    //! Return size of #subpopulations_
    size_t num_subpops() const noexcept {return subpopulations_.size();}

    //! Individual array for each subpopulation
    std::vector<std::vector<std::shared_ptr<Individual>>> subpopulations_;
    //! First-year individuals separated for #sample().
    //! Note that it becomes empty in #migrate().
    std::vector<std::vector<std::shared_ptr<Individual>>> juveniles_subpops_;
    //! Counts of juveniles of the year: [[number for each location] for each season]
    std::vector<std::vector<uint_fast32_t>> juveniles_demography_ = {};
    //! Samples: [{capture_year => individuals} for each location]
    std::vector<std::map<int_fast32_t, std::vector<std::shared_ptr<Individual>>>> loc_year_samples_ = {};
    //! (year, season) => [[count for each age] for each location]
    std::map<std::pair<int_fast32_t, int_fast32_t>, std::vector<std::vector<uint_fast32_t>>> demography_ = {};

    //! @ingroup params
    //!@{

    //! \f$K\f$: carrying capacity used in reproduce()
    const double carrying_capacity_;
    //! \f$r\f$: coefficient used in reproduce()
    const double recruitment_coef_;
    //! \f$k \in (0, \infty)\f$ for overdispersion in reproduce().
    //! Equivalent to Poisson when \f$k \to \infty\f$ (or \f$k<0\f$ for convience).
    const double k_nbinom_;
    //!@}
    //! Current time.
    int_fast32_t year_ = 0;
    //! Uniform Random Bit Generator
    std::unique_ptr<URBG> engine_;
};

} // namespace pbf

#endif /* PBT_POPULATION_HPP_ */
