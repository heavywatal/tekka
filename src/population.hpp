/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

#include "parameters.hpp"

#include <array>
#include <cstdint>
#include <iosfwd>
#include <vector>
#include <map>
#include <memory>
#include <random>

namespace pcglite {
  template <class UIntType> class permuted_congruential_engine;
  using pcg64 = permuted_congruential_engine<uint64_t>;
}

namespace pbf {

using URBG = pcglite::pcg64;

class Individual;
using ShPtrIndividual = std::shared_ptr<Individual>;

enum class Sex { F, M };

class SubPopulation {
  public:
    std::vector<ShPtrIndividual>& operator[](Sex x) {
        return adults[static_cast<uint_fast8_t>(x)];
    }
    const std::vector<ShPtrIndividual>& operator[](Sex x) const {
        return adults[static_cast<uint_fast8_t>(x)];
    }
    size_t size() const {
        return adults[0].size() + adults[1].size();
    }
    void clear() {
        adults[0].clear();
        adults[1].clear();
    }
    //! Individual array for each subpopulation
    std::array<std::vector<ShPtrIndividual>, 2> adults{};
    //! First-year individuals separated for #sample().
    //! Note that it becomes empty in #migrate().
    std::vector<ShPtrIndividual> juveniles{};
    //! Samples: {capture_year => individuals}
    std::map<int_fast32_t, std::vector<ShPtrIndividual>> samples{};
    //! [[[count for each age] for season] for year]
    std::vector<std::array<std::vector<uint_fast32_t>, 4>> demography{};
};

/*! @brief Population class
*/
class Population {
    //! No individual lives longer than this.
    constexpr static inline int_fast32_t MAX_AGE = 80;

  public:
    //! Constructor
    Population(const Parameters&);
    ~Population();

    //! Main iteration
    void run();

    //! Construct and write tree from #SubPopulation::samples.
    //! Serial IDs are assigned to sampled individuals.
    //! IDs larger than sample size are assigned to unsampled ancestors.
    std::ostream& write_sample_family(std::ostream& ost) const;
    //! Write #SubPopulation::demography.
    std::ostream& write_demography(std::ostream&) const;
    //! Write #subpopulations_ for debugging.
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Population&);

  private:
    //! give birth to children
    void reproduce();
    //! Reproduce according to constant lognormal distribution.
    //! Total recruitment is calculated first, and then split it.
    void reproduce_lognormal();
    //! Reproduce according to population size and carrying capacity.
    void reproduce_logistic();
    //! Append new individuals to #SubPopulation::juveniles.
    //! The expected number of children for each adult is proportional to its weight.
    //! All the females are evaluated for recruitment,
    //! whereas males are stochastically chosen.
    void reproduce_impl(SubPopulation&, const std::vector<uint_fast32_t>& litter_sizes);
    //! Calculate stochastic litter sizes proportional to female weight
    std::vector<uint_fast32_t>
    litter_sizes_logistic(const std::vector<ShPtrIndividual>& females, size_t popsize);

    //! Evaluate survival.
    void survive(int_fast32_t season);
    //! Finite death rate per quarter year: \f$ d = 1 - \exp(- M - eF) \f$
    double death_rate(const int_fast32_t age, const int_fast32_t season) const;

    //! Evaluate migration. Merge #SubPopulation::juveniles to #SubPopulation::adults.
    void migrate();
    //! Generate random number for new location.
    uint_fast32_t destination(int_fast32_t age, uint_fast32_t loc);

    //! Sample individuals.
    void sample(SubPopulation& subpops, size_t n_adults, size_t n_juveniles);
    //! Implementation of sample().
    void sample(std::vector<ShPtrIndividual>& src,
                std::vector<ShPtrIndividual>& dst, size_t n);

    //! Initialize #SubPopulation::demography
    void init_demography(int_fast32_t duration);
    //! Record current state to #SubPopulation::demography
    void record_demography(int_fast32_t season);

    //! For male selection
    std::vector<double> weights(const std::vector<ShPtrIndividual>&) const;

    //! Call all the init* functions
    void propagate_params();
    //! Prepare #MIGRATION_DESTINATION_
    void init_migration();
    //! Prepare #NATURAL_MORTALITY_, #FISHING_MORTALITY_, #FISHING_COEF_
    void init_mortality();
    //! Prepare #WEIGHT_FOR_AGE_
    void init_weight();
    //! Test if dependent variables are ready.
    bool is_ready() const;

    //! @ingroup params
    //!@{

    Parameters params_;

    //! elongated version of #Parameters::natural_mortality
    std::vector<double> NATURAL_MORTALITY_{};
    //! elongated version of #Parameters::fishing_mortality
    std::vector<double> FISHING_MORTALITY_{};
    //! Fluctuation of fishing mortality by year
    std::vector<double> FISHING_COEF_{};
    //! year version of #Parameters::weight_for_age
    std::vector<double> WEIGHT_FOR_AGE_{};
    //! Alias for readability
    using PairDestDist=std::pair<uint_fast32_t, std::discrete_distribution<uint_fast32_t>>;
    //! Discrete distributions for migration
    std::vector<std::vector<PairDestDist>> MIGRATION_DESTINATION_{};

    //!@}

    //!
    std::vector<SubPopulation> subpopulations_{};
    //! Current time.
    int_fast32_t year_{0};
    //! Uniform Random Bit Generator
    std::unique_ptr<URBG> engine_{nullptr};
};

} // namespace pbf

#endif /* PBT_POPULATION_HPP_ */
