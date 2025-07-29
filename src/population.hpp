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

//! @cond
namespace pcglite {
  template <class UIntType> class permuted_congruential_engine;

  using pcg64 = permuted_congruential_engine<uint64_t>;
}
//! @endcond

namespace pbf {

//! @cond

//! Type alias of random number generator
using URBG = pcglite::pcg64;

class Individual;
using ShPtrIndividual = std::shared_ptr<Individual>;

enum class Sex { F, M };

//! @endcond

//! Component of Population; Container of Individual
class SubPopulation {
  public:
    //! Getter of #adults.
    std::vector<ShPtrIndividual>& operator[](Sex x) {
        return adults[static_cast<uint_fast8_t>(x)];
    }
    //! Getter of #adults.
    const std::vector<ShPtrIndividual>& operator[](Sex x) const {
        return adults[static_cast<uint_fast8_t>(x)];
    }
    //! The number of #adults.
    size_t size() const {
        return adults[0].size() + adults[1].size();
    }
    //! Remove all #adults.
    void clear() {
        adults[0].clear();
        adults[1].clear();
    }
    //! Array of shared pointers to Individual's.
    //! Females and males are separately stored for #Population::reproduce().
    std::array<std::vector<ShPtrIndividual>, 2> adults{};
    //! First-year individuals separated for #Population::sample().
    std::vector<ShPtrIndividual> juveniles{};
    //! Samples: {capture_year => individuals}
    std::map<int_fast32_t, std::vector<ShPtrIndividual>> samples{};
    //! [[[count for each age] for season] for year]
    std::vector<std::array<std::vector<int_fast32_t>, 4>> demography{};
};

//! Main class that implements simulation
class Population {
    //! No individual lives longer than this.
    constexpr static inline int_fast32_t MAX_AGE = 80;

  public:
    //! Constructor
    Population(const Parameters&);
    ~Population();

    //! The main loop of simulation.
    //! 1. reproduce()
    //! 2. record_demography() and survive() for each season
    //! 3. sample()
    //! 4. migrate()
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
    //! Call reproduce_impl() via reproduce_lognormal() or reproduce_logistic()
    //! depending on #Parameters::med_recruitment.
    //! Currently, the breeding places are hardcoded to subpopulations 0 and 1.
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
    //! All the #SubPopulation::adults are evaluated regardless of age.
    void reproduce_impl(SubPopulation&, const std::vector<int_fast32_t>& litter_sizes);
    //! Calculate stochastic litter sizes proportional to female weight
    std::vector<int_fast32_t>
    litter_sizes_logistic(const std::vector<ShPtrIndividual>& females, int_fast32_t popsize);

    //! Evaluate survival using death_rate().
    void survive(int_fast32_t season);
    //! Finite death rate per quarter year: \f$ d = 1 - \exp(- M - eF) \f$
    double death_rate(const int_fast32_t age, const int_fast32_t season) const;

    //! Evaluate migration. Merge #SubPopulation::juveniles to #SubPopulation::adults.
    void migrate();
    //! Generate random number for new location.
    int_fast32_t destination(int_fast32_t age, int_fast32_t loc);

    //! Move individuals to SubPopulation::samples.
    //! "Juveniles" are sampled from #SubPopulation::juveniles,
    //! first-year individuals before migration.
    void sample(SubPopulation& subpops, int_fast32_t n_adults, int_fast32_t n_juveniles);
    //! Implementation of sample().
    int_fast32_t sample(std::vector<ShPtrIndividual>& src,
                        std::vector<ShPtrIndividual>& dst, int_fast32_t n);

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

    //! Store the values from command line and json.
    Parameters params_;

    //! elongated version of #Parameters::natural_mortality
    std::vector<double> NATURAL_MORTALITY_{};
    //! elongated version of #Parameters::fishing_mortality
    std::vector<double> FISHING_MORTALITY_{};
    //! elongated version of #Parameters::fishing_coef
    std::vector<double> FISHING_COEF_{};
    //! year version of #Parameters::weight_for_age
    std::vector<double> WEIGHT_FOR_AGE_{};
    //! Alias for readability
    using PairDestDist=std::pair<int_fast32_t, std::discrete_distribution<int_fast32_t>>;
    //! Discrete distributions for #migrate() based on #Parameters::migration_matrices.
    std::vector<std::vector<PairDestDist>> MIGRATION_DESTINATION_{};

    //! Array of SubPopulation.
    std::vector<SubPopulation> subpopulations_{};
    //! Current time.
    int_fast32_t year_{0};
    //! Uniform Random Bit Generator
    std::unique_ptr<URBG> engine_{nullptr};
};

} // namespace pbf

#endif /* PBT_POPULATION_HPP_ */
