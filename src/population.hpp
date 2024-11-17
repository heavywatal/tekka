/*! @file population.hpp
    @brief Interface of Population class
*/
#pragma once
#ifndef PBT_POPULATION_HPP_
#define PBT_POPULATION_HPP_

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

enum class Sex { F, M };

class SubPopulation {
  public:
    std::vector<std::shared_ptr<Individual>>& operator[](Sex x) {
        return adults[static_cast<int>(x)];
    }
    const std::vector<std::shared_ptr<Individual>>& operator[](Sex x) const {
        return adults[static_cast<int>(x)];
    }
    size_t size() const {
        return adults[0].size() + adults[1].size();
    }
    void clear() {
        adults[0].clear();
        adults[1].clear();
    }
    //! Individual array for each subpopulation
    std::array<std::vector<std::shared_ptr<Individual>>, 2> adults{};
    //! First-year individuals separated for #sample().
    //! Note that it becomes empty in #migrate().
    std::vector<std::shared_ptr<Individual>> juveniles{};
    //! Samples: {capture_year => individuals}
    std::map<int_fast32_t, std::vector<std::shared_ptr<Individual>>> samples{};
    //! [[[count for each age] for season] for year]
    std::vector<std::array<std::vector<uint_fast32_t>, 4>> demography{};
};

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

    //! Append new individuals to #SubPopulation::juveniles at `location`.
    //! The expected number of children for each adult is proportional to its weight.
    //! All the females are evaluated for recruitment,
    //! whereas males are stochastically chosen.
    void reproduce(uint_fast32_t location, size_t popsize);

    //! Evaluate survival.
    void survive(int_fast32_t season);

    //! Evaluate migration. Merge #SubPopulation::juveniles to #SubPopulation::adults.
    void migrate();

    //! Sample individuals.
    void sample(SubPopulation& subpops, size_t n_adults, size_t n_juveniles);
    //! Implementation of sample().
    void sample(std::vector<std::shared_ptr<Individual>>& src,
                std::vector<std::shared_ptr<Individual>>& dst, size_t n);

    //! Initialize #SubPopulation::demography
    void init_demography(int_fast32_t duration);
    //! Record current state to #SubPopulation::demography
    void record_demography(int_fast32_t season);

    //! For male selection
    std::vector<double> weights(const std::vector<std::shared_ptr<Individual>>&) const;

    //!
    std::vector<SubPopulation> subpopulations_{};

    //! @ingroup params
    //!@{

    //! \f$K\f$: carrying capacity used in reproduce()
    const double carrying_capacity_{};
    //! \f$r\f$: coefficient used in reproduce()
    const double recruitment_coef_{};
    //! \f$k \in (0, \infty)\f$ for overdispersion in reproduce().
    //! Equivalent to Poisson when \f$k \to \infty\f$ (or \f$k<0\f$ for convience).
    const double k_nbinom_{};
    //!@}
    //! Current time.
    int_fast32_t year_{0};
    //! Uniform Random Bit Generator
    std::unique_ptr<URBG> engine_{nullptr};
};

} // namespace pbf

#endif /* PBT_POPULATION_HPP_ */
