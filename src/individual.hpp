/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include "parameters.hpp"

#include <cstdint>
#include <iosfwd>
#include <limits>
#include <memory>
#include <vector>
#include <random>
#include <unordered_map>

namespace pbf {

/*! @brief Individual class
*/
class Individual {
  public:
    //! for initial population
    Individual() = default;
    //! for sexual reproduction
    Individual(std::shared_ptr<Individual> father,
               std::shared_ptr<Individual> mother,
               int_fast32_t year) noexcept
    : father_(std::move(father)),
      mother_(std::move(mother)),
      birth_year_(year) {}
    Individual(const Individual&) = delete;
    Individual(Individual&&) = delete;
    ~Individual() = default;

    //! Finite death rate per quarter year: \f$ d = 1 - \exp(- M - eF) \f$
    double death_rate(const int_fast32_t year, const int_fast32_t season) const {
        const auto q_age = 4 * age(year) + season;
        return DEATH_RATE(q_age, year);
    }

    //! Generate random number for new location.
    template <class URBG>
    uint_fast32_t destination(uint_fast32_t loc, int_fast32_t year, URBG& engine) const {
        auto& [dest, dist] = MIGRATION_DESTINATION_[age(year)][loc];
        return (dest == std::numeric_limits<uint_fast32_t>::max()) ? dist(engine) : dest;
    }

    //! Write ancestors recursively.
    void trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>& ids,
                    uint_fast32_t loc, int_fast32_t year) const;
    //! Write all the data members in TSV
    std::ostream& write(std::ostream&) const;
    //! Write all the data members in TSV with translated IDs
    std::ostream& write(std::ostream&, const std::unordered_map<const Individual*, uint_fast32_t>&) const;
    //! Write column names for trace_back()
    static std::ostream& write_trace_back_header(std::ostream&);
    //! Write column names for write()
    static std::ostream& write_names(std::ostream&);
    friend std::ostream& operator<<(std::ostream&, const Individual&);

    //! No individual lives longer than this.
    constexpr static inline int_fast32_t MAX_AGE = 80;
    //! Set static variables that depend on others.
    static void set_dependent_static(const Parameters&, const uint_fast32_t years);
    //! Test if static variables are ready.
    static bool is_ready(const uint_fast32_t years) {
        return (
          FISHING_COEF_.size() >= years
          && NATURAL_MORTALITY_.size() >= 4u * MAX_AGE
          && FISHING_MORTALITY_.size() >= 4u * MAX_AGE
          && WEIGHT_FOR_AGE_.size() >= MAX_AGE
          && MIGRATION_DESTINATION_.size() >= MAX_AGE
        );
    }
    //! Finite death rate per quarter year
    static double DEATH_RATE(const int_fast32_t q_age, const int_fast32_t year) {
        const auto m = NATURAL_MORTALITY_[q_age];
        const auto f = FISHING_MORTALITY_[q_age] * FISHING_COEF_[year];
        return 1.0 - std::exp(-m - f);
    }

    //! @name Getter functions
    //!@{

    //! Age in the given year
    int_fast32_t age(const int_fast32_t year) const noexcept {
        return year - birth_year_;
    }
    //! Weight for year age
    double weight(int_fast32_t year) const noexcept {
        return WEIGHT_FOR_AGE_[age(year)];
    }
    //!@}

  private:
    //! Alias for readability
    using PairDestDist=std::pair<uint_fast32_t, std::discrete_distribution<uint_fast32_t>>;
    //! Prepare #MIGRATION_DESTINATION_
    static void set_static_migration(const Parameters&);
    //! Prepare mortality vectors
    static void set_static_mortality(const Parameters&);
    //! Prepare #WEIGHT_FOR_AGE_
    static void set_static_weight(const Parameters&);
    static inline std::vector<double> NATURAL_MORTALITY_{};
    static inline std::vector<double> FISHING_MORTALITY_{};
    //! Fluctuation of fishing mortality by year
    static inline std::vector<double> FISHING_COEF_{};
    //! year version of #Parameters::weight_for_age
    static inline std::vector<double> WEIGHT_FOR_AGE_{};
    //! Discrete distributions for migration
    static inline std::vector<std::vector<PairDestDist>> MIGRATION_DESTINATION_{};

    //! Pointer to father
    const std::shared_ptr<Individual> father_{nullptr};
    //! Pointer to mother
    const std::shared_ptr<Individual> mother_{nullptr};
    //! Year of birth
    int_fast32_t birth_year_{-4};
};

} // namespace pbf

#endif /* PBT_INDIVIDUAL_HPP_ */
