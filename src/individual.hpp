/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include "random_fwd.hpp"

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <vector>
#include <unordered_map>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbf {

//! @brief Parameters for Individual class (JSON file)
/*! @ingroup params
*/
class IndividualJson {
    friend class Individual;
  public:
    //! Constructor
    IndividualJson();
    //! Read class variables from stream in json
    void read(std::istream&);
    //! Write class variables to stream in json
    std::ostream& write(std::ostream&) const;
  private:
    //! Alias for readability
    using RowMatrix = std::vector<std::vector<double>>;
    //! @ingroup params
    //!@{

    //! Array of \f$M\f$ for quarter age: instantaneous mortality due to natural causes
    std::vector<double> natural_mortality = {};
    //! Array of \f$F\f$ for quarter age: instantaneous mortality due to fishing activities
    std::vector<double> fishing_mortality = {};
    //! Array of \f$e\f$ by year: coefficient of fishing mortality.
    //! The last part is used for the last years if its length differs from `--years` option.
    std::vector<double> fishing_coef = {};
    //! Weight in kg for quarter age
    std::vector<double> weight_for_age = {};
    //! Transition matrix for migration
    std::vector<RowMatrix> migration_matrices = {};
    //!@}
};

/*! @brief Individual class
*/
class Individual {
  public:
    //! Alias
    Individual() = delete;
    //! for initial population
    explicit Individual(bool is_male) noexcept: is_male_(is_male) {}
    //! for sexual reproduction
    Individual(const std::shared_ptr<Individual>& father,
               const std::shared_ptr<Individual>& mother,
               int_fast32_t year, bool is_male) noexcept
    : father_(father), mother_(mother),
      birth_year_(year), is_male_(is_male) {}
    Individual(const Individual&) = delete;
    ~Individual() = default;

    //! Finite death rate per quarter year: \f$ d = 1 - \exp(- M - eF) \f$
    double death_rate(const int_fast32_t year, const int_fast32_t season) const {
        const auto q_age = 4 * age(year) + season;
        return DEATH_RATE(q_age, year);
    }

    //! Generate random number for new location.
    int_fast32_t migrate(uint_fast32_t loc, int_fast32_t year, URBG&);

    //! Write ancestors recursively.
    void trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>* ids,
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
    //! Parameters shared among instances (JSON file).
    static inline IndividualJson JSON;
    //! Set static variables that depend on others.
    static void set_dependent_static(const uint_fast32_t years);
    //! Test if static variables are ready.
    static bool is_ready(const uint_fast32_t years) {
        return (
          FISHING_COEF_.size() >= years
          && JSON.natural_mortality.size() >= 4u * MAX_AGE
          && JSON.fishing_mortality.size() >= 4u * MAX_AGE
          && WEIGHT_FOR_AGE_.size() >= MAX_AGE
          && MIGRATION_DESTINATION_.size() >= MAX_AGE
        );
    }
    //! Finite death rate per quarter year
    static double DEATH_RATE(const int_fast32_t q_age, const int_fast32_t year) {
        const auto m = JSON.natural_mortality[q_age];
        const auto f = JSON.fishing_mortality[q_age] * FISHING_COEF_[year];
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
    //! Just returns #is_male_
    bool is_male() const noexcept {return is_male_;}
    //!@}

  private:
    //! Alias for readability
    using PairDestDist=std::pair<int_fast32_t, std::discrete_distribution<int_fast32_t>>;
    //! Prepare #MIGRATION_DESTINATION_
    static void set_static_migration();
    //! Prepare mortality vectors
    static void set_static_mortality();
    //! Prepare #WEIGHT_FOR_AGE_
    static void set_static_weight();
    //! Fluctuation of fishing mortality by year
    static inline std::vector<double> FISHING_COEF_;
    //! year version of #IndividualJson::weight_for_age
    static inline std::vector<double> WEIGHT_FOR_AGE_;
    //! Discrete distributions for migration
    static inline std::vector<std::vector<PairDestDist>> MIGRATION_DESTINATION_;

    //! Pointer to father
    const std::shared_ptr<Individual> father_ = nullptr;
    //! Pointer to mother
    const std::shared_ptr<Individual> mother_ = nullptr;
    //! Year of birth
    int_fast32_t birth_year_ = -4;
    //! Sex
    const bool is_male_;
};

} // namespace pbf

#endif /* PBT_INDIVIDUAL_HPP_ */
