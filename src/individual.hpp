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
#include <functional>
#include <limits>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbf {

//! @brief Parameters for Individual class (command-line)
/*! @ingroup params
*/
struct IndividualParams {
    //! @ingroup params
    //@{
    //! \f$r\f$:  used in Individual::recruitment()
    double RECRUITMENT_COEF = 2.0;
    //! \f$K\f$: carrying capacity used in Individual::recruitment()
    double CARRYING_CAPACITY = 1e+3;
    //! \f$k \in (0, \infty)\f$ for overdispersion in recruitment().
    //! Equivalent to Poisson when \f$k \to \infty\f$ (or \f$k<0\f$ for convience).
    double NEGATIVE_BINOM_K = -1.0;
    //@}
};

//! @brief Parameters for Individual class (JSON file)
/*! @ingroup params
*/
class IndividualJson {
    friend class Individual;
  public:
    //! @cond
    IndividualJson();
    //! @endcond
    //! Read class variables from stream in json
    void read(std::istream&);
    //! Write class variables to stream in json
    std::ostream& write(std::ostream&) const;
  private:
    //! @ingroup params
    //@{
    //! mortality due to natural causes
    std::vector<double> natural_mortality;
    //! mortality due to fishing activities
    std::vector<double> fishing_mortality;
    //! fluctuation of fishing mortality by year
    std::vector<double> fishing_coef;
    //! precalculated values (quarter age)
    std::vector<double> weight_for_age;
    //! transition matrix for migration
    std::vector<std::vector<std::vector<double>>> migration_matrices;
    //@}
};

/*! @brief Individual class
*/
class Individual {
  public:
    //! Alias
    using param_type = IndividualParams;
    Individual() = delete;
    //! for initial population
    explicit Individual(bool is_male): is_male_(is_male) {}
    //! for sexual reproduction
    Individual(const std::shared_ptr<Individual>& father,
               const std::shared_ptr<Individual>& mother, int_fast32_t year, bool is_male)
    : father_(father), mother_(mother),
      birth_year_(year), is_male_(is_male) {}

    //! finite death rate per year
    double death_rate(const int_fast32_t year) const {
        return death_rate(year - birth_year_, year);
    }

    //! number of juveniles
    uint_fast32_t recruitment(int_fast32_t year, double density_effect, URBG&) const noexcept;

    //! return new location
    uint_fast32_t migrate(uint_fast32_t loc, int_fast32_t year, URBG&);

    //! collect ancestral IDs
    void trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>* ids,
                    uint_fast32_t loc, int_fast32_t year) const;
    //! write all the data members in TSV
    std::ostream& write(std::ostream&) const;
    //! write all the data members in TSV with translated IDs
    std::ostream& write(std::ostream&, const std::unordered_map<const Individual*, uint_fast32_t>&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);
    //! column names for write()
    static std::vector<std::string> names();

    //! Parameters shared among instances (JSON file)
    static inline IndividualJson JSON;
    //! Set static variables that depend on others
    static void set_dependent_static(const uint_fast32_t years);
    static bool is_ready(const uint_fast32_t years) {
        return FISHING_COEF_.size() >= years;
    }
    //! finite death rate per year
    static double death_rate(const int_fast32_t age, const int_fast32_t year) {
        const auto f = FISHING_MORTALITY_[age] * FISHING_COEF_[year];
        return 1.0 - std::exp(-NATURAL_MORTALITY_[age] - f);
    }

    //! Set #PARAM_
    static void param(const param_type& p) {PARAM_ = p;}
    //! Get #PARAM_
    static const param_type& param() {return PARAM_;}

    //! @name Getter functions
    //@{
    //! weight for year age
    double weight(int_fast32_t year) const noexcept {
        return WEIGHT_FOR_AGE_[year - birth_year_];
    }
    //! !#father_
    bool is_first_gen() const noexcept {return !father_;}
    //! @cond
    const Individual* father_get() const noexcept {return father_.get();}
    const Individual* mother_get() const noexcept {return mother_.get();}
    int_fast32_t birth_year() const noexcept {return birth_year_;}
    bool is_male() const noexcept {return is_male_;}
    //! @endcond
    //@}

  private:
    static void set_static_migration();
    static void set_static_mortality();
    static void set_static_weight();
    constexpr static inline int_fast32_t MAX_AGE = 80;
    //! Parameters shared among instances (command-line)
    static inline param_type PARAM_;
    //! natural mortality per year
    static inline std::vector<double> NATURAL_MORTALITY_;
    //! fishing mortality per year
    static inline std::vector<double> FISHING_MORTALITY_;
    //! fluctuation of fishing mortality by year
    static inline std::vector<double> FISHING_COEF_;
    //! precalculated values (age)
    static inline std::vector<double> WEIGHT_FOR_AGE_;
    //! discrete distributions for migration
    static inline std::vector<std::vector<std::function<uint_fast32_t(URBG&)>>> MIGRATION_DISTRIBUTIONS_;

    //! father
    const std::shared_ptr<Individual> father_ = nullptr;
    //! mother
    const std::shared_ptr<Individual> mother_ = nullptr;
    //! year of birth
    int_fast32_t birth_year_ = -4;
    //! sex
    const bool is_male_;
};

} // namespace pbf

#endif /* PBT_INDIVIDUAL_HPP_ */
