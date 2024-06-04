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
struct IndividualJson {
    //! @cond
    IndividualJson();
    void set_dependent_static();
    void read(std::istream&);
    void write(std::ostream&) const;
    //! @endcond

    //! @ingroup params
    //@{
    //! mortality due to natural causes
    std::vector<double> NATURAL_MORTALITY;
    //! mortality due to fishing activities
    std::vector<double> FISHING_MORTALITY;
    //! precalculated values (quater age)
    std::vector<double> WEIGHT_FOR_AGE;
    //! transition matrix for migration
    std::vector<std::vector<std::vector<double>>> MIGRATION_MATRICES;
    //@}

    //! finite death rate per year
    std::vector<double> DEATH_RATE;
    //! precalculated values (age)
    std::vector<double> WEIGHT_FOR_YEAR_AGE;
    //! discrete distributions for migration
    std::vector<std::vector<std::function<uint_fast32_t(URBG&)>>> MIGRATION_DISTRIBUTIONS;
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

    //! evaluate survival
    bool is_dead(const int_fast32_t year, URBG&) const;

    //! number of juveniles
    uint_fast32_t recruitment(int_fast32_t year, double density_effect, URBG&) const noexcept;

    //! return new location
    uint_fast32_t migrate(uint_fast32_t loc, int_fast32_t year, URBG&);

    //! collect ancestoral IDs
    void trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>* ids,
                    uint_fast32_t loc, int_fast32_t year) const;
    //! write all the data members in TSV
    std::ostream& write(std::ostream&) const;
    //! write all the data members in TSV with translated IDs
    std::ostream& write(std::ostream&, const std::unordered_map<const Individual*, uint_fast32_t>&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);
    //! column names for write()
    static std::vector<std::string> names();
    //! Read class variables from stream in json
    static void read_json(std::istream& ist) {JSON_.read(ist);}
    //! Write class variables to stream in json
    static void write_json(std::ostream& ost) {JSON_.write(ost);}
    //! Set #PARAM_
    static void param(const param_type& p) {PARAM_ = p;}
    //! Get #PARAM_
    static const param_type& param() {return PARAM_;}

    //! @name Getter functions
    //@{
    //! IndividualJson.DEATH_RATE
    static const std::vector<double>&
    death_rate() {return JSON_.DEATH_RATE;}
    //! IndividualJson.WEIGHT_FOR_YEAR_AGE
    double weight(int_fast32_t year) const noexcept {
        return JSON_.WEIGHT_FOR_YEAR_AGE[year - birth_year_];
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
    //! Parameters shared among instances (command-line)
    static param_type PARAM_;
    //! Parameters shared among instances (JSON file)
    static IndividualJson JSON_;

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
