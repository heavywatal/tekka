/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include "common.hpp"

#include <iosfwd>
#include <cstdint>
#include <cmath>
#include <random>
#include <memory>

#include <boost/program_options.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

/*! @brief Individual class
*/
class Individual {
  public:
    //! default constructor
    Individual() = default;
    //! sexual reproduction
    Individual(const std::shared_ptr<Individual>& father, const std::shared_ptr<Individual>& mother, const uint_fast32_t year)
    : father_(father), mother_(mother), id_(++LAST_ID_),
      birth_year_(year), location_(mother->location()) {}

    //! evaluate survival
    bool has_survived(const uint_fast32_t year, const uint_fast32_t quarter, urbg_t&) const;

    //! evaluate maturity
    bool is_in_breeding_place() const {
        return location_ < 2u;
    }

    //! number of juveniles
    uint_fast32_t recruitment(const uint_fast32_t year, urbg_t& g) const {
        const double mean = RECRUITMENT_COEF_ * weight(year);
        std::poisson_distribution<uint_fast32_t> poisson(mean);
        return poisson(g);
    }

    //! change #location_
    void migrate(const uint_fast32_t year, urbg_t& g) {
        location_ = MIGRATION_DISTRIBUTIONS_[year - birth_year_][location_](g);
    }

    void trace_back(std::map<uint_fast32_t, Individual*>* nodes) {
        if (nodes->emplace(id_, this).second && father_) {
            father_->trace_back(nodes);
            mother_->trace_back(nodes);
        }
    }

    //! write all the data members in TSV
    std::ostream& write(std::ostream&) const;
    //! write ids in TSV
    std::ostream& write_ids(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);
    //! column names for write()
    static std::vector<std::string> names();
    //! options description for Individual class
    static boost::program_options::options_description options_desc();
    //! set static variables from config.hpp
    static void set_default_values();
    //! set class variables from json
    static void from_json(const json::json&);
    //! encode class variables to json
    static void to_json(json::json&);

    //! number of locations
    static size_t num_locations() {
        return MIGRATION_MATRICES_[0].size();
    }
    //! getter of #WEIGHT_FOR_AGE_
    double weight(const uint_fast32_t year) const {
        return WEIGHT_FOR_AGE_[4u * (year - birth_year_)];
    }
    //! getter of #id_
    uint_fast32_t id() const {return id_;}
    //! getter of #location_
    uint_fast32_t location() const {return location_;}

  private:
    //! set static variables that depend on other variables
    static void set_dependent_static();
    //! maximum age to consider
    constexpr static uint_fast32_t MAX_AGE_ = 80u;
    //! \f$K\f$ in weight()
    constexpr static double GROWTH_RATE_ = 0.08;
    //! \f$L\f$ in weight()
    constexpr static double MAX_WEIGHT_ = 500.0;
    //! parameter for recruitment()
    static double RECRUITMENT_COEF_;
    //! mortality due to natural causes
    static std::vector<double> NATURAL_MORTALITY_;
    //! mortality due to fishing activities
    static std::vector<double> FISHING_MORTALITY_;
    //! survival rate per quater year
    static std::vector<double> SURVIVAL_RATE_;
    //! precalculated values
    static std::vector<double> WEIGHT_FOR_AGE_;
    //! transition matrix for migration
    static std::vector<std::vector<std::vector<double>>> MIGRATION_MATRICES_;
    //! discrete distributions for migration
    static std::vector<std::vector<std::discrete_distribution<uint_fast32_t>>> MIGRATION_DISTRIBUTIONS_;
    //! ID for a new instance
    static uint_fast32_t LAST_ID_;

    //! father
    const std::shared_ptr<Individual> father_ = nullptr;
    //! mother
    const std::shared_ptr<Individual> mother_ = nullptr;
    //! ID
    const uint_fast32_t id_ = 0;
    //! year of birth
    uint_fast32_t birth_year_ = 0;
    //! current location
    uint_fast32_t location_ = 0;
};

} // namespace pbt

#endif /* PBT_INDIVIDUAL_HPP_ */
