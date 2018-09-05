/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <vector>
#include <map>

namespace wtl {class sfmt19937_64;}

namespace boost {namespace program_options {class options_description;}}

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

//! alias of uniform random bit generator
using URBG = wtl::sfmt19937_64;

/*! @brief Individual class
*/
class Individual {
  public:
    //! for initial population
    Individual() = default;
    //! for sexual reproduction
    Individual(const std::shared_ptr<Individual>& father,
               const std::shared_ptr<Individual>& mother, uint_fast32_t year)
    : father_(father), mother_(mother),
      birth_year_(year), location_(mother ? mother->location() : 0u) {}

    //! evaluate survival
    bool has_survived(const uint_fast32_t year, const uint_fast32_t quarter, URBG&) const;

    //! evaluate maturity
    bool is_in_breeding_place() const noexcept {
        return location_ < num_breeding_places();
    }

    //! number of juveniles
    uint_fast32_t recruitment(const uint_fast32_t year, URBG&) const noexcept;

    //! change #location_
    void migrate(const uint_fast32_t year, URBG&);

    //! compare Individual pointer
    struct less {
        //! compare Individual pointer
        bool operator()(const Individual* const px, const Individual* const py) const noexcept {
            return px < py;
        }
    };

    //! collect ancestoral IDs
    void trace_back(std::map<const Individual*, uint_fast32_t, less>* nodes, uint_fast32_t year=0u) const noexcept {
        if (nodes->emplace(this, year).second && father_) {
            father_->trace_back(nodes);
            mother_->trace_back(nodes);
        }
    }

    //! write all the data members in TSV
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);
    //! column names for write()
    static std::vector<std::string> names();
    //! options description for Individual class
    static boost::program_options::options_description options_desc();
    //! set static variables from config.hpp
    static void set_default_values();
    //! Read class variables from stream in json
    static void read_json(std::istream&);
    //! Write class variables to stream in json
    static void write_json(std::ostream&);

    //! number of locations
    static size_t num_locations() noexcept {
        return MIGRATION_MATRICES_[0u].size();
    }
    //! number of breeding places
    static constexpr size_t num_breeding_places() noexcept {return 2u;}

    //! @name Getter functions
    //@{
    double weight(const uint_fast32_t year) const noexcept {
        return WEIGHT_FOR_AGE_[4u * (year - birth_year_)];
    }
    bool is_first_gen() const noexcept {return !father_;}
    const Individual* father_get() const noexcept {return father_.get();}
    const Individual* mother_get() const noexcept {return mother_.get();}
    uint_fast32_t birth_year() const noexcept {return birth_year_;}
    uint_fast32_t location() const noexcept {return location_;}
    //@}

  private:
    //! set static variables that depend on other variables
    static void set_dependent_static();
    //! parameter for recruitment()
    //! @ingroup params
    static double RECRUITMENT_COEF_;
    //! \f$k\f$ for overdispersion in recruitment()
    //! @ingroup params
    static double NEGATIVE_BINOM_K_;
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

    //! father
    const std::shared_ptr<Individual> father_ = nullptr;
    //! mother
    const std::shared_ptr<Individual> mother_ = nullptr;
    //! year of birth
    uint_fast32_t birth_year_ = 0u;
    //! current location
    uint_fast32_t location_ = 0u;
};

} // namespace pbt

#endif /* PBT_INDIVIDUAL_HPP_ */
