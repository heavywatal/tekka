/*! @file individual.hpp
    @brief Interface of Individual class
*/
#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include <iosfwd>
#include <cstdint>
#include <random>

#include <boost/program_options.hpp>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

/*! @brief Individual class
*/
class Individual {
  public:
    //! default constructor
    Individual()
    : id_(++LAST_ID_) {}
    //! sexual reproduction
    Individual(const Individual& father, const Individual& mother, const uint_fast32_t time)
    : id_(++LAST_ID_), father_id_(father.id()), mother_id_(mother.id()), birth_date_(time) {}

    //! evaluate survival
    bool has_survived() const;

    //! evaluate maturity
    bool is_matured(const uint_fast32_t time) const {
        return time > birth_date_ + AGE_OF_MATURATION_;
    }

    //! random number of eggs
    template <class URBG>
    uint_fast32_t clutch_size(URBG& g) const {return POISSON_CLUTCH_SIZE_(g);}

    //! write
    std::ostream& write(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream&, const Individual&);
    //! options description for Individual class
    static boost::program_options::options_description options_desc();
    //! initialize static members after po::notify(vm);
    static void init_static_members();
    //! unit test
    static void test();

    //! getter of #id_
    uint_fast32_t id() const {return id_;}
  private:
    //! by the quater-year
    constexpr static uint_fast32_t AGE_OF_MATURATION_ = 3 * 4;
    //! parameter for #POISSON_CLUTCH_SIZE_
    static uint_fast32_t MEAN_CLUTCH_SIZE_;
    //! distribution for clutch_size()
    static std::poisson_distribution<uint_fast32_t> POISSON_CLUTCH_SIZE_;
    //! mortality due to natural causes
    static double NATURAL_MORTALITY_;
    //! mortality due to fishing activities
    static double FISHING_MORTALITY_;
    //! ID for a new instance
    static uint_fast32_t LAST_ID_;

    //! ID
    const uint_fast32_t id_ = 0;
    //! father's ID
    const uint_fast32_t father_id_ = 0;
    //! mother's ID
    const uint_fast32_t mother_id_ = 0;
    //! by the quater-year
    uint_fast32_t birth_date_ = 0;
};

} // namespace pbt

#endif /* PBT_INDIVIDUAL_HPP_ */
