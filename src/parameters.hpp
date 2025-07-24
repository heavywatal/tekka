#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace pbf {

//! @defgroup parameters Parameters
//! List of commandline options.
//! Tables are at the bottom of the page.
//!
//! The longer ones can be used as keys in JSON config file,
//! e.g., `{"years": 80}` for `-y,--years` option.

//! @brief Parameters for Population class
struct Parameters {
    //! Constructor
    Parameters();
    //! Read json-only variables from stream
    void read(std::istream&);
    //! Write json-only variables to stream
    std::ostream& write(std::ostream&) const;

//! @addtogroup parameters
//! @{

    //! @name Global
    //!@{

    //! RNG seed; 32-bit signed integer for R
    int32_t seed{0};
    //! Output directory
    std::string outdir{};
    //! Duration of simulation
    int_fast32_t years{80};
    //! Initial number of juveniles at location 0 in season 0 in year 0.
    //! It falls back to \f$\operatorname{med}(R)\f$ or \f$K\f$ if 0.
    int_fast32_t origin{0};
    //!@}

    //! @name Reproduction
    //!@{

    //! \f$\operatorname{med}(R) = \exp(\mu_{R})\f$ for lognormal distribution
    double med_recruitment{0.0};
    //! \f$\sigma_{R}\f$ for lognormal distribution
    double sigma_recruitment{1.0};
    //! \f$K\f$: carrying capacity used in Population::reproduce()
    double carrying_capacity{2e3};
    //! \f$r\f$: coefficient used in Population::reproduce()
    double recruitment{2.0};
    //! \f$k \in (0, \infty)\f$ for overdispersion in Population::reproduce().
    //! Equivalent to Poisson when \f$k \to \infty\f$ (or \f$k<0\f$ for convience).
    double overdispersion{-1.0};
    //!@}

    //! @name Sampling
    //!@{

    //! Sample last _ years
    int_fast32_t last{3};
    //! per location
    std::vector<int_fast32_t> sample_size_adult{10, 10};
    //! per location
    std::vector<int_fast32_t> sample_size_juvenile{10, 10};
    //!@}

    //! @name Configurable only via JSON
    //!@{

    //! Alias for readability
    using RowMatrix = std::vector<std::vector<double>>;
    //! Array of \f$M\f$ for quarter age: instantaneous mortality due to natural causes
    std::vector<double> natural_mortality{};
    //! Array of \f$F\f$ for quarter age: instantaneous mortality due to fishing activities
    std::vector<double> fishing_mortality{};
    //! Array of \f$e\f$ by year: coefficient of fishing mortality.
    //! The last element is used for the last years if its length differs from #years.
    std::vector<double> fishing_coef{};
    //! Weight in kg for quarter age
    std::vector<double> weight_for_age{};
    //! Transition matrix for Population::migrate()
    std::vector<RowMatrix> migration_matrices{};
    //!@}

//!@}

};

} // namespace pbf
