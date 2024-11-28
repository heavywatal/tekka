#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace pbf {

//! @brief Parameters for Population class
/*! @ingroup params
*/
struct Parameters {
    //! Constructor
    Parameters();
    //! Read json-only variables from stream
    void read(std::istream&);
    //! Write json-only variables to stream
    std::ostream& write(std::ostream&) const;

    //! Initial population size relative to K
    double origin{0.2};
    //! Duration of simulation
    int years{100};
    //! Sample last _ years
    int last{3};
    //! per location
    std::vector<size_t> sample_size_adult{10u, 10u};
    //! per location
    std::vector<size_t> sample_size_juvenile{10u, 10u};
    //! Output directory
    std::string outdir{};
    //! RNG seed; 32-bit signed integer for R
    int32_t seed{0};
    //! \f$K\f$: carrying capacity used in Population::reproduce()
    double carrying_capacity{1e3};
    //! \f$r\f$: coefficient used in Population::reproduce()
    double recruitment{2.0};
    //! \f$k \in (0, \infty)\f$ for overdispersion in Population::reproduce().
    //! Equivalent to Poisson when \f$k \to \infty\f$ (or \f$k<0\f$ for convience).
    double overdispersion{-1.0};

    //! @ingroup only-json
    //!@{

    //! Alias for readability
    using RowMatrix = std::vector<std::vector<double>>;
    //! Array of \f$M\f$ for quarter age: instantaneous mortality due to natural causes
    std::vector<double> natural_mortality{};
    //! Array of \f$F\f$ for quarter age: instantaneous mortality due to fishing activities
    std::vector<double> fishing_mortality{};
    //! Array of \f$e\f$ by year: coefficient of fishing mortality.
    //! The last part is used for the last years if its length differs from `--years` option.
    std::vector<double> fishing_coef{};
    //! Weight in kg for quarter age
    std::vector<double> weight_for_age{};
    //! Transition matrix for migration
    std::vector<RowMatrix> migration_matrices{};
    //!@}
};

} // namespace pbf
