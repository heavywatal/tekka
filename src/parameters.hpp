#pragma once

#include <iosfwd>
#include <vector>

namespace pbf {

//! @brief Parameters for Individual class (JSON file)
/*! @ingroup params
*/
class Parameters {
    friend class Individual;
    friend class Population;
  public:
    //! Constructor
    Parameters();
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
