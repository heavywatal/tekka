#pragma once
#ifndef PBT_INDIVIDUAL_HPP_
#define PBT_INDIVIDUAL_HPP_

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <unordered_map>

namespace pbf {

/*! @brief Component of SubPopulation
*/
class Individual {
  public:
    //! The origin without parents for Population initialization
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

    //! Write ancestors recursively.
    void trace_back(std::ostream& ost, std::unordered_map<const Individual*, int_fast32_t>& ids,
                    int_fast32_t loc, int_fast32_t year) const;
    //! Write all the data members in TSV
    std::ostream& write(std::ostream&) const;
    //! Write all the data members in TSV with translated IDs
    std::ostream& write(std::ostream&, const std::unordered_map<const Individual*, int_fast32_t>&) const;
    //! Write column names for trace_back()
    static std::ostream& write_trace_back_header(std::ostream&);
    //! Write column names for write()
    static std::ostream& write_names(std::ostream&);
    friend std::ostream& operator<<(std::ostream&, const Individual&);

    //! Age in the given year
    int_fast32_t age(const int_fast32_t year) const noexcept {
        return year - birth_year_;
    }

  private:
    //! Pointer to father
    const std::shared_ptr<Individual> father_{nullptr};
    //! Pointer to mother
    const std::shared_ptr<Individual> mother_{nullptr};
    //! Year of birth
    int_fast32_t birth_year_{-4};
};

} // namespace pbf

#endif /* PBT_INDIVIDUAL_HPP_ */
