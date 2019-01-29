/*! @file segment.hpp
    @brief Interface of Segment class
*/
#pragma once
#ifndef PBT_SEGMENT_HPP_
#define PBT_SEGMENT_HPP_

#include "individual.hpp"

#include <set>
#include <type_traits>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbf {

/*! @brief Segment class
*/
class Segment {
  public:
    //! @name Constructors
    //@{
    //! Constructor
    Segment(const Individual* i, bool b, bool s=false)
    : individual_(i), is_from_father_(b), is_sampled_(s) {}
    Segment() = delete;
    Segment(const Segment&) = delete;
    //! Move constructor
    Segment(Segment&&) = default;
    //@}

    //! @name Setter functions
    //@{
    //! Set #mutations_
    void set_mutations(std::vector<double>&& v) const noexcept {
        mutations_ = std::move(v);
    }
    //! Set #ancestral_segment_
    void set_ancestral_segment(Segment* p) const noexcept {
        ancestral_segment_ = p;
    }
    //@}
    //! return pointer to father or mother
    const Individual* ancestor() const noexcept {
        if (is_from_father_) {
            return individual_->father_get();
        } else {
            return individual_->mother_get();
        }
    }
    //! convert #mutations_ to boolean vector
    std::vector<bool> binary_genotype(const std::vector<double>& positions) const noexcept {
        std::vector<bool> output;
        output.reserve(positions.size());
        std::set<double> genotype;
        accumulate(&genotype);
        for (double x: positions) {
            output.emplace_back(genotype.find(x) != genotype.end());
        }
        return output;
    }
    //! operator for std::set
    bool operator<(const Segment& other) const noexcept {
        return (&individual_ < &other.individual_) ||
               ((&individual_ == &other.individual_) &&
                (is_from_father_ < other.is_from_father_));
    }

    //! pointer to holder
    const Individual* individual_;
    //! true if #ancestral_segment_ is in father
    const bool is_from_father_;
    //! being sampled or an ancestor of a sample
    const bool is_sampled_;

  private:
    //! accumulate #mutations_ recursively
    void accumulate(std::set<double>* genotype) const noexcept {
        genotype->insert(mutations_.begin(), mutations_.end());
        if (ancestral_segment_) {
            ancestral_segment_->accumulate(genotype);
        }
    }

    //! positions of mutated sites
    mutable std::vector<double> mutations_;
    //! pointer to ancestral segment
    mutable Segment* ancestral_segment_ = nullptr;
};

static_assert(!std::is_default_constructible<Segment>{}, "");
static_assert(!std::is_copy_constructible<Segment>{}, "");
static_assert(std::is_nothrow_move_constructible<Segment>{}, "");

} // namespace pbf

#endif /* PBT_SEGMENT_HPP_ */
