/*! @file segment.hpp
    @brief Interface of Segment class
*/
#pragma once
#ifndef PBT_SEGMENT_HPP_
#define PBT_SEGMENT_HPP_

#include "individual.hpp"

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////

namespace pbt {

/*! @brief Segment class
*/
class Segment {
  public:
    //! constructor
    Segment(const Individual* i, bool b, bool s=false)
    : individual_(i), is_from_father_(b), is_sampled_(s) {}

    //! setter of #mutations_
    void set_mutations(std::vector<double>&& v) const {
        mutations_ = v;
    }
    //! setter of #ancestral_segment_
    void set_ancestral_segment(Segment* p) const {
        ancestral_segment_ = p;
    }
    //! return pointer to father or mother
    const Individual* ancestor() const {
        if (is_from_father_) {
            return individual_->father_get();
        } else {
            return individual_->mother_get();
        }
    }
    //! convert #mutations_ to boolean vector
    std::vector<bool> binary_genotype(const std::vector<double>& positions) const {
        std::vector<bool> output;
        output.reserve(positions.size());
        std::set<double> genotype;
        accumulate(&genotype);
        for (double x: positions) {
            output.emplace_back(genotype.find(x) != genotype.end());
        }
        return output;
    }

    //! pointer to holder
    const Individual* individual_;
    //! true if #ancestral_segment_ is in father
    const bool is_from_father_;
    //! being sampled or an ancestor of a sample
    const bool is_sampled_;

  private:
    //! accumulate #mutations_ recursively
    void accumulate(std::set<double>* genotype) const {
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

//! Compare Individual pointer
class less_Segment {
  public:
    //! Compare Individual pointer
    bool operator()(const Segment& x, const Segment& y) const {
        return Individual::less{}(x.individual_, y.individual_) || (x.is_from_father_ < y.is_from_father_);
    }
};

} // namespace pbt

#endif /* PBT_SEGMENT_HPP_ */