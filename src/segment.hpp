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
    Segment(const Individual* i, bool b)
    : individual(i), is_from_father(b) {}

    struct less {
        bool operator()(const Segment& x, const Segment& y) const {
            return Individual::less{}(x.individual, y.individual) || (x.is_from_father < y.is_from_father);
        }
    };

    void set_mutations(std::vector<double>&& v) const {
        mutations_ = v;
    }
    void set_ancestral_segment(Segment* p) const {
        ancestral_segment_ = p;
    }
    const Individual* ancestor() const {
        if (is_from_father) {
            return individual->father_get();
        } else {
            return individual->mother_get();
        }
    }
    void accumulate(std::set<double>* genotype) const {
        genotype->insert(begin(), end());
        if (ancestral_segment_) {
            ancestral_segment_->accumulate(genotype);
        }
    }

    bool is_first_gen() const {return individual->is_first_gen();}
    std::vector<double>::const_iterator begin() const {return mutations_.begin();}
    std::vector<double>::const_iterator end() const {return mutations_.end();}

    const Individual* individual;
    const bool is_from_father;
  private:
    mutable std::vector<double> mutations_;
    mutable Segment* ancestral_segment_ = nullptr;
};

} // namespace pbt

#endif /* PBT_SEGMENT_HPP_ */
