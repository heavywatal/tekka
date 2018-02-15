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
    Segment(const Individual* i, bool b, bool s=false)
    : individual_(i), is_from_father_(b), is_sampled_(s) {}

    struct less {
        bool operator()(const Segment& x, const Segment& y) const {
            return Individual::less{}(x.individual_, y.individual_) || (x.is_from_father_ < y.is_from_father_);
        }
    };

    void set_mutations(std::vector<double>&& v) const {
        mutations_ = v;
    }
    void set_ancestral_segment(Segment* p) const {
        ancestral_segment_ = p;
    }
    const Individual* ancestor() const {
        if (is_from_father_) {
            return individual_->father_get();
        } else {
            return individual_->mother_get();
        }
    }
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

    bool is_first_gen() const {return individual_->is_first_gen();}
    std::vector<double>::const_iterator begin() const {return mutations_.begin();}
    std::vector<double>::const_iterator end() const {return mutations_.end();}
    bool is_sampled() const {return is_sampled_;}

  private:
    void accumulate(std::set<double>* genotype) const {
        genotype->insert(begin(), end());
        if (ancestral_segment_) {
            ancestral_segment_->accumulate(genotype);
        }
    }

    const Individual* individual_;
    const bool is_from_father_;
    const bool is_sampled_;
    mutable std::vector<double> mutations_;
    mutable Segment* ancestral_segment_ = nullptr;
};

} // namespace pbt

#endif /* PBT_SEGMENT_HPP_ */
