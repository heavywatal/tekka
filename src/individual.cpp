/*! @file individual.cpp
    @brief Implementation of Individual class
*/
#include "individual.hpp"

#include <wtl/debug.hpp>

#include <type_traits>

namespace pbf {

static_assert(std::is_nothrow_default_constructible_v<Individual>);
static_assert(!std::is_copy_constructible_v<Individual>);
static_assert(!std::is_move_constructible_v<Individual>);
static_assert(std::is_nothrow_destructible_v<Individual>);

void Individual::trace_back(std::ostream& ost, std::unordered_map<const Individual*, uint_fast32_t>& ids,
                            uint_fast32_t loc, int_fast32_t year) const {
    if (!ids.emplace(this, static_cast<uint_fast32_t>(ids.size())).second && (year == 0)) return;
    if (father_) father_->trace_back(ost, ids, loc, 0);
    if (mother_) mother_->trace_back(ost, ids, loc, 0);
    write(ost, ids);
    if (year > 0) {
        ost << "\t" << loc << "\t" << year << "\n";
    } else {
        ost << "\t\t\n";
    }
}

std::ostream& Individual::write_trace_back_header(std::ostream& ost) {
   return write_names(ost) << "\tlocation\tcapture_year\n";
}

std::ostream& Individual::write_names(std::ostream& ost) {
    return ost << "id\t"
               << "father_id\t"
               << "mother_id\t"
               << "birth_year";
}

std::ostream& Individual::write(std::ostream& ost) const {
    return ost << this << "\t"
               << father_ << "\t"
               << mother_ << "\t"
               << birth_year_;
}

std::ostream& Individual::write(std::ostream& ost, const std::unordered_map<const Individual*, uint_fast32_t>& ids) const {
    return ost << ids.at(this) << "\t"
               << ids.at(father_.get()) << "\t"
               << ids.at(mother_.get()) << "\t"
               << birth_year_;
}

//! Shortcut of Individual::write
std::ostream& operator<<(std::ostream& ost, const Individual& x) {
    return x.write(ost);
}

} // namespace pbf
