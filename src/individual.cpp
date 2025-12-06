#include "individual.hpp"

#include <fmt/format.h>
#include <wtl/debug.hpp>

#include <iterator>

namespace pbf {

using OutputIt = std::back_insert_iterator<fmt::memory_buffer>;

template <>
OutputIt Individual::format_to(OutputIt out) const {
    return fmt::format_to(
      out, "{}\t{}\t{}\t{}",
      fmt::ptr(this), fmt::ptr(father_.get()), fmt::ptr(mother_.get()), birth_year_
    );
}

template <>
OutputIt Individual::format_to(OutputIt out,
  const std::unordered_map<const Individual*, int_fast32_t>& ids) const {
    return fmt::format_to(
      out, "{}\t{}\t{}\t{}",
      ids.at(this), ids.at(father_.get()), ids.at(mother_.get()), birth_year_
    );
}

template <>
void Individual::trace_back(OutputIt out,
  std::unordered_map<const Individual*, int_fast32_t>& ids,
  int_fast32_t loc, int_fast32_t year) const {
    if (!ids.emplace(this, static_cast<int_fast32_t>(ids.size())).second && (year == 0)) return;
    if (father_) father_->trace_back(out, ids, loc, 0);
    if (mother_) mother_->trace_back(out, ids, loc, 0);
    format_to(out, ids);
    if (year > 0) {
        fmt::format_to(out, "\t{}\t{}\n", loc, year);
    } else {
        fmt::format_to(out, "\t\t\n");
    }
}

std::string Individual::str() const {
    fmt::memory_buffer buffer;
    format_to(std::back_inserter(buffer));
    return fmt::to_string(buffer);
}

std::string Individual::trace_back_header() {
    return fmt::format("{}\tlocation\tcapture_year", header());
}

std::string_view Individual::header() noexcept {
    return "id\tfather_id\tmother_id\tbirth_year";
}

} // namespace pbf
