pairwise_parent_offspring = function(.tbl, min_adult_age) {
  .tbl = dplyr::filter(.tbl, !is.na(.data$capture_year))
  adults = .tbl %>%
    dplyr::filter(.data$capture_year - .data$birth_year >= min_adult_age) %>%
    dplyr::transmute(.data$id, .data$capture_year, adult_birth_year = .data$birth_year)
  juveniles = .tbl %>%
    dplyr::transmute(.data$id, cohort = .data$birth_year, .data$mother_id, .data$father_id)
  tidyr::crossing(id = .tbl$id, adult_id = adults$id) %>%
    dplyr::filter(.data$id > .data$adult_id) %>%
    dplyr::left_join(juveniles, by = "id") %>%
    dplyr::left_join(adults, by = c(adult_id = "id")) %>%
    dplyr::filter(.data$cohort - .data$adult_birth_year >= min_adult_age) %>%
    dplyr::mutate(is_pop = (.data$adult_id == .data$mother_id | .data$adult_id == .data$father_id))
}

summarize_pop = function(.tbl) {
  .tbl %>%
    dplyr::mutate(capture_age = .data$capture_year - .data$adult_birth_year) %>%
    dplyr::group_by(.data$cohort, .data$capture_year, .data$capture_age) %>%
    dplyr::summarise(pops = sum(.data$is_pop), comps = n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete_(c(~cohort, ~capture_year, ~capture_age), fill = list(pops = 0L, comps = 0L))
}

#' Convert a result data frame to POP format
#' @param .tbl result data frame
#' @param min_adult_age integer
#' @return tibble
#' @rdname pop
#' @export
as_pop = function(.tbl, min_adult_age=5L) {
  .tbl = .tbl %>%
    pairwise_parent_offspring(min_adult_age) %>%
    summarize_pop()
  class(.tbl) = c(class(.tbl), "pop")
  .tbl
}

#' Write a POP data frame to a file
#' @inheritParams readr::write_tsv
#' @rdname pop
#' @export
write_pop = function(x, path="pop.txt") {
  stopifnot(inherits(x, "pop"))
  readr::write_file("# ckdat:\n", path)
  readr::write_tsv(x, path, na = "", append = TRUE, col_names = FALSE)
}
