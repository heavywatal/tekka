pairwise_half_sibling = function(.tbl) {
  .tbl = dplyr::filter(.tbl, !is.na(.data$capture_year))
  tsv_i = .tbl %>%
    dplyr::transmute(.data$id, mother_i = .data$mother_id, father_i = .data$father_id, cohort_i = .data$birth_year)
  tsv_j = .tbl %>%
    dplyr::transmute(.data$id, mother_j = .data$mother_id, father_j = .data$father_id, cohort_j = .data$birth_year)
  tidyr::crossing(id_i = .tbl$id, id_j = .tbl$id) %>%
    dplyr::filter(.data$id_i < .data$id_j) %>%
    dplyr::left_join(tsv_i, by = c(id_i = "id")) %>%
    dplyr::left_join(tsv_j, by = c(id_j = "id")) %>%
    dplyr::mutate(is_hsp = (.data$mother_i == .data$mother_j) | (.data$father_i == .data$father_j))
}

summarize_hsp = function(.tbl) {
  .tbl %>%
    dplyr::group_by(.data$cohort_i, .data$cohort_j) %>%
    dplyr::summarise(comps = dplyr::n(), hsps = sum(.data$is_hsp)) %>%
    dplyr::ungroup()
}

#' Convert a result data frame to HSP format
#' @param .tbl A result data frame
#' @return tibble with additional class name "hsp"
#' @rdname hsp
#' @export
as_hsp = function(.tbl) {
  .tbl = .tbl %>%
    pairwise_half_sibling() %>%
    summarize_hsp()
  class(.tbl) = c(class(.tbl), "hsp")
  .tbl
}

#' Write a HSP data frame to a file
#' @param x An outcome of as_hsp()
#' @param path A file name or connection to write to
#' @rdname hsp
#' @export
write_hsp = function(x, path="hsp.txt") {
  stopifnot(inherits(x, "hsp"))
  lines = c(
    "# HSP 1-false-negative ratio",
    "1.0",
    "# number of HSP data points",
    nrow(x),
    "# HSP data"
  )
  readr::write_lines(lines, path)
  readr::write_tsv(x, path, na = "", append = TRUE, col_names = FALSE)
}
