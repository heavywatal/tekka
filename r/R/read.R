#' Run C++ simulation
#' @inheritParams readr::read_tsv
#' @return tibble
#' @rdname read
#' @export
read_result = function(file) {
  readr::read_tsv(file, col_types = readr::cols(capture_year = "i"))
}
