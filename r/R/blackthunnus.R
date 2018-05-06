#' Run C++ simulation
#' @param args command line arguments as a string vector
#' @return tibble
#' @rdname blackthunnus
#' @export
blackthunnus = function(args=character(0L)) {
  if (length(args) == 1L) {
    args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
  }
  args = c("blackthunnus", "-q", args)
  message(paste(args, collapse = " "))
  .out = cpp_blackthunnus(args)
  if (length(.out) > 0L) read_result(.out)
}
