#' Run C++ simulation
#' @param args command line arguments as a string vector or list of strings
#' @return tibble
#' @rdname blackthunnus
#' @export
blackthunnus = function(args=character(0L)) {
  if (length(args) == 1L) {
    args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
  }
  args = c("blackthunnus", "-q", args)
  message(paste(args, collapse = " "))
  cpp_blackthunnus(args) %>% read_result()
}
