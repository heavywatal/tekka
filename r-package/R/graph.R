#' convert result to igraph
#' @param .tbl result tibble
#' @return igraph
#' @rdname graph
#' @export
as_igraph = function(.tbl) {
  .tbl %>%
    dplyr::filter(.data$father_id > 0L) %>%
    gather_chromosome() %>%
    dplyr::select(.data$parent_id, .data$id) %>%
    igraph::graph_from_data_frame()
}
