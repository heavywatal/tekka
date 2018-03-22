#' convert result to igraph
#' @param .tbl result tibble
#' @return igraph
#' @rdname graph
#' @export
as_igraph = function(.tbl) {
  .tbl %>%
    gather_chromosome() %>%
    dplyr::select(.data$parent_id, .data$id) %>%
    igraph::graph_from_data_frame()
}

#' return vertices with zero degree
#' @inheritParams igraph::degree
#' @return vertices
#' @rdname graph
#' @export
leaf_V = function(graph, mode = "out") {
  igraph::V(graph)[igraph::degree(graph, mode = mode) < 1]
}

#' find kinship below max_degree
#' @inheritParams igraph::ego
#' @return tibble
#' @rdname graph
#' @export
find_kinship = function(graph, nodes, order = 4L) {
  stopifnot(is.character(nodes))
  purrr::map_dfr(seq_len(order), ~{
    tibble::tibble(
      degree = .x,
      from = nodes,
      to = igraph::ego(graph, order = .x, nodes = nodes, mode = "all", mindist = .x) %>%
        purrr::map(., ~.x$name[.x$name %in% nodes])
    )
  }) %>%
    tidyr::unnest() %>%
    dplyr::filter(.data$from < .data$to) %>%
    dplyr::transmute(
      .data$from,
      .data$to,
      .data$degree,
      paths = purrr::pmap(., function(from, to, ...) {
        igraph::all_shortest_paths(graph, from, to, mode = "all", weights = NA)$res %>%
          purrr::map(~.$name)
      })
    ) %>%
    dplyr::arrange(.data$from, .data$to) %>%
    tidyr::unnest()
}
