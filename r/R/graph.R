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

#' find kinship below given order
#' @param order integer
#' @return tibble
#' @rdname graph
#' @export
find_kinship = function(.tbl, order = 4L) {
  graph = as_igraph(.tbl)
  sampled_nodes = dplyr::filter(.tbl, !is.na(.data$capture_year))$id
  birth_year = stats::setNames(.tbl$birth_year, .tbl$id)
  find_kinship_impl(graph, sampled_nodes, order = order) %>%
    find_shortest_paths(graph) %>%
    dplyr::mutate(direction = purrr::map(.data$path, ~as.integer(diff(birth_year[.x]) < 0L))) %>%
    dplyr::mutate(backward = purrr::map_int(.data$direction, sum))
}

# find pairs
find_kinship_impl = function(graph, nodes, order) {
  stopifnot(is.character(nodes))
  tibble::tibble(
    from = nodes,
    to = igraph::ego(graph, order = order, nodes = nodes, mode = "all") %>%
      purrr::map(., ~.x$name[.x$name %in% nodes])
  ) %>%
    tidyr::unnest() %>%
    dplyr::filter(.data$from < .data$to) %>%
    dplyr::arrange(.data$from, .data$to)
}

# add columns: path and degree
find_shortest_paths = function(kinship, graph) {
  kinship %>%
    dplyr::mutate(
      path = purrr::pmap(., function(from, to, ...) {
        igraph::all_shortest_paths(graph, from, to, mode = "all", weights = NA)$res %>%
          purrr::map(~.$name)
      })
    ) %>%
    tidyr::unnest() %>%
    dplyr::filter(!purrr::map_lgl(.data$path, ~"0x0" %in% .x)) %>%
    dplyr::mutate(degree = lengths(.data$path) - 1L)
}
