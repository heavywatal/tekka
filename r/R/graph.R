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
    dplyr::mutate(path = purrr::map(.data$path, ~as.integer(diff(birth_year[.x]) < 0L))) %>%
    dplyr::filter(purrr::map_lgl(.data$path, ~length(rle(.x)$lengths) < 3L && !identical(.x, c(0L, 1L)))) %>%
    label_kinship()
}

# find pairs
find_kinship_impl = function(graph, nodes, order) {
  stopifnot(is.character(nodes))
  .ego = igraph::ego(graph, order = order, nodes = nodes, mode = "all")
  tibble::tibble(
    from = nodes,
    to = purrr::map(.ego, ~.x$name[.x$name %in% nodes])
  ) %>%
    tidyr::unnest() %>%
    dplyr::filter(.data$from < .data$to) %>%
    dplyr::arrange(.data$from, .data$to)
}

# add columns: path and degree
find_shortest_paths = function(kinship, graph) {
  .path = purrr::pmap(kinship, function(from, to, ...) {
    igraph::all_shortest_paths(graph, from, to, mode = "all", weights = NA)$res %>%
      purrr::map(igraph::as_ids)
  })
  kinship %>%
    dplyr::mutate(path = .path) %>%
    tidyr::unnest() %>%
    dplyr::filter(!purrr::map_lgl(.data$path, ~"0x0" %in% .x))
}

label_kinship = function(kinship) {
  kinship %>%
    dplyr::mutate(path = purrr::map_chr(.data$path, paste0, collapse="")) %>%
    dplyr::count(.data$from, .data$to, .data$path) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      degree = nchar(.data$path),
      path = paste(.data$path, .data$n, sep="_"),
      n = NULL,
      label = kinlabels[.data$path]
    )
}
