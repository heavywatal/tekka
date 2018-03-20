#' Make snp table from a result
#' @param .tbl result tibble
#' @param lambda mutation rate per haploid per generation
#' @param rho recombination rate per generation
#' @rdname snp
#' @export
make_snp = function(.tbl, lambda, rho=0) {
  if (!missing(rho)) {
    stop('Recombination has not been implemented yet')
  }
  .tbl %>%
    gather_chromosome() %>%
    sprinkle_mutations(lambda=lambda) %>%
    accumulate_mutations() %>%
    complete_genotype()
}

# row chromosome
gather_chromosome = function(.tbl) {
  .tbl %>%
    dplyr::transmute(.data$id, .data$father_id, .data$mother_id, is_sampled = !is.na(.data$capture_year)) %>%
    tidyr::gather('chr', 'parent_id', -'id', -'is_sampled') %>%
    dplyr::mutate(chr = c(father_id=1L, mother_id=2L)[.data$chr]) %>%
    dplyr::arrange(.data$id, .data$chr)
}

# choose random mutation positions [0, 1)
sprinkle_mutations = function(.tbl, lambda) {
  stopifnot(all(c('chr', 'parent_id') %in% names(.tbl)))
  dplyr::mutate(.tbl,
    num_mutation = stats::rpois(n(), lambda),
    pos = purrr::map(.data$num_mutation, stats::runif, min=0.0, max=1.0),
    num_mutation = NULL
  )
}

.accumulate_mutations_f = function(.id, .chr, .tbl) {
  .row = dplyr::filter(.tbl, .data$id == .id, .data$chr == .chr)
  .pos = .row$pos
  .id = .row$parent_id
  while (.id != "0x0") {
    .chr = sample.int(2L, 1L)
    .row = dplyr::filter(.tbl, .data$id == .id, .data$chr == .chr)
    .pos = c(.pos, .row$pos)
    .id = .row$parent_id
  }
  purrr::flatten_dbl(.pos)
}

# trace back ancestors to collect mutation positions
accumulate_mutations = function(.tbl) {
  .tbl %>%
    dplyr::filter(.data$is_sampled) %>%
    dplyr::transmute(
      .data$id,
      .data$chr,
      pos = purrr::map2(.data$id, .data$chr, .accumulate_mutations_f, .tbl=.tbl)
    )
}

# complete genotypes for all the positions
complete_genotype = function(.tbl) {
  .tbl %>%
    dplyr::mutate(genotype=1L) %>%
    tidyr::unnest() %>%
    tidyr::complete_(c('id', 'chr', 'pos'), fill=list(genotype=0L))
}

# spread data frame to wide format
spread_genotype = function(.tbl) {
  .tbl %>%
    dplyr::mutate(pos = .data$pos %>% as.factor() %>% as.integer() %>% {sprintf('snp%04d', .)}) %>%
    tidyr::spread('pos', 'genotype', fill=0L)
}
