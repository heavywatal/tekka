library(tidyverse)

min_adult_age = 5L
tsv = read_tsv('sample_family.tsv.gz', col_types = cols(capture_year='i')) %>%
  dplyr::filter(!is.na(capture_year)) %>%
  print()

# #######1#########2#########3#########4#########5#########6#########7#########
# POPs

pairwise_parent_offspring = function(.tbl) {
    adults = .tbl %>%
      dplyr::filter(capture_year - birth_year >= min_adult_age) %>%
      dplyr::transmute(id, capture_year, adult_birth_year= birth_year)
    juveniles = .tbl %>%
      dplyr::transmute(id, cohort= birth_year, mother_id, father_id)
    tidyr::crossing(id= .tbl$id, adult_id= adults$id) %>%
    dplyr::filter(id > adult_id) %>%
    dplyr::left_join(juveniles, by='id') %>%
    dplyr::left_join(adults, by=c(adult_id='id')) %>%
    dplyr::filter(cohort - adult_birth_year >= min_adult_age) %>%
    dplyr::mutate(is_pop = (adult_id == mother_id | adult_id == father_id))
}

summarize_pop = function(.tbl) {
    .tbl %>%
    dplyr::group_by(cohort, capture_year, capture_age= capture_year - adult_birth_year) %>%
    dplyr::summarise(pops= sum(is_pop), comps= n()) %>%
    dplyr::ungroup() %>%
    tidyr::complete(cohort, capture_year, capture_age, fill=list(pops=0L, comps=0L))
}

write_pop = function(.tbl, path='pop.txt') {
    readr::write_file('# ckdat:\n', path)
    readr::write_tsv(.tbl, path, na='', append=TRUE, col_names=FALSE)
}

pop_data = tsv %>%
    pairwise_parent_offspring() %>% print() %>%
    summarize_pop() %>% print()

write_pop(pop_data)

# #######1#########2#########3#########4#########5#########6#########7#########
# HSPs

pairwise_half_sibling = function(.tbl) {
    tsv_i = .tbl %>%
      dplyr::transmute(id, mother_i= mother_id, father_i= father_id, cohort_i= birth_year)
    tsv_j = .tbl %>%
      dplyr::transmute(id, mother_j= mother_id, father_j= father_id, cohort_j= birth_year)
    tidyr::crossing(id_i= .tbl$id, id_j= .tbl$id) %>%
    dplyr::filter(id_i < id_j) %>%
    dplyr::left_join(tsv_i, by=c(id_i='id')) %>%
    dplyr::left_join(tsv_j, by=c(id_j='id')) %>%
    dplyr::mutate(is_hsp= (mother_i == mother_j) | (father_i == father_j))
}

summarize_hsp = function(.tbl) {
    .tbl %>%
    dplyr::group_by(cohort_i, cohort_j) %>%
    dplyr::summarise(comps= n(), hsps= sum(is_hsp)) %>%
    dplyr::ungroup()
}

write_hsp = function(.tbl, path='hsp.txt') {
    lines = c('# HSP 1-false-negative ratio',
      '1.0',
      '# number of HSP data points',
      nrow(.tbl),
      '# HSP data')
    readr::write_lines(lines, path)
    readr::write_tsv(.tbl, path, na='', append=TRUE, col_names=FALSE)
}

hsp_data = tsv %>%
  pairwise_half_sibling() %>% print() %>%
  summarize_hsp() %>% print()

write_hsp(hsp_data)
