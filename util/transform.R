library(tidyverse)

min_adult_age = 3L
tsv = read_tsv('output.tsv')

# #######1#########2#########3#########4#########5#########6#########7#########
# POPs

adults = tsv %>%
  dplyr::filter(capture_year - birth_year >= min_adult_age) %>%
  dplyr::transmute(id, capture_year, adult_birth_year= birth_year) %>%
  print()

juveniles = tsv %>%
  dplyr::transmute(id, cohort= birth_year, mother_id, father_id) %>%
  print()

pre_pop = tidyr::crossing(id= tsv$id, adult_id= adults$id) %>%
  dplyr::filter(id > adult_id) %>%
  dplyr::left_join(juveniles, by='id') %>%
  dplyr::left_join(adults, by=c(adult_id='id')) %>%
  dplyr::filter(cohort - adult_birth_year >= min_adult_age) %>%
  dplyr::mutate(is_pop = (adult_id == mother_id | adult_id == father_id)) %>%
  print()

pre_pop %>% dplyr::filter(is_pop)

pop_data = pre_pop %>%
  dplyr::group_by(cohort, capture_year, capture_age= capture_year - adult_birth_year) %>%
  dplyr::summarise(pops= sum(is_pop), comps= n()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(cohort, capture_year, capture_age, fill=list(pops=0L, comps=0L)) %>%
  print()

# #######1#########2#########3#########4#########5#########6#########7#########
# HSPs

tsv_i = tsv %>%
  dplyr::transmute(id, mother_i= mother_id, father_i= father_id, cohort_i= birth_year)
tsv_j = tsv %>%
  dplyr::transmute(id, mother_j= mother_id, father_j= father_id, cohort_j= birth_year)

pre_hsp = tidyr::crossing(id_i= tsv$id, id_j= tsv$id) %>%
  dplyr::filter(id_i < id_j) %>%
  dplyr::left_join(tsv_i, by=c(id_i='id')) %>%
  dplyr::left_join(tsv_j, by=c(id_j='id')) %>%
  dplyr::mutate(is_hsp= (mother_i == mother_j) | (father_i == father_j)) %>%
  print()

pre_hsp %>% dplyr::filter(is_hsp)

hsp_data = pre_hsp %>%
  dplyr::group_by(cohort_i, cohort_j) %>%
  dplyr::summarise(comps= n(), hsps= sum(is_hsp)) %>%
  dplyr::ungroup() %>%
  print()
