library(tekkamaki)

result = blackthunnus('-y100 -l5') %>% print()
sampled = dplyr::filter(result, !is.na(capture_year)) %>% print()
graph = as_igraph(result) %>% print()
kinship = find_kinship(graph, sampled$id, 8L) %>% print()
