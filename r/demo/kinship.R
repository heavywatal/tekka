library(tekkamaki)

result = blackthunnus('-y100 -l5') %>% print()
sample_family = result$sample_family[[1L]]
kinship = find_kinship(sample_family, 4L) %>% print()

location_map = sample_family %>% dplyr::select(id, location)
kinship %>%
  dplyr::left_join(location_map, by = c(from="id")) %>%
  dplyr::rename(from_location = location) %>%
  dplyr::left_join(location_map, by = c(to="id")) %>%
  dplyr::rename(to_location = location) %>%
  dplyr::filter(to_location == 1L, from_location == 1L)
