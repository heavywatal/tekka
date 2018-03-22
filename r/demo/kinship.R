library(tekkamaki)

result = blackthunnus('-y100 -l5') %>% print()
kinship = find_kinship(result, 8L) %>% print()
