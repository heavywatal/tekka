kinlabels = c(
  `1_1` = "PO",
  `0_1` = "PO",
  `10_1` = "HS",
  `10_2` = "FS",
  `11_1` = "GG",
  `00_1` = "GG",
  `110_1` = "HUN",
  `100_1` = "HUN",
  `110_2` = "UN",
  `100_2` = "UN",
  `1100_1` = "HC",
  `1100_2` = "C"
)

devtools::use_data(kinlabels, overwrite = TRUE, internal = TRUE)
