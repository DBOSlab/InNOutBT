## Code to prepare the package's associated data

library(treeio)
library(usethis)
library(utils)

vatsData <- utils::read.csv("inst/extdata/vataireoids_1610_25May2022_BayesTraits.csv")
vatsTree <- treeio::read.beast("inst/extdata/vatcomb30_24May22_Yule.con.tree")

# Adding geo/climate data
usethis::use_data(vatsData, overwrite = TRUE)

# Adding tree data
usethis::use_data(vatsTree, overwrite = TRUE)
