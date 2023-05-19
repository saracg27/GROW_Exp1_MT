######################################################################

################### Functions to into the environment ###############
# This are functions that I created to keep the R code in the #######
# project cleaner. They can be copy pasted and modified as needed. ##


######################################################################

#map through all the elements of times vector, subset and write table

subExpMeta <- function(meta, times, output_dir) {
   purrr::map(times, ~ {
      to_keep <- meta$Time == .x
      meta_subset <- meta[to_keep, ]
      file_name <- paste0("meta", .x, ".csv")
      file_path <- file.path(output_dir, file_name)
      write.csv(meta_subset, file = file_path, row.names = FALSE)
   })
}
# meta: object to subset
# times: vector to use for mapping and subsetting
# output_dir: directory to store the subsetting output