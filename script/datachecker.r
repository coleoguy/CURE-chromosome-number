# Heath Blackmon
# coleoguy@gmail.com
# This file runs and check that it can find a tree for each dataset 
# in the chrome folder. For each clade it records the number of chromosome
# records, the number of tips on tree, the type of tree 
# (ultrametric, 1 or posterior) and the number of overlapping species.

# get a list of the chromosome files
file_list <- list.files(path = "../data/chrome/", pattern = "\\.csv$", full.names = TRUE)

# Initialize list to store data
all_data <- list()
## Reading the data
for (file in file_list) {
  df <- read.csv(file)
  if (!"haploid" %in% names(df)) {
    cat("Missing 'haploid' column in:", basename(file), "\n")
    next
  }
  file_name <- tools::file_path_sans_ext(basename(file))
  temp_df <- data.frame(clade = file_name, haploid = df$haploid)
  all_data[[length(all_data) + 1]] <- temp_df
}
dat <- do.call(rbind, all_data)
rm(list=ls()[-2])

# Fix the haploid column by parsing multiple formats (commas, dashes, ranges)
dat_clean <- dat %>%
  mutate(
    haploid = str_trim(haploid),
    haploid = map_dbl(haploid, ~ {
      vals <- str_split(.x, ",|–|-")[[1]] |> str_trim()
      nums <- suppressWarnings(as.numeric(vals))
      nums <- nums[!is.na(nums) & nums > 0]
      if (length(nums) == 0) return(NA_real_)
      if (length(nums) == 1) return(nums)
      if (length(nums) == 2 && diff(nums) > 0) return(runif(1, nums[1], nums[2]))
      return(sample(nums, 1))  # fallback: sample any valid value
    })
  ) %>%
  filter(is.finite(haploid), haploid > 0)
dat <- dat_clean
rm(dat_clean)



