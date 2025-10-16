# Heath Blackmon
# coleoguy@gmail.com
# This file runs and check that it can find a tree for each dataset 
# in the chrome folder. For each clade it records the number of chromosome
# records, the number of tips on tree, the type of tree 
# (ultrametric, 1 or posterior) and the number of overlapping species.

library(dplyr)

# get a list of the chromosome files
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = TRUE)

# Initialize list to store data
dat <- list()
## Reading the data
for (file in file_list) {
  df <- read.csv(file)
  if (!"haploid" %in% names(df)) {
    cat("Missing 'haploid' column in:", basename(file), "\n")
    next
  }
  file_name <- tools::file_path_sans_ext(basename(file))
  temp_df <- data.frame(clade = file_name, species = df$species, haploid = df$haploid)
  dat[[length(dat) + 1]] <- temp_df
}
all.dat <- do.call(rbind, dat)
col.dat <- nrow(all.dat)
rm(list=ls()[-c(1,2,3, 6)])

# Fix the haploid column by parsing multiple formats (commas, dashes, ranges)
all.dat <- all.dat %>%
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
use.dat <- nrow(all.dat)


# now lets start going through each clade and evaluating
# the available data. First we set up three vectors to track
tree.tips <- chrome.tips1 <- clades <-
  chrome.tips2 <- overlap.tips <- c()
res <- list()

for(i in 1:length(file_list)){
  tree.tips[i] <- NA
  cur.clade <- strsplit(basename(file_list[i]), 
                        split=".", fixed=T)[[1]][1]
  clades[i] <- cur.clade
  chrome.tips1[i] <- nrow(dat[[i]])
  chrome.tips2[i] <- nrow(all.dat[all.dat$clade==cur.clade,])
  fnew <- paste0("../data/tree/", cur.clade, ".new")
  if(file.exists(fnew)){
    res[[i]] <- ape::read.tree(fnew)
    tree.tips[i] <- length(res[[i]]$tip.label)
  }else{
    fnew <- paste0("../data/tree/", cur.clade, ".nex")
    if(file.exists(fnew)){
      res[[i]] <- ape::read.nexus(fnew)
      tree.tips[i] <- length(res[[i]]$tip.label)
    }
  }
  if(is.na(tree.tips[i])){
    tree.tips[i] <- "tree missing"
  }
}
trees <- res
restab <- data.frame(clades, chrome.tips1, chrome.tips2, tree.tips)
colnames(restab) <- c("clade","collected.data","filtered.data","tree.size")

rm(list=ls()[-c(7,13,15)])



sum(dat[[9]]$species %in% trees[[9]]$tip.label)


READY TO WRITE THE PART THAT DOES THE TREE MATCHING





cat(paste("Total Records Collected:", 
          col.dat, "\nTotal Usable Records:", 
          use.dat))




