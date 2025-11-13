##########################################################
# Chromosome Data Cleaning + Continuous Trait Mapping
# Author: [Your Name]
# Date: [Today’s Date]
##########################################################

# ---- Libraries ----
library(ape)
library(phytools)
library(phangorn)

##########################################################
#                Helper Functions
##########################################################

# Read and clean a tree for a given clade name
GetTree <- function(clade) {
  new_path <- paste0("../data/trees/", clade, ".new")
  nex_path <- paste0("../data/trees/", clade, ".nex")
  
  if (file.exists(new_path)) {
    tree <- read.tree(new_path)
  } else if (file.exists(nex_path)) {
    tree <- read.nexus(nex_path)
  } else {
    stop("No tree file found for clade: ", clade)
  }
  
  duptips <- which(duplicated(tree$tip.label))
  if (length(duptips) > 0)
    tree <- drop.tip(tree, duptips)
  
  return(tree)
}

# Clean haploid column (handles ranges, commas, decimals)
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  
  for (i in 1:nrow(df)) {
    x <- df$haploid[i]
    if (is.numeric(x)) {
      if (!x == round(x))
        df$haploid[i] <- stoch_round(x)
    } else {
      nums <- as.numeric(unlist(strsplit(x, "[^0-9]+")))
      nums <- nums[!is.na(nums)]
      if (length(nums) == 0) next
      if (length(nums) > 1)
        df$haploid[i] <- sample(nums, 1)
      else
        df$haploid[i] <- nums
    }
  }
  df$haploid <- as.numeric(df$haploid)
  return(df)
}

# Ensure bifurcating, ultrametric, scaled tree
AdjustTree <- function(tree, is_plant = FALSE) {
  if (!is.ultrametric(tree)) {
    if (is_plant) {
      tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE)
    } else {
      tree <- chronos(tree)
    }
  }
  tree$edge.length <- tree$edge.length / max(branching.times(tree))
  tree$edge.length[tree$edge.length < 0] <- 1e-9
  if (!is.binary(tree)) tree <- multi2di(tree)
  return(tree)
}

##########################################################
#                Main Script
##########################################################

# ---- User input ----
clade <- "testudines"
data_path <- paste0("../data/chrome/", clade, ".csv")

# ---- Load and clean data ----
dat <- read.csv(data_path, stringsAsFactors = FALSE)
dat <- CleanHaploid(dat)

# ---- Load and clean tree ----
tree <- GetTree(clade)
#tree <- tree[[#]]
#plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "lilaceae")
tree <- AdjustTree(tree, clade %in% plants)

# ---- Match tree and data ----
matched <- intersect(tree$tip.label, dat$species)
tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
dat <- dat[dat$species %in% matched, ]

# ---- Prepare haploid vector ----
hap_vec <- dat$haploid
names(hap_vec) <- dat$species
hap_vec <- hap_vec[tree$tip.label]

# ---- Continuous trait mapping ----
cm <- contMap(tree, hap_vec, plot = FALSE)
cm$tree$tip.label <- sub("^(.)", "\\U\\1", cm$tree$tip.label, perl = TRUE)

# ---- Plot ----
graphics.off()
par(mar = c(1, 1, 1, 1))
par("usr")[3:4]      # current y-limits
par("usr")[1:2]

plotSimmap(
  #type = "fan",
  cm$tree,
  cm$cols,
  lwd = 2,
  fsize = 0.0000008,
  offset = 0.5,
  ylim = c(-7, 126.84)
)

# ---- Add color bar ----
r <- max(nodeHeights(cm$tree))
add.color.bar(
  col   = cm$cols,
  lims  = cm$lims,
  horiz = TRUE,
  fsize = 0.8,
  leg   = r/2,
  x     = 0,
  y     = -5,
  scale = 1.0 * r,
  prompt = FALSE,
  title = "Haploid Chromosome Number"
)

##########################################################
#                     End Script
##########################################################
