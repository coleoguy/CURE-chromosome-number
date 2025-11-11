### Megan Copeland

## Load libraries
library(ape)           
library(diversitree)   
library(chromePlus)
library(phangorn)

## Functions

# function to read and return a tree object for a given clade name
GetTree <- function(x){
  tree <- NULL
  new_path <- paste0("../data/trees/", x, ".new")  # possible Newick file path
  nex_path <- paste0("../data/trees/", x, ".nex")  # possible Nexus file path
  
  # check which file exists and read accordingly
  if (file.exists(new_path)) {
    tree <- read.tree(new_path)
  } else if (file.exists(nex_path)) {
    tree <- read.nexus(nex_path)
  } else {
    return(NULL)  # skip if no matching tree file found
  }
  
  # clean up any duplicate tip labels
  tree <- GetCleanTree(tree)
  return(tree)
}

# function to remove duplicated tip labels from a phylo or multiPhylo object
GetCleanTree <- function(tree){
  if (inherits(tree, "phylo")) {                   # if single tree
    duptips <- which(duplicated(tree$tip.label))   # find duplicate tips
    if (length(duptips) > 0)
      tree <- drop.tip(tree, duptips)              # remove them
  } else if (inherits(tree, "multiPhylo")) {       # if multiple trees
    tree <- lapply(tree, function(tr) {
      duptips <- which(duplicated(tr$tip.label))
      if (length(duptips) > 0)
        tr <- drop.tip(tr, duptips)
      tr
    })
  }
  return(tree)
}

# function to clean haploid column (to handle ranges, comma seperated values, multiple entries)
CleanHaploid <- function(df) {
  if (!"haploid" %in% names(df)) return(df)
  
  df$haploid <- sapply(df$haploid, function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    x <- as.character(x)
    
    # extract all numeric values (handles ranges, commas, etc.)
    nums <- as.numeric(unlist(strsplit(x, "[^0-9]+")))
    nums <- nums[!is.na(nums)]
    if (length(nums) == 0) return(NA_real_)
    
    # if two numbers and separated by a dash (e.g. "6-8"), sample within range
    if (length(nums) == 2 && grepl("-", x)) 
      return(sample(seq(min(nums), max(nums)), 1))
    
    # if multiple distinct values, use mode if unique, else random
    if (length(nums) > 1) {
      tab <- table(nums)
      modeval <- as.numeric(names(tab)[tab == max(tab)])
      return(if (length(modeval) == 1) modeval else sample(nums, 1))
    }
    
    # single value
    return(nums)
  })
  
  # ensure numeric (sapply sometimes returns character)
  df <- df[!is.na(df$haploid),]
  df$haploid <- as.numeric(df$haploid)
  
  return(df)
}

## Main Loop

# get all plant clades that are using the angiosperm phylo
plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "lilaceae")

# list all chromosome data files
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = TRUE)

# initialize an empty list to store MCMC results
results <- list()

# iterate through each chrome data file
for (i in 1:length(file_list)) {
  f <- file_list[[i]]
  print(f)
  clade <- tools::file_path_sans_ext(basename(f))  # extract clade name from filename
  
  # read chromosome data for this clade
  dat <- read.csv(f)
  dat <- CleanHaploid(dat)
  
  # read the corresponding tree
  tree <- GetTree(clade)
  
  # match tree and data by species names
  matched <- intersect(tree$tip.label, dat$species)

  # drop non-matching taxa from tree and data
  tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
  
  dat <- dat[dat$species %in% matched, ]
  dat <- dat[, c("species", "haploid")]
  
  # convert data to matrix format required by diversitree
  mat <- datatoMatrix(x = dat, buffer = 1, hyper = FALSE)
  

  # ensure the tree is ultrametric
  if (!is.ultrametric(tree)) {
    if (clade %in% plants) {
      tree <- nnls.tree(cophenetic(tree), tree, rooted = T)
    } else {
      tree <- chronos(tree)
    }
  }
  
  # scale tree to unit length
  tree$edge.length <- tree$edge.length / max(branching.times(tree))
  
  # replace any small negative branch lengths with tiny positive value
  tree$edge.length[tree$edge.length < 0] <- 1e-9
  
  # ensure tree is fully bifurcating 
  if (!is.binary(tree)) {
    tree <- multi2di(tree)
  }
  
  lik <- make.mkn(tree = tree, states = mat, k=ncol(mat), 
                  strict=F, control=list(method="ode",root=ROOT.OBS))
  argnames(lik)
  
  #  costrain our model to match the dynamics of chromosome evolution
  conlik <- constrainMkn(data = mat, lik = lik, hyper=F,
                         polyploidy = T, verbose=F)
  
  argnames(conlik)
  
  # run mcmc
  res <- mcmc(lik=conlik, x.init=runif(length(argnames(conlik))),
              prior=make.prior.exponential(2), nsteps=1, w=1)
  
  # store results 
  if (!is.null(res)) {
    results[[clade]] <- res
  }
}
