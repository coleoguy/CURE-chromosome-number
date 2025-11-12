# This script is meant to run the final analysis of
# your data set. The final result will be a CSV file
# that saves to your computer.

# Your job today is to read this document line by
# line and do what it says and document any problems
# come up in this process.

# These are the packages that we will use today
# if you do not have these installed you will have
# to install prior to runnning the rest of the code.
library(ape)           
library(diversitree)   
library(chromePlus)
library(phangorn)

# We have written a number of functions that the
# analysis script will use. This is the preffered
# method of doing big analyses like this. You 
# will often write functions that do big chuncks
# of the analysis that you can test individually
# and then combine together to do something even
# more complex.

# in R when you see a long series of # it can
# allow you to collapse or expand sections of code
# in the lines below look for a small downward 
# pointing triangle. If you press that it will
# collapse the code till it reaches the next
# occurence of the same number of # symbols.
# Try pressing the button below on line 33.

###### Functions ######
# function to read and return a tree object for a given clade name
GetTree <- function(x){
  tree <- NULL
  new_path <- paste0(x, ".new")  # possible Newick file path
  nex_path <- paste0(x, ".nex")  # possible Nexus file path
  
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
  if(class(tree) == "phylo"){                   # if single tree
    duptips <- which(duplicated(tree$tip.label))   # find duplicate tips
    if (length(duptips) > 0)
      tree <- drop.tip(tree, duptips)              # remove them
  }
  if(class(tree) == "multiPhylo"){       # if multiple trees
    for(i in 1:length(tree)){
      tr <- tree[[i]]
      duptips <- which(duplicated(tr$tip.label))
      if(length(duptips) > 0)
        tree[[i]] <- drop.tip(tr, duptips)
    }
  }
  return(tree)
}
# function to clean haploid column (to handle ranges, comma seperated values, multiple entries)
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  for(i in 1:nrow(df)){
    x <- df$haploid[i]
    if(is.numeric(x)){
      if(!x == round(x)){
        df$haploid[i] <- stoch_round(x)
      }
    }else{
      nums <- as.numeric(strsplit(x, split="-", fixed=T)[[1]])
      if(length(nums) > 2){
        stop(paste("Looks like you have too many dashes! Check row:",i))
      }
      if(length(nums) == 2){
        df$haploid[i] <- sample(nums[1]:nums[2], 1)
      }
      nums <- as.numeric(strsplit(x, split=",", fixed=T)[[1]])
      if(length(nums)>1){
        df$haploid[i] <- sample(nums, 1)
      }
    }
  }
  df$haploid <- stoch_round(as.numeric(df$haploid))
  return(df)
}
# function to make sure trees are bifurcating
# and ultrametric.
AdjustTree <- function(tree, plants){
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
  return(tree)
}


###### End Functions ######

# This allows us to hide lines 33 to 117 of
# the code. Highlight lines 32 to 118 and hit
# run. When you do this you should see 4 
# functions appear in your environment.

# Main Analysis

####################################
# In the lines below there are     #
# two spots where you will enter.  #
# information for your clade       #
####################################

clade <- "YOUR CLADE NAME HERE"
# read chromosome data 
dat <- read.csv("YOUR FILE NAME HERE")
dat <- CleanHaploid(dat)

# read the corresponding tree
tree <- GetTree(clade)

# get all plant clades that are using the angiosperm phylo
plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "lilaceae")

# match tree and data by species names
if(class(tree) == "phylo"){
  matched <- intersect(tree$tip.label, dat$species)
  tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
}
if(class(tree) == "multiPhylo"){
  result <- list()
  for(i in 1:length(tree)){
    tr <- tree[[i]]
    matched <- intersect(tr$tip.label, dat$species)
    result[[i]] <- keep.tip(tr, matched)
  }
  tree <- result
  class(tree) <- "multiPhylo"
}

# drop non-matching taxa from tree and data
dat <- dat[dat$species %in% matched, ]
dat <- dat[, c("species", "haploid")]

# convert data to matrix format required by diversitree
mat <- datatoMatrix(x = dat, buffer = 1, hyper = FALSE)

if("phylo" %in% class(tree)){
  tree <- AdjustTree(tree, plants)
}
if("multiPhylo" %in% class(tree)){
  for(j in 1:length(tree)){
    tree[[j]] <- AdjustTree(tree[[j]])
    print(paste("working on tree", j))
  }
}
if("phylo" %in% class(tree)){
  lik <- make.mkn(tree = tree, states = mat, k=ncol(mat), 
                  strict=F, control=list(method="ode",root=ROOT.OBS))
  argnames(lik)
  #  costrain our model to match the dynamics of chromosome evolution
  conlik <- constrainMkn(data = mat, lik = lik, hyper=F,
                         polyploidy = F, verbose=F)
  print(argnames(conlik))
}
if("multiPhylo" %in% class(tree)){
  res <- list()
  for(j in 1:100){
    lik <- make.mkn(tree = tree[[j]], states = mat, k=ncol(mat), 
                    strict=F, control=list(method="ode",root=ROOT.OBS))
    argnames(lik)
    #  costrain our model to match the dynamics of chromosome evolution
    conlik <- constrainMkn(data = mat, lik = lik, hyper=F,
                           polyploidy = F, verbose=F)
    print(argnames(conlik))
  }
}

