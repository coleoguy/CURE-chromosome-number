library(chromePlus)
library(diversitree)

# get all the chromosome data files
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = TRUE)
GetTree <- function(x){
  tree <- NULL
  if(file.exists(paste0("../data/trees/",x,".new"))){
    tree <- read.tree(paste0("../data/trees/",x,".new"))
  }
  if(file.exists(paste0("../data/trees/",x, ".nex"))){
    tree <- read.nexus(paste0("../data/trees/",x,".nex"))
  }
  tree <- GetCleanTree(tree)
  return(tree)
}


GetCleanTree <- function(tree){
  duptips <- c()
  if(class(tree) == "phylo"){
    valnames <- unique(tree$tip.label)
    for(i in 1:length(valnames)){
      duptips <- c(duptips, which(valnames[i] == tree$tip.label)[-1])
    }
    tree <- drop.tip(tree, duptips)
  }else{
    trees <- list()
    for(j in 1:length(tree)){
      valnames <- unique(tree[[j]]$tip.label)
      for(i in 1:length(valnames)){
        duptips <- c(duptips, which(valnames[i] == tree[[j]]$tip.label)[-1])
      }
      trees[[j]] <- drop.tip(tree[[j]], duptips)
    }
    tree <- trees
  }
  return(tree)
}



# this loop will go through all of the clades and 
# check various items
for(i in 1:length(file_list)){
  dat <- read.csv(file_list[i])
  checker <- all(c("species","haploid") %in% colnames(dat))
  clade <- strsplit(basename(file_list[i]), split=".", fixed=T)[[1]][1]
  if(!checker){
    #print(clade)
  }
  tree <- GetTree(clade)
  if(class(tree) == "phylo"){
    n <- sum(tree$tip.label %in% dat$species)
  }else{
    n <- sum(tree[[1]]$tip.label %in% dat$species)
  }
  print(paste(clade, ": ", i, " matches: ", n))
}


for(i in 1:length(tree$tip.label)){
  cursp <- tree$tip.label[i]
  tree$tip.label[i] <- strsplit(cursp, split = ".", fixed=T)[[1]][2]
}
write.tree(tree, file="bryophytes.new")

# now get the chromosome data
c(1,3,4,5,8,9,13,15-19,21,28,29,32,34,38,39,40,42,53,56)

# if newick tree
tree <- read.tree(tree.file)
# if nexus tree
tree <- read.nexus(tree.file)

# Enter your match number on your worksheet
sum(dat$species %in% tree$tip.label)

# now we must prune the tree and the data that do not match
# prune the tree first
drop <- !tree$tip.label %in% dat$species
tree <- drop.tip(phy=tree, tip=tree$tip.label[drop])
# prune data next
dat <- dat[dat$species %in% tree$tip.label,]

# this code deals with any duplicate records in the chromosome data
# leaves us with dat as our primary data file
chroms <- c()
for(i in 1:length(tree$tip.label)){
  cursp <- tree$tip.label[i]
  hap <- dat$haploid[dat$species == cursp]
  if(length(hap)>1){
    hap <- sample(hap, 1)
  }
  chroms[i] <- hap
}
names(chroms) <- tree$tip.label
dat <- data.frame(species=names(chroms), haploid=chroms)

# if there are any commas or dashes it still creates a problem


# convert to the matrix type that works with diversitree
mat <- datatoMatrix(x=dat, buffer=1, hyper=F)

# check and if necessary ultrametricize our tree
if(!is.ultrametric(tree)){
  tree <- chronos(phy=tree)
  plot(tree)
}

# scale tree to unit length
tree$edge.length <- tree$edge.length/max(branching.times(tree))
# check for any tiny negatives
tree$edge.length[tree$edge.length<0] <- .000000001
is.ultrametric(tree)
# this is the step where we create a likelihood function
lik <- make.mkn(tree = tree, states = mat, k=ncol(mat), 
                strict=F, control=list(method="ode",root=ROOT.OBS))
argnames(lik)

# lets costrain our model to match the dynamics of chromosome evolution
conlik <- constrainMkn(data = mat, lik = lik, hyper=F,
                       polyploidy = F, verbose=F)
# checking argumnet names
argnames(conlik)
# checking that functino works
conlik(pars=c(.1,.1,.1,.1))

res <- mcmc(lik=conlik, x.init=runif(4),prior=make.prior.exponential(2),nsteps=1, w=1)


tree <- read.nexus("data/trees/accipitriformes.nex")
write.csv(tree$tip.label, file="birdtips.csv")
foo <- read.csv("birdtips.csv",header=F)[,1]




dat <-read.csv("../data/chrome/carabidae.csv")
x <- sample(dat$haploid)
vars <- c()
for(i in 3:777){
  vars[i] <- var(x[1:i])/mean(x[1:i])
}
plot(vars)



