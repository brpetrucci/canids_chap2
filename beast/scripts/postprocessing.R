### 
# load packages
library(ape)
library(paleobuddy)
library(FossilSim)

### 
# load files

# directories
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/"
output_dir <- paste0(base_dir, "beast/output/")

# read trees
trees <- read.nexus(file = paste0(output_dir, "tree.trees"))

# read trace
trace <- read.delim(file = paste0(output_dir, "canidae_run.log"), 
                    sep = "\t", comment.char = "#")

# read occurrences csv
occurrences <- read.csv(paste0(base_dir, 
                               "beast/data/raw_data/canidae_occurrences.csv"))

# read diet data
#diet <- read.csv(file = "filename")

###
# make occurrences into a paleobuddy fossil-record object

# delete unnecessary columns
occs <- occurrences[, -c(1, 3, 6:7)]

# add extant/extinct column
occs$Extant <- unlist(lapply(1:nrow(occurrences), 
                             function(x) any(occurrences[occurrences$taxon == occurrences$taxon[x], ]$late_age == 0)))

# reorder columns
occs <- occs[, c(1, 4, 2, 3)]

# change column names
colnames(occs) <- c("Species", "Extant", "MaxT", "MinT")

# number of fossils per species
n_occs <- lapply(unique(occs$Species), function(x) sum(occs$Species == x))
names(n_occs) <- unique(occs$Species)
n_occs <- t(as.data.frame(n_occs))

###
# manipulate maximum posterior tree and fossils data frame

# max posterior sample
mps <- which(trace$posterior == max(trace$posterior))

# then need to manually look for the tree
# on trees.tree and copy it to a file since ape
# deletes the text in the file that we need

# read MP tree file
mptree <- read.nexus(file = paste0(output_dir, "mptree.nex"))

# find tips that are first ages
fa_tips <- mptree$tip.label[grepl("*first", mptree$tip.label)]

# vector for first ages
fa_sps <- c()

# make a list of first ages for these species
for (i in 1:length(fa_tips)) {
  # find the tip number for this tip
  tip_num <- which(mptree$tip.label == fa_tips[i])
  
  # get age for this tip
  fa_sp <- max(node.depth.edgelength(mptree)) - 
    node.depth.edgelength(mptree)[tip_num]
  
  # add it to vector
  fa_sps <- c(fa_sps, fa_sp)
}

# name it with species names
names(fa_sps) <- gsub("_first", "", fa_tips)

# remove them from tree
mptree <- drop.tip(mptree, fa_tips)

# remove last from tip labels
mptree$tip.label <- gsub("_last", "", mptree$tip.label)

# create final fossils data frame
fossils <- data.frame(matrix(nrow = 0, ncol = 5))

# remove fossils that are already in the tree, i.e. the last occurrence
for (i in 1:length(unique(occs$Species))) {
  # get occurrences for this species
  sp_occs <- occs[occs$Species == unique(occs$Species)[i], ]
  
  # find min occurrence
  min_t <- min(sp_occs$MinT)
  
  # which to delete (last one)
  occ_del <- max(which(sp_occs$MinT == min_t))
  
  # delete from sp_occs if species is extinct
  if (sp_occs$Extant[1] == FALSE) {
    sp_occs <- sp_occs[-occ_del, ]
  }
  
  # vector of sampling times
  samp_times <- c()
  
  if (nrow(sp_occs) > 0) {
    # for each occurrence, sample a time until it is 
    # less than first age, more than last age
    for (occ in 1:nrow(sp_occs)) {
      # which tip is this in the tree
      tip_num <- which(mptree$tip.label == sp_occs$Species[1])
      
      # get fa and la
      fa <- fa_sps[sp_occs$Species[1]]
      la <- max(node.depth.edgelength(mptree)) - 
        node.depth.edgelength(mptree)[tip_num]
      
      # initialize
      occ_time <- Inf
      
      # iterate until it's correct
      while (occ_time > fa || occ_time < la) {
        # draw
        occ_time <- runif(1, la, fa)
      }
      
      # add it to samp_times
      samp_times <- c(samp_times, occ_time)
    }
    
    # make a SampT field in sp_occs
    sp_occs$SampT <- samp_times
  }
  
  # add sp_occs to fossils
  fossils <- rbind(fossils, sp_occs)
}

# delete MaxT and MinT
fossils <- fossils[, -c(2, 3, 4)]

###
# function to identify which edge each fossil should go in
fossil_edges <- function(tree, fossils) {
  # time of the root of the tree
  root_time <- max(node.depth.edgelength(tree))
  
  # initialize results vector (nodes after)
  res <- c()
  
  # iterate through fossil occurrences
  for (i in 1:nrow(fossils)) {
    # get the time for this fossil
    t_fossil <- fossils$SampT[i]
    
    # if it happened before the root, give first non-tip node
    if (t_fossil > root_time) {
      res <- c(res, length(tree$tip.label) + 1)
      
    } else {
      # find node for the current 
      cur_node <- which(tree$tip.label == fossils$Species[i])
      t_tip <- root_time - node.depth.edgelength(tree)[cur_node]
      t_fossil <- t_fossil - t_tip
      
      cur_edge <- which(tree$edge[, 2] == cur_node)
      cur_edge_len <- tree$edge.length[cur_edge]
      
      while (cur_edge_len < t_fossil) {
        t_fossil <- t_fossil - cur_edge_len
        cur_node <- tree$edge[cur_edge, 1]
        cur_edge <- which(tree$edge[, 2] == cur_node)
        cur_edge_len <- tree$edge.length[cur_edge]
      }
      
      res <- c(res, cur_node)
    }
  }
  
  return(res)
}

###
# final manipulation of tree and fossils, and then
# use FossilSim to place fossils on tree

# create new tree 
fs_tree <- mptree
names_tipnum <- data.frame(names = mptree$tip.label,
                           tipnum = 1:length(mptree$tip.label))
fs_tree$tip.label <- paste0("t", unlist(lapply(fs_tree$tip.label, 
                        function(x) names_tipnum[names_tipnum$names == x, 2])))

f_edges <- fossil_edges(mptree, fossils)
fs_fossils <- data.frame(sp = unlist(lapply(fossils$Species, 
                        function(x) names_tipnum[names_tipnum$names == x, 2])),
                         edge = f_edges,
                         hmin = fossils$SampT,
                         hmax = fossils$SampT)
fs_fossils <- as.fossils(fs_fossils, from.taxonomy = TRUE)

satree <- SAtree.from.fossils(fs_tree, fs_fossils)$tree

# need to 1. switch names of tips to tN where N is the number on tip.label;
# 2. do the same for fossils but with just the number; 3. run the edges
# function to get the node before which a fossil is to be placed; 4. make 
# fossils into a FossilSim fossils object; 5. run SAtree.from.fossils to get
# the SA tree; 6. switch the names on the tree back to before, so t1_1 should
# be outgroup_1