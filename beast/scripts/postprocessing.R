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
trees <- read.nexus(file = paste0(output_dir, "tree_big.trees"))

# read trace
trace <- read.delim(file = paste0(output_dir, "canidae_run_big.log"), 
                    sep = "\t", comment.char = "#")

# read occurrences csv
occurrences <- read.csv(paste0(base_dir, 
                               "beast/data/raw_data/canidae_occurrences.csv"))

# read morphological data to get extant species
mol <- read.nexus.data(paste0(base_dir, "beast/data/raw_data/canidae_mol.nex"))
ext_taxa <- names(mol)[-which(names(mol) %in% c("Aenocyon_dirus", 
                                                "Dusicyon_australis"))]
# read diet data
#diet <- read.csv(file = "filename")
diet_raw <- read.delim(paste0(base_dir, "revbayes/data/diet.txt"),
                   sep = " ")
diet <- diet_raw[, 1] - 1
names(diet) <- rownames(diet_raw)

###
# make occurrences into a paleobuddy fossil-record object

# delete unnecessary columns
occs <- occurrences[, -c(1, 3, 6:7)]

# add extant/extinct column
occs$Extant <- unlist(lapply(1:nrow(occurrences), 
                             function(x) occurrences$taxon[x] %in% ext_taxa))

# reorder columns
occs <- occs[, c(1, 4, 2, 3)]

# change column names
colnames(occs) <- c("Species", "Extant", "MaxT", "MinT")

# number of fossils per species
n_occs <- lapply(unique(occs$Species), function(x) sum(occs$Species == x))
names(n_occs) <- unique(occs$Species)
n_occs <- as.data.frame(n_occs)

# add extant species that aren't already there
no_fossil_taxa <- ext_taxa[!(ext_taxa %in% colnames(n_occs))]
for (taxa in no_fossil_taxa) {
  n_occs[, no_fossil_taxa] <- 0
}
n_occs <- t(n_occs)

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
      # find node corresponding to the tip of this species
      cur_node <- which(tree$tip.label == fossils$Species[i])
      
      # check the time for the tip
      t_tip <- root_time - node.depth.edgelength(tree)[cur_node]
      
      # if it is not extant, subtract t_tip to have 
      # a relative time for the fossil
      t_fossil <- t_fossil - t_tip
      
      # find edge leading to this tip
      cur_edge <- which(tree$edge[, 2] == cur_node)
      
      # and get its length
      cur_edge_len <- tree$edge.length[cur_edge]
      
      # if cur_edge_len > t_fossil, the fossil goes on this edge, so we just
      # use cur_node for this fossil; if not, we need to find the right edge
      while (cur_edge_len < t_fossil) {
        # get new relative fossil time 
        t_fossil <- t_fossil - cur_edge_len
        
        # get new node (node right before edge that leads to previous)
        cur_node <- tree$edge[cur_edge, 1]
        
        # find edge again
        cur_edge <- which(tree$edge[, 2] == cur_node)
        
        # and get its length again
        # repeat until cur_edge_len > t_fossil
        cur_edge_len <- tree$edge.length[cur_edge]
      }
      
      # add the last node we visited to the res vector
      res <- c(res, cur_node)
    }
  }
  
  return(res)
}

###
# final manipulation of tree and fossils, and then
# use FossilSim to place fossils on tree

# create new tree to avoid modifying mptree
fs_tree <- mptree

# data frame of which numbered tip corresponds to which name
names_tipnum <- data.frame(names = mptree$tip.label,
                           tipnum = 1:length(mptree$tip.label))

# change tip labels to tN to match FossilSim expectation
fs_tree$tip.label <- paste0("t", unlist(lapply(fs_tree$tip.label, 
                        function(x) names_tipnum[names_tipnum$names == x, 2])))

# run fossil_edges to find which edges each fossil goes on
f_edges <- fossil_edges(mptree, fossils)

# make FossilSim fossil object with species as numbers and edge numbers
# for each fossil specimen
fs_fossils <- data.frame(sp = unlist(lapply(fossils$Species, 
                        function(x) names_tipnum[names_tipnum$names == x, 2])),
                         edge = f_edges,
                         hmin = fossils$SampT,
                         hmax = fossils$SampT)
fs_fossils <- as.fossils(fs_fossils, from.taxonomy = TRUE)

# get tree with fossils
satree_obj <- SAtree.from.fossils(fs_tree, fs_fossils)
satree <- satree_obj$tree

# vector to hold diet values for each tip
diet_satree <- c()

# manipulate tip labels back to named species
for (i in 1:length(satree$tip.label)) {
  # get tip label
  label <- satree$tip.label[i]
  
  # get species number
  label_num <- as.numeric(gsub("_[0-9]+", "", gsub("t", "", label)))
  
  # get species name
  label_species <- names_tipnum[label_num, 1]
  
  # add diet info to vector
  diet_satree <- c(diet_satree, diet[label_species])
  
  # substitute number in label for species name
  satree$tip.label[i] <- gsub("t[0-9]+_", paste0(label_species, "_"), label)
}

# name diet vector
names(diet_satree) <- satree$tip.label 

# # find tips that don't have data
# nodiet_tipnums <- which(is.na(diet_satree))

# # get final tree that doesn't have the tips without data
# satree_final <- drop.tip(satree, nodiet_tipnums)
# 
# # do the same with diet
# diet_final <- diet_satree[-nodiet_tipnums]
# 
# correct edge lengths for extant species so they're all at 0
for (i in 1:length(ext_taxa[which(ext_taxa %in% names(diet))])) {
  # check which tip it is
  tipN <- which(satree$tip.label == paste0(ext_taxa[i], "_", n_occs[ext_taxa[i], ] + 1))

  # check age of tip
  tip_age <- max(node.depth.edgelength(satree)) -
    node.depth.edgelength(satree)[tipN]

  # if it is not at 0
  if (tip_age != 0) {
    # find edge for this tip
    tip_edge <- which(satree$edge[, 2] == tipN)

    # increase length by tip_age to make it 0
    satree$edge.length[tip_edge] <- satree$edge.length[tip_edge] + tip_age
  }
}

# save diet data
write.nexus.data(diet_satree, paste0(base_dir, "revbayes/data/canidae_diet.nex"),
                 format = "standard")

# correct node labels
satree$node.label <- (length(satree$tip.label) + 1):(2*length(satree$tip.label) - 1)

# and save tree
write.nexus(satree, file = paste0(base_dir, "revbayes/data/canidae_tree.nex"),
            translate = TRUE)

# get extant tree
satree_extant <- drop.fossil(satree)
write.nexus(satree_extant, file = paste0(base_dir, "revbayes/data/canidae_tree_extant.nex"),
            translate = TRUE)

# select those data
diet_extant <- diet_satree[satree_extant$tip.label]
write.nexus.data(diet_extant, paste0(base_dir, "revbayes/data/canidae_diet_extant.nex"),
                 format = "standard")
