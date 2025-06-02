###
# packages

# devtools
library(devtools)

# rBt
library(tracerer)

# ape
library(ape)

# phytools
library(phangorn)

# treeio
library(treeio)

###
# prep

# set directories
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/"

# source clades R script
source(paste0(base_dir, "clades.R"))

###
# make new sets of posterior trees taking out extra metadata and terminals

# function to make SRFBD and FBDS tree outputs 
# into normal outputs (one terminal/species)
make_new_posterior <- function(dir, infile, outfile, 
                               mcc_file = NULL, extant = FALSE) {
  # read the Nexus file as lines
  tree_lines <- readLines(paste0(dir, infile))
  
  # parse translate block
  translate_start <- grep("^\\s*TRANSLATE", tree_lines, ignore.case = TRUE)
  translate_end <- grep(";", tree_lines[(translate_start + 1):length(tree_lines)], 
                        fixed = TRUE)[1] + translate_start
  translate_block <- tree_lines[translate_start:translate_end]
  
  # parse tree strings
  tree_strings <- grep("^\\s*tree ", tree_lines, value = TRUE)
  
  tree_names <- sub("tree ([^ ]+) =.*", "\\1",
                    grep("^\\s*TREE\\s", tree_lines, 
                         ignore.case = TRUE, value = TRUE))
  
  # extract list of treedata 
  tree_list <- lapply(seq_along(tree_strings), function(i) {
    # get string for this tree
    tree_str <- tree_strings[i]
    
    # create a temporary file
    tmp_file <- tempfile(fileext = ".nex")
    
    # save tree to temporary file
    cat("#NEXUS\n\nBEGIN TREES;\n",
        paste(translate_block, collapse = "\n"), "\n",
        tree_str, "\nEND;\n",
        file = tmp_file)
    
    # read the tree from the temporary file
    treedata <- read.beast(tmp_file)
    
    # reorder data based on the nodes column
    treedata@data <- treedata@data[order(as.numeric(treedata@data$node)), ]
    
    # add range information to tips
    #treedata@data$range[1:length(treedata@phylo$tip.label)] <-
    #  sub('_first|_last|_ext', '', treedata@phylo$tip.label)
    
    # return treedata
    treedata
  }) 
  
  # get the tips to drop for the first tree (it will be the same for all)
  tree <- tree_list[[1]]
  
  # start a vector to hold tips to be dropped
  if (extant) {
    # check if fbds analysis
    if (any(grepl("_ext", tree@phylo$tip.label))) {
      # if so, need to get extant species as _ext
      species <- c(paste0(extant_species, "_ext"),
                   "Dusicyon_australis_last", "Aenocyon_dirus_last")
    } else {
      # if not, need to check which ones are singletons
      ext_singles <- recent_species[which(!(paste0(recent_species, "_last") %in% 
                                              tree@phylo$tip.label))]
      
      # if not, just need the recent species with _last
      species <- c(paste0(recent_species[!(recent_species %in% ext_singles)], 
                          "_last"),
                   paste0(ext_singles, "_first"))
    }
    
    # make drop_tips all other tips
    drop_tips <- tree@phylo$tip.label[!(tree@phylo$tip.label %in% species)]
  } else {
    # otherwise, gotta check deeper on what to drop
    drop_tips <- unlist(lapply(seq_along(tree@phylo$tip.label), function(i) {
      # get tip
      tip <- tree@phylo$tip.label[i]
      
      # check whether it's a first or last occurrence
      first <- grepl('_first', tip)
      last <- grepl('_last', tip)
      
      # get the taxon for this tip
      tx <- gsub('(_first|_last|_ext)', '', tip)
      
      # if it is the first and there's a last or ext, or if
      # it is the last and there's an ext, drop it
      if ((first && ((paste0(tx, "_last") %in% tree@phylo$tip.label) ||
                     (paste0(tx, "_ext") %in% tree@phylo$tip.label))) ||
          (last && paste0(tx, "_ext") %in% tree@phylo$tip.label)) {
        tip
      }
    }))
  }

  # get singletons
  singles <- tree@phylo$tip.label[grepl("_first", tree@phylo$tip.label) &
                                    !(tree@phylo$tip.label %in% drop_tips)]
  
  # if mcc is a file, we find and save an MCC-ish tree
  if (!is.null(mcc_file)) {
    # get multiPhylo object with trees only
    trees <- lapply(seq_along(tree_list), function(i) tree_list[[i]]@phylo)
    class(trees) <- "multiPhylo"
    
    # get mcc tree
    mcc_tree <- mcc(trees)
    
    # get index for tree from sample that has the same topology
    idx_mcc_tree <- min(which(unlist(lapply(seq_along(trees), function(i) {
      all.equal(trees[[i]], mcc_tree, use.edge.length = FALSE
                )}))))
    
    # and get the treedata
    mcc_treedata <- tree_list[[idx_mcc_tree]]
    
    # get posterior clade probabilities
    post_clades <- as.numeric(mcc_tree$node.label)
    
    # add to data
    mcc_treedata@data$posterior <- c(rep(NA, 
                                         length(mcc_treedata@phylo$tip.label)),
                                     post_clades)
    
    # make new column to name tips
    mcc_treedata@data$species <- mcc_treedata@data$range
    
    # add species info to singleton tips
    mcc_treedata@data$species[which(mcc_treedata@phylo$tip.label %in% singles)] <-
      sub("_first", "", singles)
    
    # write to mcc_file
    write.beast(mcc_treedata, paste0(dir, mcc_file),
                tree.name = "mcc_tree")
  }
  
  tree_list_new <- lapply(seq_along(tree_list), function(i) {
    # get the tree
    tree <- tree_list[[i]]
    
    # get the tree without these tips
    tree_clean <- drop.tip(tree, drop_tips)
    
    # clean the name up
    tree_clean@phylo$tip.label <- gsub('(_first|_last|_ext)', '', 
                                       tree_clean@phylo$tip.label)
    
    # rename nodelabels
    tree_clean@phylo$node.label <- 
      as.character(seq_len(tree_clean@phylo$Nnode) + 
                     length(tree_clean@phylo$tip.label))
    
    # remove range and orientation data, if there are any
    if ("range" %in% colnames(tree_clean@data)) {
      tree_clean@data <- tree_clean@data[, !(colnames(tree_clean@data) %in% 
                                               c("range", "orientation"))]
    }
    
    # return the new tree
    tree_clean
  })
  
  # make new translation block
  tip_labels <- tree_list_new[[1]]@phylo$tip.label
  tip_translation <- paste(seq_along(tip_labels), tip_labels, sep = " ", 
                           collapse = ",\n\t\t")
  translate_block <- sprintf("\tTRANSLATE\n\t\t%s;\n", tip_translation)
  
  # make the new tree_strings list
  tree_strings <- sapply(seq_along(tree_list_new), function (i) {
    # get this tree
    tree <- tree_list_new[[i]]
    
    # make the tip labels numbers
    tree@phylo$tip.label <- seq_along(tip_labels)
    
    # get the newick
    con <- textConnection("OUT", "w", local = TRUE)
    write.beast.newick(tree, file = con)
    close(con)
    tree_newick <- get("OUT")
    
    # insert it into the tree strings
    sprintf("tree %s = %s", tree_names[i], tree_newick[1])
  })
  
  # write final nexus file
  cat("#NEXUS\n\nBEGIN TREES;\n",
      paste(translate_block, collapse = "\n"), "\n",
      paste(tree_strings, collapse = "\n"), "\nEND;\n",
      file = paste0(dir, outfile))
  
  # return the cleaned up tree list
  return(tree_list_new)
}

##
# SRFBD - run that converged the most

# directory
srfbd_dir <- paste0(base_dir, "srfbd/")

# input file
srfbd_infile <- "output/new.trees"

# output file for complete trees
srfbd_complete_outfile <- "output/mod_new.trees"

# file for MCC tree
mcc_file <- "output/mcc_new.nex"

# make new posterior (complete trees)
srfbd_complete_trees <- make_new_posterior(srfbd_dir, srfbd_infile, 
                                           srfbd_complete_outfile, mcc_file)
### remember to add range info to nodes missing it (e.g. V. velox and U. cinereoargenteus)

# output file for extant trees
srfbd_extant_outfile <- "output/mod_new_ftrd_extant.trees"

# make new posterior (extant trees)
srfbd_extant_trees <- make_new_posterior(srfbd_dir, srfbd_infile, 
                                         srfbd_extant_outfile, 
                                         extant = TRUE)
 
##
# FBD specimen (both first and last occurrences)

# directory
fbds_dir <- paste0(base_dir, "fbds/")

# input file
fbds_both_infile <- "output/both.trees"

# output file for complete trees
fbds_both_complete_outfile <- "output/mod_both.trees"

# make new posterior (complete trees)
fbds_both_complete_trees <- make_new_posterior(fbds_dir, fbds_both_infile, 
                                               bds_both_complete_outfile)

# output file for extant trees
fbds_both_extant_outfile <- "output/mod_both_extant.trees"

# make new posterior (extant trees)
fbds_both_extant_trees <- make_new_posterior(fbds_dir, fbds_both_infile, 
                                             fbds_both_complete_outfile, 
                                             extant = TRUE)

##
# FBD specimen (just first occurrence)

# input file
fbds_first_infile <- "output/first.trees"

# output file for complete trees
fbds_first_complete_outfile <- "output/mod_first.trees"

# make new posterior (complete trees)
fbds_first_complete_trees <- make_new_posterior(fbds_dir, fbds_first_infile, 
                                                fbds_first_complete_outfile)

# output file for extant trees
fbds_first_extant_outfile <- "output/mod_first_extant.trees"

# make new posterior (extant trees)
fbds_first_extant_trees <- make_new_posterior(fbds_dir, fbds_first_infile,
                                              fbds_first_extant_outfile,
                                              extant = TRUE)

###
# get clade probabilities and divergence time posterior distributions

trees <- lapply(seq_along(srfbd_complete_trees), function(i) {
  srfbd_complete_trees[[i]]@phylo
})
class(trees) <- "multiPhylo"

# function to check the posterior probability of a clade for a given tree sample
clade_prob <- function(trees, clade, clades, conds) {
  # get clade it's conditioned on
  cond <- conds[clade]
  
  # take only trees where condition is met
  trees_cond <- trees[unlist(lapply(seq_along(trees), 
                                    function(i) is.monophyletic(trees[[i]], clades[[cond]])))]
  
  # get clade probability in isolation
  prob <- sum(unlist(lapply(1:length(trees), 
                    function(x) {
                      is.monophyletic(trees[[x]], clades[[clade]])
                      }))) / length(trees)
  
  # and within trees_cond trees only
  prob_cond <- sum(unlist(lapply(1:length(trees_cond), 
                              function(x) {
                               is.monophyletic(trees_cond[[x]], clades[[clade]])
                   }))) / length(trees_cond)
  
  # probability of the condition
  prob_cond_clade <- sum(unlist(lapply(1:length(trees), 
                                 function(x) {
                                   is.monophyletic(trees[[x]], clades[[cond]])
                                 }))) / length(trees)
  
  ### gotta think what I'll do about conditions...
  
}


###
# retired

# # read logs from runs and thin them
# read_logs <- function(dir, log_filename, n_logs, burnin, thin) {
#   # read first log file
#   log <- read.delim(paste0(dir, "output/1_", log_filename),
#                     comment.char = "#")
# 
#   # start final log data frame
#   final_log <- data.frame(matrix(nrow = 0, ncol = ncol(log)))
# 
#   # start index list
#   idx_list <- vector("list", n_logs)
# 
#   # iterate through number of logs
#   for (i in 1:n_logs) {
#     # read first log file
#     if (i > 1) {
#       log <- read.delim(paste0(srfbd_dir, "output/", i, "_", log_filename),
#                         comment.char = "#")
#     }
# 
#     # burn it
#     #burnt_idx <- round((burnin * nrow(log)), 0):nrow(log)
# 
#     # thin it
#     #thinned_idx <- sample(burnt_idx, round(thin * length(burnt_idx)))
# 
#     # add just those indices to the final log
#     #final_log <- rbind(final_log, log[thinned_idx, ])
#     final_log <- rbind
#     # add idx to idx_list
#     idx_list[[i]] <- thinned_idx
#   }
# 
#   # return log and idx_list
#   return(list(LOG = final_log, IDX = idx_list))
# }
# 
# # read trees with the indices from log
# read_trees <- function(dir, tree_filename, idx_list) {
#   # get a final tree sample list
#   trees_final <- c(rtree(1))
#   class(trees_final) <- "multiPhylo"
#   
#   # iterate through number of tree files
#   for (i in 1:length(idx_list)) {
#     # read trees
#     trees <- parse_beast_trees(paste0(dir, "output/", i, "_", 
#                                       tree_filename))
#     
#     # select only the ones that we want from idx_list
#     trees_iter <- trees[idx_list[[i]]]
#     
#     # add it to trees final
#     trees_final <- c(trees_final, trees_iter)
#   }
#   
#   # remove first object of trees_final since it was just a sim tree
#   trees_final <- trees_final[-1]
#   
#   # return trees_final
#   return(trees_final)
# }
# 
# # function to return trees and log for a given burnin and thin
# read_all <- function(dir, filename, tree_suffix, n_logs, burnin, thin) {
#   
#   # filename for analysis
#   log_filename <- paste0(filename, ".log")
#   
#   # get log and idx
#   log_idx <- read_logs(dir, log_filename, n_logs, burnin, thin)
#   log <- log_idx$LOG
#   idx <- log_idx$IDX
#   
#   # filename for ftrd trees
#   tree_filename <- paste0(filename, tree_suffix, ".trees")
#   
#   # get tree sample
#   trees <- read_trees(dir, tree_filename, idx)
#   
#   # return log and trees
#   return(list(LOG = log, TREES = trees, IDX = idx))
# }
# 
# # extract just extant trees
# make_extant <- function(trees, type = "srfbd") {
#   # iterate through tree object
#   for (i in 1:length(trees)) {
#     # get this tree
#     tree <- trees[[i]]
#     
#     # remove all fossils
#     tree_ext <- tree_ext <- drop.fossil(tree)
#     
#     # if type is srfbd, remove just _first
#     if (type == "srfbd") {
#       # get tips with _first
#       first_tips <- grep("_first", tree_ext$tip.label)
#       
#       # drop them
#       tree_ext <- drop.tip(tree_ext, tree_ext$tip.label[first_tips])
#       
#       # remove _last from tip labels
#       tree_ext$tip.label <- sub("_last", "", tree_ext$tip.label)
#       
#       # set trees element to this new tree
#       trees[[i]] <- tree_ext
#     }
#   }
#   
#   # name them
#   names(trees) <- paste0("tree_", 1:length(trees))
#   
#   # return trees
#   return(trees)
# }
# 
# # function to get most likely topology
# get_MAP_top <- function(trees) {
#   # get topology of first tree
#   tops <- c(trees[[1]])
#   class(tops) <- "multiPhylo"
#   
#   # number of each topology
#   nums <- c(1)
#   
#   # iterate through rest of the trees
#   for (i in 2:length(trees)) {
#     # get this tree
#     tree <- trees[[i]]
#     
#     # check if found
#     found <- FALSE
#     
#     # compare with each of the tops already found
#     for (j in 1:length(tops)) {
#       if (all.equal.phylo(tree, tops[[j]], use.edge.length = FALSE)) {
#         # increase number of tops in this
#         nums[j] <- nums[j] + 1
#         
#         # it was found
#         found <- TRUE
#         
#         break
#       }
#     }
#     
#     if (!found) {
#       # if we never found it, add it to tops
#       tops <- c(tops, tree)
#       
#       # and nums
#       nums <- c(nums, 1)
#     }
#   }
#   
#   return(list(tops = tops, nums = nums))
# }
# 
# # function to find stem lineage
# find_stem <- function(tree, sps) {
#   # get root node
#   root_node <- max(which(node.depth.edgelength(tree) == 0))
#   
#   # get two sisters from first split
#   sisters <- Children(tree, root_node)
#   
#   # check sizes of each sister
#   sister_sizes <- sapply(sisters, function(node) 
#     length(Descendants(tree, node, type = "tips")[[1]]))
#   
#   # get the smaller one, that's the stem
#   stem <- tree$tip.label[Descendants(tree, sisters[which.min(sister_sizes)], type = "tips")[[1]]]
#   
#   # return it as text
#   paste(stem, collapse = "/")
# }
# 
# # get stem for each tree in a sample
# find_stems <- function(trees, sps) {
#   # create results vector
#   stems <- c()
#   
#   # iterate through trees
#   for (i in 1:length(trees)) {
#     print(i)
#     # get this tree
#     tree <- trees[[i]]
#     
#     # find stem
#     stems <- c(stems, find_stem(tree, sps))
#     
#     # # get tips
#     # tips <- tree$tip.label
#     # 
#     # # species labels
#     # species <- unique(sub("_first", "", 
#     #                       sub("_last", "", tips)))
#     # 
#     # # iterate through tip labels
#     # for (j in 1:length(species)) {
#     #   # get this species
#     #   sp <- species[j]
#     #   
#     #   # get first tip for species
#     #   first_tip <- paste0(sp, "_first")
#     #   
#     #   # get both
#     #   both_tips <- paste0(sp, c("_first", "_last"))
#     #   
#     #   # check if stem is within this species
#     #   if (is.monophyletic(tree, tips[tips != first_tip])) {
#     #     # if so, add to stem and break
#     #     stems <- c(stems, first_tip)
#     #     break
#     #   }
#     #   # check if stem is the lineage including this species
#     #   else if (is.monophyletic(tree, tips[-which(tips %in% both_tips)])) {
#     #     # add species to stems and break
#     #     stems <- c(stems, sp)
#     #     break
#     #   }
#     # }
#   }
#   
#   # return stems
#   return(stems)
# }