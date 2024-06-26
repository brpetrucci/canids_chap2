### 
# load packages
library(ape)
library(paleobuddy)

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
diet <- read.csv(file = "filename")

###
# make occurrences into a paleobuddy fossil-record object

# delete unnecessary columns
occs <- occurrences[, -c(1, 3, 6:7)]

# add extant/extinct column
occs$Extant <- occurrences$late_age == 0

# reorder columns
occs <- occs[, c(1, 4, 2, 3)]

# change column names
colnames(occs) <- c("Species", "Extant", "MaxT", "MinT")

###
# manipulate maximum posterior tree

# max posterior sample
mps <- which(trace$posterior == max(trace$posterior))

# then need to manually look for the tree
# on trees.tree and copy it to a file since ape
# deletes the text in the file that we need

# read MP tree file
mptree <- read.nexus(file = paste0(base_dir, "beast/mptree.nex"))

# find tips that are first ages
fa_tips <- mptree$tip.label[grepl("*first", mptree$tip.label)]

# remove them from tree
mptree <- drop.tip(mptree, fa_tips)

# remove last from tip labels
mptree$tip.label <- gsub("_last", "", mptree$tip.label)
