### 
# load packages
library(ape)
library(paleobuddy)

### 
# read trees
trees <- read.nexus(file = "/Users/petrucci/Downloads/canidae_run-tree (2).trees")

# read trace
trace <- read.delim(file = "filename", delim = "delim")

# read occurrences csv
occurrences <- read.csv("/Users/petrucci/Documents/research/canids_chap2/data/raw_data/canidae_occurrences.csv")

