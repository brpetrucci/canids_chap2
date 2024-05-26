#######################################
##            Chapter 2              ##
##      Canids FBDR + SSE study      ##
##.         Data cleanup             ##
##     Bruno do Rosario Petrucci     ##
#######################################

###
# load packages

# ape
library(ape)

# PBDB
library(paleobioDB)

###
# read data

# base directory
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/data/"

# raw data directory
raw_dir <- paste0(base_dir, "raw_data/")

# fossil occurrences
occurrences <- read.csv(paste0(raw_dir, 
                               "canidae_occurrences.csv"))[, -1]

# molecular data
mol <- read.nexus.data(paste0(raw_dir, "canidae_mol.nex"))
mol_dnabin <- nexus2DNAbin(mol)

# morphological data
morpho <- read.nexus.data(paste0(raw_dir, "canidae_morpho.nex"))

###
# extract range data

# function returning max and min for each range for a given taxon name
get_range <- function(occurrences, name) {
  # get all occurrences with that name
  named_occs <- occurrences[occurrences$taxon == name, ]
  
  # get highest FA occurrences
  fa_high <- named_occs[which(named_occs$early_age ==
                                max(named_occs$early_age)), ]
  
  # max and min for fa
  fa_max <- max(fa_high$early_age)
  fa_min <- max(fa_high$late_age)
  
  # get lowest last age occurrence
  la_low <- named_occs[which(named_occs$late_age ==
                               min(named_occs$late_age)), ]
  
  # max and min for la
  la_min <- min(la_low$late_age)
  la_max <- ifelse(la_min == 0, 0, min(la_low$early_age))
  # if la_min is 0, the species is extant, so we know
  # their late age with certainty
  
  return(c(name, fa_max, fa_min, la_max, la_min))
}

# apply to all species
ranges <- lapply(unique(occurrences$taxon), 
                 function(x) get_range(occurrences, x))

# make it a data frame
ranges_df <- t(as.data.frame(ranges))
colnames(ranges_df) <- c("taxon", "fa_max", "fa_min", "la_max", "la_min")
rownames(ranges_df) <- 1:nrow(ranges_df)

# write ranges to a file
write.table(ranges_df, paste0(base_dir, "canidae_ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

###
# creating a data set without occurrences with >5my uncertainty

# filter occurrences
occurrences_filter <- (occurrences$early_age - occurrences$late_age) > 5
ftrd_occs <- occurrences[!occurrences_filter, ]

# filter morpho
ftrd_names <- unique(occurrences$taxon)[!(unique(occurrences$taxon) %in% 
                                            ftrd_occs$taxon)]
ftrd_morpho <- morpho[-which(names(morpho) %in% ftrd_names)]

# filter ranges
ftrd_ranges_df <- ranges_df[-which(ranges_df[, 1] %in% ftrd_names), ]

# write ranges and morpho data
write.table(ftrd_ranges_df, paste0(base_dir, "ftrd_canidae_ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.nexus.data(ftrd_morpho, paste0(base_dir, "ftrd_canidae_morpho.nex"),
                 format = "standard")

                                   