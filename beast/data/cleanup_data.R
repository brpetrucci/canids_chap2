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
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/data/"

# raw data directory
raw_dir <- paste0(base_dir, "raw_data/")

# fossil occurrences
occurrences <- read.csv(paste0(raw_dir, 
                               "canidae_occurrences.csv"))[, -1]

# molecular data
mol <- read.nexus.data(paste0(raw_dir, "canidae_mol.nex"))

# morphological data
morpho_raw <- read.nexus.data(paste0(raw_dir, "canidae_morpho.nex"))

# cut state 76 since it is non-variable
morpho <- lapply(morpho_raw, function(x) x[-76])

# cut Vulpes bengalensis since it has so little data
# cut outgroup since we decided we didn't need it
morpho <- morpho[-(which(names(morpho) == "Vulpes_bengalensis" 
                         | names(morpho) == "outgroup"))]

###
# extract range data

# function returning max and min for each range for a given taxon name
get_range <- function(occurrences, name) {
  # get all occurrences with that name
  named_occs <- occurrences[occurrences$taxon == name, ]
  
  # number of occurrences
  n_occs <- nrow(named_occs)
  
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
  
  return(c(name, fa_max, fa_min, la_max, la_min, n_occs))
}

# apply to all species
ranges <- lapply(unique(occurrences$taxon), 
                 function(x) get_range(occurrences, x))

# make it a data frame
ranges_df <- t(as.data.frame(ranges))
colnames(ranges_df) <- c("taxon", "fa_max", "fa_min", "la_max", "la_min", "k")
rownames(ranges_df) <- 1:nrow(ranges_df)

# add data without fossil occurrences to ranges
taxa_names <- unique(c(names(morpho), names(mol)))
nofossil_taxa <- taxa_names[!(taxa_names %in% ranges_df[, 1])]
ranges_df <- rbind(ranges_df, data.frame(taxon = nofossil_taxa,
                                         fa_max = rep(0, length(nofossil_taxa)),
                                         fa_min = rep(0, length(nofossil_taxa)),
                                         la_max = rep(0, length(nofossil_taxa)),
                                         la_min = rep(0, length(nofossil_taxa)),
                                         k = 0))
rownames(ranges_df) <- 1:nrow(ranges_df)

# write ranges to a file
write.table(ranges_df, paste0(base_dir, "canidae_ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# add ?? for molecular data for species without molecular data
nomol_taxa <- taxa_names[!(taxa_names %in% names(mol))]
nomol_mol <- lapply(nomol_taxa, function(x) rep("?", length(mol[[1]])))
names(nomol_mol) <- nomol_taxa
mol_complete <- c(mol, nomol_mol)

# write DNA to a file
write.nexus.data(mol_complete, paste0(base_dir, "canidae_mol.nex"))

# reading partitions and writing to DNA file
partitions <- readLines(paste0(base_dir, "partitions_simplified.txt"))
write(partitions, paste0(base_dir, "canidae_mol.nex"), append = TRUE)

# add ?? for morphological data for species without morphological data
nomorpho_taxa <- taxa_names[!(taxa_names %in% names(morpho))]
nomorpho_morpho <- lapply(nomorpho_taxa, function(x) rep("?", length(morpho[[1]])))
names(nomorpho_morpho) <- nomorpho_taxa
morpho_complete <- c(morpho, nomorpho_morpho)

# write DNA to a file
write.nexus.data(morpho_complete, paste0(base_dir, "canidae_morpho.nex"),
                 format = "standard")

###
# creating a data set without occurrences with >5my uncertainty

# filter occurrences
occurrences_filter <- (occurrences$early_age - occurrences$late_age) > 5
ftrd_occs <- occurrences[!occurrences_filter, ]

# filter morpho
ftrd_names <- unique(occurrences$taxon)[!(unique(occurrences$taxon) %in% 
                                            ftrd_occs$taxon)]
ftrd_morpho <- morpho_complete[-which(names(morpho_complete) %in% ftrd_names)]
ftrd_mol <- mol_complete[-which(names(mol_complete) %in% ftrd_names)]

# filter ranges
ftrd_ranges_df <- ranges_df[-which(ranges_df[, 1] %in% ftrd_names), ]
rownames(ftrd_ranges_df) <- 1:nrow(ftrd_ranges_df)

# write ranges and morpho data
write.table(ftrd_ranges_df, paste0(base_dir, "ftrd_canidae_ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.nexus.data(ftrd_morpho, paste0(base_dir, "ftrd_canidae_morpho.nex"),
                 format = "standard")
write.nexus.data(ftrd_mol, paste0(base_dir, "ftrd_canidae_mol.nex"))
write(partitions, paste0(base_dir, "ftrd_canidae_mol.nex"), append = TRUE)

