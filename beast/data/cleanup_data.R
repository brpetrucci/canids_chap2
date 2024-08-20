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

# palaeoverse
library(palaeoverse)

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

# extant taxa
ext_taxa <- names(mol)[-which(names(mol) %in% c("Aenocyon_dirus", 
                                                "Dusicyon_australis"))]

# morphological data
morpho_raw <- read.nexus.data(paste0(raw_dir, "canidae_morpho.nex"))

# cut Vulpes bengalensis since it has so little data
# cut outgroup since we decided we didn't need it
morpho <- morpho_raw[-(which(names(morpho_raw) == "Vulpes_bengalensis" 
                         | names(morpho_raw) == "outgroup"))]

# check which characters have no variation
no_var <- c()
for (i in 1:length(morpho[[1]])) {
  # get list of values
  vals <- unlist(lapply(morpho, function(x) x[i]))
  
  # get unique characters present
  chars <- unique(vals)
  
  # if there is only one non-? character, add to no_var
  if (sum(chars != "?") == 1) {
    no_var <- c(no_var, i)
  }
}

# remove the characters with no variation
morpho <- lapply(morpho, function(x) x[-no_var])

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
  
  # if species is extant, la_min and la_max are 0
  if (named_occs$taxon[1] %in% ext_taxa) {
    la_min <- la_max <- 0
  } else {
    # max and min for la
    la_min <- min(la_low$late_age)
    la_max <- max(la_low$early_age)
  }
  
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

# filter data
ftrd_names <- unique(occurrences$taxon)[!(unique(occurrences$taxon) %in% 
                                            ftrd_occs$taxon)]
ftrd_morpho_raw <- morpho_complete[-which(names(morpho_complete) %in% ftrd_names)]
ftrd_mol <- mol_complete[-which(names(mol_complete) %in% ftrd_names)]

# find which characters in morpho are now unvarying

# remove character 46 from morpho, since it doesn't vary
ftrd_morpho <- lapply(ftrd_morpho_raw, function(x) x[-46])

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

###
# plot ranges

# organize it
ranges_plot_df <- data.frame(x1 = as.numeric(ranges_df$fa_max), 
                             x2 = as.numeric(ranges_df$la_max),
                             y = 1:nrow(ranges_df))

# plot
plot(1, type = "n", xlab = "", axes = FALSE,
     ylab = "", xlim = c(40, 0),  
     ylim = c(0, 158)) 
segments(x0 = ranges_plot_df$x1, y0 = ranges_plot_df$y, 
         x1 = ranges_plot_df$x2, y1 = ranges_plot_df$y,
         lwd = 2)
axis_geo(side = 1, intervals = "epochs")
