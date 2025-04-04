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
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/"

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
get_range <- function(occurrences, name, ext_as_zero = TRUE) {
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
  la_max <- max(la_low$early_age)
  
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

# final ranges (change extants to 0 0)
ranges_final <- ranges_df

# iterate through extant taxa
ranges_final$la_max[ranges_final$taxon %in% ext_taxa] <- 0
ranges_final$la_min[ranges_final$taxon %in% ext_taxa] <- 0

# write ranges to a file
write.table(ranges_final, paste0(base_dir, "srfbd/ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# add ?? for molecular data for species without molecular data
nomol_taxa <- taxa_names[!(taxa_names %in% names(mol))]
nomol_mol <- lapply(nomol_taxa, function(x) rep("?", length(mol[[1]])))
names(nomol_mol) <- nomol_taxa
mol_complete <- c(mol, nomol_mol)

# write DNA to a file
write.nexus.data(mol_complete, paste0(base_dir, "srfbd/mol.nex"))

# reading partitions and writing to DNA file
partitions <- readLines(paste0(base_dir, "partitions_simplified.txt"))
write(partitions, paste0(base_dir, "srfbd/mol.nex"), append = TRUE)

# add ?? for morphological data for species without morphological data
nomorpho_taxa <- taxa_names[!(taxa_names %in% names(morpho))]
nomorpho_morpho <- lapply(nomorpho_taxa, function(x) rep("?", length(morpho[[1]])))
names(nomorpho_morpho) <- nomorpho_taxa
morpho_complete <- c(morpho, nomorpho_morpho)

# write DNA to a file
write.nexus.data(morpho_complete, paste0(base_dir, "srfbd/morpho.nex"),
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

# final ranges (change extants to 0 0)
ftrd_ranges_final <- ftrd_ranges_df

# iterate through extant taxa
ftrd_ranges_final$la_max[ftrd_ranges_final$taxon %in% ext_taxa] <- 0
ftrd_ranges_final$la_min[ftrd_ranges_final$taxon %in% ext_taxa] <- 0

# write ranges and morpho data
write.table(ftrd_ranges_final, paste0(base_dir, "srfbd/ftrd_ranges.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.nexus.data(ftrd_morpho, paste0(base_dir, "srfbd/ftrd_morpho.nex"),
                 format = "standard")
write.nexus.data(ftrd_mol, paste0(base_dir, "srfbd/ftrd_mol.nex"))
write(partitions, paste0(base_dir, "srfbd/ftrd_mol.nex"), append = TRUE)

###
# make specimen level data

# make a function for it
specimen_data <- function(ranges_df, mol, morpho, prefix = "") {
  # create data frame for occurrences
  fbds_taxa_first <- fbds_taxa_last <- fbds_taxa_both <- 
    data.frame(matrix(nrow = 0, ncol = 3))
  
  # and for taxa ages
  fbds_first_ages <- fbds_last_ages <- fbds_both_ages <-
    data.frame(matrix(nrow = 0, ncol = 2))
  
  # create molecular and morpho lists
  mol_first <- mol_last <- mol_both <- morpho_first <- morpho_last <-
    morpho_both <- list()
  
  # iterate through ranges
  for (i in 1:nrow(ranges_df)) {
    # extract this row
    range <- ranges_df[i, -6]
    
    # check if the species is extant
    if (range$taxon %in% ext_taxa) {
      # add extant occurrence to occs
      ext_occ <- c(paste0(range$taxon, "_ext"), 0, 0)
      
      # add ext_occ to the taxa data frames
      fbds_taxa_first <- rbind(fbds_taxa_first, ext_occ)
      fbds_taxa_last <- rbind(fbds_taxa_last, ext_occ)
      fbds_taxa_both <- rbind(fbds_taxa_both, ext_occ)
      
      # add 0 to ages
      fbds_first_ages <- rbind(fbds_first_ages, ext_occ[1:2])
      fbds_last_ages <- rbind(fbds_last_ages, ext_occ[1:2])
      fbds_both_ages <- rbind(fbds_both_ages, ext_occ[1:2])
      
      # add data to mol
      mol_first <- c(mol_first, list(mol[[range$taxon]]))
      mol_last <- c(mol_last, list(mol[[range$taxon]]))
      mol_both <- c(mol_both, list(mol[[range$taxon]]))
      
      # and to morpho
      morpho_first <- c(morpho_first, list(morpho[[range$taxon]]))
      morpho_last <- c(morpho_last, list(morpho[[range$taxon]]))
      morpho_both <- c(morpho_both, list(morpho[[range$taxon]]))
      
      # if it is a singleton, just skip to the next species
      if (all(range[2:length(range)] == 0)) next
    }
    
    # get first and last occurrences
    fa <- c(paste0(range$taxon, "_first"), range$fa_max, range$fa_min)
    la <- c(paste0(range$taxon, "_last"), range$la_max, range$la_min)
    
    # add fa to first and both
    fbds_taxa_first <- rbind(fbds_taxa_first, fa)
    fbds_taxa_both <- rbind(fbds_taxa_both, fa)
    
    # get an age from a uniform draw
    fa_age <- runif(1, as.numeric(fa[3]), as.numeric(fa[2]))
    
    # add age to age data frames
    fbds_first_ages <- rbind(fbds_first_ages, c(fa[1], fa_age))
    fbds_both_ages <- rbind(fbds_both_ages, c(fa[1], fa_age))
    
    # add data to mol
    mol_first <- c(mol_first, list(mol[[range$taxon]]))
    mol_both <- c(mol_both, list(mol[[range$taxon]]))
    
    # and to morpho
    morpho_first <- c(morpho_first, list(morpho[[range$taxon]]))
    morpho_both <- c(morpho_both, list(morpho[[range$taxon]]))
    
    # if it is a singleton, add to last as well
    if (all(fa[-1] == la[-1])) {
      fbds_taxa_last <- rbind(fbds_taxa_last, fa)
      fbds_last_ages <- rbind(fbds_last_ages, c(fa[1], fa_age))
      
      # and add to mol and morpho
      mol_last <- c(mol_last, list(mol[[range$taxon]]))
      morpho_last <- c(morpho_last, list(morpho[[range$taxon]]))
    } else {
      # if not add la to last and both
      fbds_taxa_last <- rbind(fbds_taxa_last, la)
      fbds_taxa_both <- rbind(fbds_taxa_both, la)
      
      # get an age from a uniform draw
      la_age <- runif(1, as.numeric(la[3]), as.numeric(la[2]))
      
      # add age to age data frames
      fbds_last_ages <- rbind(fbds_last_ages, c(la[1], la_age))
      fbds_both_ages <- rbind(fbds_both_ages, c(la[1], la_age))
      
      # add data to mol
      mol_last <- c(mol_last, list(mol[[range$taxon]]))
      mol_both <- c(mol_both, list(mol[[range$taxon]]))
      
      # and to morpho
      morpho_last <- c(morpho_last, list(morpho[[range$taxon]]))
      morpho_both <- c(morpho_both, list(morpho[[range$taxon]]))
    }
  }
  
  # name occurrences data frames
  colnames(fbds_taxa_first) <- colnames(fbds_taxa_last) <- 
    colnames(fbds_taxa_both) <- c("taxon", "max_age", "min_age")
  
  # name ages data frames
  colnames(fbds_first_ages) <- colnames(fbds_last_ages) <-
    colnames(fbds_both_ages) <- c("taxon", "age")
  
  # and list of molecular and morpho data
  names(mol_first) <- names(morpho_first) <- fbds_taxa_first$taxon
  names(mol_last) <- names(morpho_last) <- fbds_taxa_last$taxon
  names(mol_both) <- names(morpho_both) <- fbds_taxa_both$taxon
  
  # write ranges to files
  write.table(fbds_taxa_first, paste0(base_dir, "fbds/data/", prefix,
                                      "taxa_first.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fbds_taxa_last, paste0(base_dir, "fbds/data/", prefix,
                                     "taxa_last.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(fbds_taxa_both, paste0(base_dir, "fbds/data/", prefix,
                                     "taxa_both.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # write ages to files
  write.table(fbds_first_ages, paste0(base_dir, "fbds/data/", prefix,
                                      "ages_first.dat"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(fbds_last_ages, paste0(base_dir, "fbds/data/", prefix,
                                     "ages_last.dat"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(fbds_both_ages, paste0(base_dir, "fbds/data/", prefix,
                                     "ages_both.dat"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # write all molecular data to files
  write.nexus.data(mol_first, paste0(base_dir, "fbds/data/", prefix,
                                     "mol_first.nex"))
  write.nexus.data(mol_last, paste0(base_dir, "fbds/data/", prefix,
                                    "mol_last.nex"))
  write.nexus.data(mol_both, paste0(base_dir, "fbds/data/", prefix,
                                    "mol_both.nex"))
  
  # append partitions to it
  write(partitions, paste0(base_dir, "fbds/data/", prefix,
                           "mol_first.nex"), append = TRUE)
  write(partitions, paste0(base_dir, "fbds/data/", prefix,
                           "mol_last.nex"), append = TRUE)
  write(partitions, paste0(base_dir, "fbds/data/", prefix,
                           "mol_both.nex"), append = TRUE)
  
  # write all morpho data to files
  write.nexus.data(morpho_first, paste0(base_dir, "fbds/data/", prefix,
                                        "morpho_first.nex"),
                   format = "standard")
  write.nexus.data(morpho_last, paste0(base_dir, "fbds/data/", prefix,
                                       "morpho_last.nex"),
                   format = "standard")
  write.nexus.data(morpho_both, paste0(base_dir, "fbds/data/", prefix,
                                       "morpho_both.nex"),
                   format = "standard")
  
  return(list(FRANGES = fbds_taxa_first, LRANGES = fbds_taxa_last,
              BRANGES = fbds_taxa_both, FAGES = fbds_first_ages,
              LAGES = fbds_last_ages, BAGES = fbds_both_ages,
              FMOL = mol_first, LMOL = mol_last, 
              BMOL = mol_both, FMORPHO = morpho_first, 
              LMORPHO = morpho_last, BMORPHO = morpho_both))
}

# get specimen data for full dataset
full_fbds_data <- specimen_data(ranges_df, mol_complete, morpho_complete)

# and for ftrd dataset
ftrd_fbds_data <- specimen_data(ftrd_ranges_df, ftrd_mol, ftrd_morpho,
                                prefix = "ftrd_")

###
# saving extant data

# get extant mol and morpho data
mol_extant <- mol_complete[ext_taxa]
morpho_extant <- morpho_complete[ext_taxa]

# find which morphological characters are now constant
morpho_constant <- which(unlist(lapply(1:length(morpho_extant[[1]]), 
              function(x) length(unique(unlist(lapply(1:length(morpho_extant), 
            function(y) morpho_extant[[y]][x])))))) < 3)

# select only the characters that vary
morpho_extant_var <- lapply(morpho_extant, function(x) x[-morpho_constant])

# write molecular data
write.nexus.data(mol_extant, paste0(base_dir, "extant/mol.nex"), 
                 format = "dna")

# write partitions
write(partitions, paste0(base_dir, "extant/mol.nex"), append = TRUE)

write.nexus.data(morpho_extant_var, paste0(base_dir, "extant/morpho.nex"),
                 format = "standard")

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
