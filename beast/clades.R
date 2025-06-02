###
# set up

# ape
library(ape)

# set directories
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/"
srfbd_dir <- paste0(base_dir, "srfbd/")

###
# define clades

# get all tips from the most complete tree
all_tips <- read.nexus(paste0(srfbd_dir, "output/1_canidae_run.trees"))[[1]]$tip.label

# species (no _first or _last)
all_species <- unique(sub("_first", "", sub("_last", "", all_tips)))

# singletons
singletons <- all_species[!(paste0(all_species, "_last") %in% all_tips)]

# get all genera
all_genera_names <- unique(sub("_([a-zA-Z]+)", "", all_species))
genera <- lapply(all_genera_names, function(x) {
  all_species[grep(x, all_species)]
})
names(genera) <- all_genera_names

# select just the extant genera
extant_genera_names <- c("Urocyon", "Otocyon", "Nyctereutes", "Vulpes",
                         "Speothos", "Chrysocyon", "Cerdocyon", 
                         "Atelocynus", "Lycalopex", "Lycaon", "Cuon",
                         "Lupulella", "Canis")
extant_genera <- lapply(extant_genera_names, function(x) genera[[x]])
names(extant_genera) <- extant_genera_names

# and then species
extant_species <- sort(unname(unlist(extant_genera)))[-c(3:5, 7:10, 12, 14:15,
                                                        17:19, 36:37, 39:40, 49)]
recent_species <- c(extant_species, "Aenocyon_dirus", "Dusicyon_australis")

###
# get all the higher taxonomic levels defined by genera

# caninae - basic groups
caninae_basal_genera_names <- c("Leptocyon", "Urocyon", "Metalopex")
vulpini_genera_names <- c("Otocyon", "Nyctereutes", "Vulpes")
canina_genera_names <- c("Eucyon", "Aenocyon", "Xenocyon", "Lupulella",
                         "Lycaon", "Cuon", "Canis")
cerdocyonina_genera_names <- c("Speothos", "Chrysocyon", "Cerdocyon", 
                               "Atelocynus", "Lycalopex", "Dusicyon")
canini_genera_names <- c(canina_genera_names, cerdocyonina_genera_names)
caninae_genera_names <- c(caninae_basal_genera_names,
                          vulpini_genera_names, canini_genera_names)

# caninae - groups of interest
caninae_nolepto_genera_names <- caninae_genera_names[-1]
urometa_genera_names <- c("Urocyon", "Metalopex")
vulpinicanini_genera_names <- c(canini_genera_names, vulpini_genera_names)
otonyc_genera_names <- c("Otocyon", "Nyctereutes")
canina_noeuc_genera_names <- canina_genera_names[-1]
canina_nolup_genera_names <- canina_noeuc_genera_names[-3]
canina_noaen_genera_names <- canina_noeuc_genera_names[-1]
canis_aenocyon_genera_names <- c("Canis", "Aenocyon")
xlcc_genera_names <- c("Canis", "Xenocyon", "Lycaon", "Cuon")
xlc_genera_names <- xlcc_genera_names[-1]
chryspedusi_genera_names <- c("Chrysocyon", "Speothos", "Dusicyon")
cerdocyonina_crown_genera_names <- c("Atelocynus", "Cerdocyon", "Lycalopex")
atelocerdo_genera_names <- c("Atelocynus", "Cerdocyon")

# borophaginae - basic groups
borophaginae_basal_genera_names <- c("Archaeocyon", "Otarocyon", "Rhizocyon",
                                     "Oxetocyon")
phlaocyonini_genera_names <- c("Cynarctoides", "Phlaocyon")
borophagina_genera_names <- c("Paratomarctus", "Carpocyon", "Protepicyon",
                              "Epicyon", "Borophagus")
aelurodontina_genera_names <- c("Tomarctus", "Aelurodon")
cynarctina_genera_names <- c("Paracynarctus", "Cynarctus")
borophagini_genera_names <- c("Cormocyon", "Desmocyon", "Metatomarctus", 
                              "Euoplocyon", "Psalidocyon", "Microtomarctus",
                              "Protomarctus", "Tephrocyon",
                              borophagina_genera_names, 
                              aelurodontina_genera_names,
                              cynarctina_genera_names)
borophaginae_genera_names <- c(borophaginae_basal_genera_names,
                               phlaocyonini_genera_names,
                               borophagini_genera_names)

# borophaginae - groups of interest
borophaginae_noarchaeo_genera_names <- borophaginae_genera_names[-1]
phlaocyoniniborophagini_genera_names <- borophaginae_noarchaeo_genera_names[-c(1:3)]
borophagini_nobasal_genera_names <- borophagini_genera_names[-c(1, 2)]
borophagini_nocynarc_genera_names <- borophagini_nobasal_genera_names[-c(14:15)]
aeluroboro_genera_names <- c(aelurodontina_genera_names,
                             borophagina_genera_names)

# hesperocyoninae - basic groups
hesperocyoninae_basal_genera_names <- c("Prohesperocyon", "Hesperocyon")
paracaedo_genera_names <- c("Paraenhydrocyon", "Caedocyon")
mse_genera_names <- c("Mesocyon", "Cynodesmus", "Sunkahetanka", 
                      "Philotrox", "Enhydrocyon")
hesperocyoninae_genera_names <- c(hesperocyoninae_basal_genera_names,
                                  paracaedo_genera_names,
                                  mse_genera_names,
                                  "Osbornodon", "Ectopocynus")

# hesperocyoninae - groups of interest
noprohesp_genera_names <- hesperocyoninae_genera_names[-1]
mse_osbornodon_genera_names <- c("Osbornodon", mse_genera_names)

### make clades and conditions list

# remove monotypic genera from clades list
mono_gen <- which(unlist(lapply(genera, function(x) length(x) == 1)))
clades <- genera[-mono_gen]

# function to get a list of genera and return a list of species
species_from_genera <- function(genera) {
  unlist(lapply(genera, function(x) all_species[grep(x, all_species)]))
}

# add caninae basic groups
clades$Canina <- species_from_genera(canina_genera_names)
clades$Cerdocyonina <- species_from_genera(cerdocyonina_genera_names)
clades$Canini <- species_from_genera(canini_genera_names)
clades$Vulpini <- species_from_genera(vulpini_genera_names)
clades$Caninae <- species_from_genera(caninae_genera_names)

# add caninae special groups
clades$NoLeptocyon <- species_from_genera(caninae_nolepto_genera_names)
clades$UrocyonMetalopex <- species_from_genera(urometa_genera_names)
clades$VulpiniCanini <- species_from_genera(vulpinicanini_genera_names)
clades$OtocyonNyctereutes <- species_from_genera(otonyc_genera_names)
clades$NoEucyon <- species_from_genera(canina_noeuc_genera_names)
clades$NoLupulella <- species_from_genera(canina_nolup_genera_names)
clades$NoAenocyon <- species_from_genera(canina_noaen_genera_names)
clades$CanisAenocyon <- species_from_genera(canis_aenocyon_genera_names)
clades$XenocyonLycaonCuonCanis <- species_from_genera(xlcc_genera_names)
clades$XenocyonLycaonCuon <- species_from_genera(xlc_genera_names)
clades$NoTexanus <- sort(clades$Cerdocyonina)[-2]
clades$ChrysocyonDusicyonSpeothos <- species_from_genera(chryspedusi_genera_names)
clades$CerdocyoninaCrown <- species_from_genera(cerdocyonina_crown_genera_names)
clades$AtelocynusCerdocyon <- species_from_genera(atelocerdo_genera_names)

# add conditional clades to conds
conds <- c("Canini", "Canini", "VulpiniCanini", "VulpiniCanini", "Canidae",
           "Caninae", "NoLeptocyon", "NoLeptocyon", "Vulpini", 
           "Canina", "NoEucyon", "NoEucyon", "NoLupulella", 
           "NoLupulella", "XenocyonLycaonCuonCanis", "Cerdocyonina", 
           "NoTexanus", "NoTexanus", "CerdocyoninaCrown")

# add borophaginae basic groups
clades$Phlaocyonini <- species_from_genera(phlaocyonini_genera_names)
clades$Borophagina <- species_from_genera(borophagina_genera_names)
clades$Aelurodontina <- species_from_genera(aelurodontina_genera_names)
clades$Cynarctina <- species_from_genera(cynarctina_genera_names)
clades$Borophagini <- species_from_genera(borophagini_genera_names)
clades$Borophaginae <- species_from_genera(borophaginae_genera_names)

# add borophaginae special groups
clades$NoArchaeocyon <- species_from_genera(borophaginae_noarchaeo_genera_names)
clades$PhlaocyoniniBorophagini <- species_from_genera(phlaocyoniniborophagini_genera_names)
clades$BorophaginiNoBasal <- species_from_genera(borophagini_nobasal_genera_names)
clades$BorophaginiCrown <- species_from_genera(borophagini_nocynarc_genera_names)
clades$AelurodontinaBorophagina <- species_from_genera(aeluroboro_genera_names)

# add conditional clades to conds
conds <- c(conds, "Borophaginae", "Borophagini", "Borophagini", "Borophagini",
           "Borophaginae", "Canidae", "Borophaginae", "NoArchaeocyon", 
           "Borophagini", "BorophaginiNoBasal", "BorophaginiCrown")

# add hesperocyoninae basic groups
clades$ParaCaedo <- species_from_genera(paracaedo_genera_names)
clades$MSE <- species_from_genera(mse_genera_names)
clades$Hesperocyoninae <- species_from_genera(hesperocyoninae_genera_names)

# add hesperocyoninae special groups
clades$NoProhesperocyon <- species_from_genera(noprohesp_genera_names)
clades$NoGregarius <- sort(clades$Hesperocyoninae)[-25]
clades$CanidaeNoBasal <- c(sort(clades$Hesperocyoninae)[-c(12, 25)],
                              clades$Borophaginae, clades$Caninae)
clades$BorophaginaeCaninae <- c(clades$Borophaginae, clades$Caninae)
clades$ColoradensisMSEOsbornodon <- c("Hesperocyon_coloradensis",
                                         species_from_genera(mse_genera_names),
                                         species_from_genera(c("Osbornodon")))
clades$MSEOsbornodon <- clades$ColoradensisMSEOsbornodon[-1]

# add conditional clades for hespero
conds <- c(conds, "Canidae", "Canidae", "Canidae",
           "Canidae", "Canidae", "Canidae", "CanidaeNoBasal", 
           "CanidaeNoBasal", "ColoradensisMSEOsbornodon")

# add conds for non-monotypic genera
conds <- c("Canidae", "MSE", "MSE", "MSE", "MSEOsbornodon", "CanidaeNoBasal", 
           "CanidaeNoBasal", "Borophaginae", "NoArchaeocyon", "Phlaocyonini",
           "Phlaocyonini", "Borophagini", "Borophagini", "Cynarctina", 
           "Cynarctina", "BorophaginiNoBasal", "Aelurodontina", 
           "Aelurodontina", "Borophagina", "Borophagina", "Borophagina",
           "Borophagina", "Caninae", "Vulpini", "NoLeptocyon",
           "NoLeptocyon", "Canina", "NoEucyon", "NoEucyon", "CerdocyoninaCrown",
           "CerdocyoninaCrown", conds)

# some sanity checks
length(conds) == length(clades)
unique(conds)[!(unique(conds) %in% names(clades))]

# add Canidae to both
clades$Canidae <- c(clades$Hesperocyoninae, clades$Borophaginae, clades$Caninae)
conds <- c(conds, "")

# add names to conds
names(conds) <- names(clades)
