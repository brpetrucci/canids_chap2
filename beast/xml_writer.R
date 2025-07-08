###########################################
## XML writing script for SRFBD and FBDS ##
## Chapter 2 - Bruno do Rosario Petrucci ##
###########################################

###
# packages

# ape
library(ape)

###
# load in some helpful files

# directory
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/"

# previous scripts (to use when we just need a set, static block)
srfbd_script <- readLines(paste0(base_dir, "srfbd/new_scripts/smorpho_new_1_canidae_run.xml"))
fbds_script <- readLines(paste0(base_dir, "fbds/new_scripts/smorpho_1_both.xml"))

###
# little auxiliary functions

# add spaces easily
add_spaces <- function(n) {
  paste(rep(' ', n), collapse = '')
}

###
# functions to write scripts

# break up data into first and last occurrences
break_ranges <- function(ranges, mol, morpho, attach) {
  # variables to hold results
  new_ranges <- data.frame(matrix(nrow = 0, ncol = 4))
  new_mol <- new_morpho <- list()
  
  # make ages vector
  ages <- c()
  
  # empty mol and morpho vectors for non-attached data
  mol_empty <- paste(rep("?", length(mol[[1]])))
  morpho_empty <- paste(rep("?", length(morpho[[1]])))
  
  # iterate through ranges
  for (i in 1:nrow(ranges)) {
    # this range
    range <- ranges[i, ]
    
    # get taxon
    tx <- range$taxon
    
    # add first specimen to new_ranges
    new_ranges <- rbind(new_ranges, c(paste0(tx, "_first"), tx,
                                      range$fa_max, range$fa_min))
    
    # draw age and add it to vector
    ages <- c(ages, runif(1, range$fa_min, range$fa_max))
    
    # add mol and morpho
    new_mol[[paste0(tx, "_first")]] <- mol_empty
    new_morpho[[paste0(tx, "_first")]] <- morpho[[tx]]
    
    # if it is not a singleton
    if (!(range$fa_max == range$la_max && range$fa_min == range$la_min)) {
      # add last specimen to new_ranges
      new_ranges <- rbind(new_ranges, c(paste0(tx, "_last"), tx,
                                        range$la_max, range$la_min))
      
      # draw age and add it to vector
      ages <- c(ages, runif(1, range$la_min, range$la_max))
      
      # add mol and morpho
      new_mol[[paste0(tx, "_last")]] <- mol[[tx]]
      new_morpho[[paste0(tx, "_last")]] <- ifelse(rep(attach == "both", 
                                                      length(morpho[[1]])),
                                                  morpho[[tx]],
                                                  morpho_empty)
    }
  }
  
  # name everything
  colnames(new_ranges) <- c("specimen", "taxon", "max_age", "min_age")
  names(new_mol) <- names(new_morpho) <- names(ages) <- new_ranges$specimen
  
  # return
  return(list(RANGES = new_ranges, AGES = ages, 
              MOL = new_mol, MORPHO = new_morpho))
}

# write data to script
write_data <- function(data, type) {
  # start data block
  data_block <- paste0(add_spaces(4), '<data')
  
  # add id and datatype if necessary
  data_block <- c(data_block,
                  paste0('id="', type, '"'),
                  'spec="Alignment"')
  if (type == "morpho") data_block <- c(data_block,
                                        'dataType="standard"')
  data_block[length(data_block)] <- paste0(data_block[length(data_block)],
                                           ">")
  # get number of states
  n_states <- sum(!unique(unlist(data)) %in% strsplit("-?nyrswkm", "")[[1]])
  
  # iterate through data
  for (i in 1:length(data)) {
    # start line
    data_line <- paste0(add_spaces(8),
                        '<sequence id="',
                        type,
                        '_seq_',
                        names(data)[i],
                        '" spec="Sequence" taxon="',
                        names(data)[i],
                        '" totalcount="',
                        n_states,
                        '" value="',
                        paste(c(data[[i]], 0:(n_states - 1)), collapse = ""),
                        '"/>')
    
    # add to data block
    data_block <- c(data_block, data_line)
  }
  
  # return data block
  return(c(data_block, ""))
}

# write one whole script
write_script <- function(model, template, maps,
                         ranges, mol, morpho,
                         mol_partitions,
                         attach, gens, store_every, log_every) {
  # variable to hold full script, started with just first line of template
  script <- c(template[1], "")
  
  # break ranges and data up into first and last occurrences
  broken_ranges <- break_ranges(ranges, mol, morpho, attach)
  
  # separate that into each object
  new_ranges <- broken_ranges$RANGES
  ages <- broken_ranges$AGES
  new_mol <- broken_ranges$MOL
  new_morpho <- broken_ranges$MORPHO
  
  # get list of taxa from ranges
  taxa <- unique(new_ranges$taxon)
  
  # add molecular and morphological data
  script <- c(script, 
              write_data(new_mol, "mol"),
              write_data(new_morpho, "morpho"))
  
  # if model is SRFBD, add taxonset certain
  if (model == "SRFBD") {
    # add taxonset declaration
    script <- c(script,
                paste0(add_spaces(4),
                       '<taxonset id="TaxonSet.certain" spec="TaxonSet">'))
    
    # iterate through new ranges
    for (i in 1:nrow(new_ranges)) {
      # add taxon to taxonset
      script <- c(script, paste0(add_spaces(4),
                                 '<taxon id="',
                                 new_ranges$specimen[i],
                                 '" spec="Taxon"/>'))
    }
    
    # finish taxonset
    script <- c(script, paste0(add_spaces(4),
                               '</taxonset>'), "")
  }
  
  # add maps
  script <- c(script, maps, "")
  
  # start MCMC block
  script <- c(script,
              paste0(add_spaces(4),
                     '<run id="mcmc" spec="MCMC" chainLength="',
                     gens, 
                     '" numInitializationAttempts="100000">'),
              paste0(add_spaces(8),
                     '<state id="state" spec="State" storeEvery="',
                     store_every, '">'),
              paste0(add_spaces(12),
                     '<tree id="Tree.t:tree" spec="',
                     ifelse(model == "SRFBD",
                            'sr.evolution.tree.SRTree" nodetype="sr.evolution.tree.SRNode" ',
                            'beast.base.evolution.tree.TraitSet" '),
                     'name="stateNode">'),
              paste0(add_spaces(16),
                     '<trait id="dateTrait.t:tree" ',
                     'spec="beast.base.evolution.tree.TraitSet" ',
                     'traitname="date-backward" value="'))
  
  # iterate through ages
  for (i in 1:length(ages)) {
    # add age
    script <- c(script,
                paste0(add_spaces(16),
                       names(ages)[i],
                       ' = ', ages[i],
                       ifelse(i == length(ages), '">', ',')))
  }
  
  # add taxonset for first DNA partition
  script <- c(script,
              paste0(add_spaces(20),
                     '<taxa id="TaxonSet.',
                     names(mol_partitions)[1],
                     '" spec="TaxonSet">'),
              paste0(add_spaces(24),
                     '<alignment id="',
                     names(mol_partitions)[1],
                     '" spec="FilteredAlignment" filter="',
                     mol_partitions[1], '">'),
              paste0(add_spaces(28),
                     '<data idref="mol"/>'),
              paste0(add_spaces(24), '</alignment>'),
              paste0(add_spaces(20), '</taxa>'),
              paste0(add_spaces(16), '</trait>'),
              paste0(add_spaces(16), 
                     '<taxonset idref="TaxonSet.',
                     names(mol_partitions)[1], '">'))
  
  # if model is SRFBD, need to add sranges to tree
  if (model == "SRFBD") {
    # iterate through taxa
    for (i in 1:length(taxa)) {
      # taxon
      tx <- taxa[i]
      
      # first and last occurrences
      first_occ <- paste0('@', tx, '_first')
      last_occ <- ifelse(any(grepl(paste0(tx, '_last'), new_ranges$specimen)),
                         paste0('@', tx, '_last'),
                         first_occ)
      
      # add to script
      script <- c(script,
                  paste0(add_spaces(16),
                         '<stratigraphicRange id="r', i,
                         '" spec="StratigraphicRange" ',
                         'firstOccurrence="',
                         first_occ,
                         '" lastOccurrence="',
                         last_occ, '"/>'))
    }
  }
  
  #### parameter list--and here things get HARD
  #### init, with stratigraphic ranges if SRFBD, and feast statements for both
  #### priors
    #### fossil ages, if FBDS
  #### likelihood
    #### alignment filters, branch and site models etc
  #### operators
    #### tip dates sampler if FBDS
    #### sample node date random walker if SRFBD
  #### loggers
  
  
}

# write all scripts
write_n_scripts <- function(model, template, n_scripts)

###
# write scripts

# data
ranges <- read.delim(paste0(base_dir, "srfbd/data/ranges.tsv"))
mol <- read.nexus.data(paste0(base_dir, "srfbd/data/mol.nex"))
morpho <- read.nexus.data(paste0(base_dir, "srfbd/data/morpho.nex"))

# partitions
mol_partitions <- c(nDNA = "1-14742",
                    mt12DNA = "14743-15426\3,15427-16570\3,14744-15426\3,15428-16570\3",
                    mt3DNA = "14745-15426\3,15429-16570\3")

template <- srfbd_script
maps <- readLines(paste0(base_dir, "maps.xml"))
