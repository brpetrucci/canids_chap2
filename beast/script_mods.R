scripts_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/srfbd/scripts/"
fbds_scripts_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/fbds/scripts/"

change_morpho_data <- function(script_files) {
  for (j in 1:length(script_files)) {
    script <- readLines(paste0(script_files[j]))
    
    idx_morpho <- min(grep("<sequence", script))
    idx_morpho_end <- min(grep("userDataType", script)) - 1
    
    idx_morpho_sub <- grep("seq_Vulpes_ferrilata", script)
    idx_morpho_sub <- idx_morpho_sub[max(which(idx_morpho_sub <= idx_morpho_end))]
    morpho_sub <- sub('.*value="([^"]+)".*', '\\1', script[idx_morpho_sub])
    
    forl <- c("first", "last")
    
    for (fl in forl) {
      script_new <- script
      
      idx_subs <- idx_morpho - 1 + 
        grep(paste0("_", fl), script_new[idx_morpho:idx_morpho_end])
      
      for (i in 1:length(idx_subs)) {
        if (fl == "first") {
          tx <- sub('.*taxon="([^"]+)_first".*', '\\1', script_new[idx_subs[i]])
          
          if (!any(grepl(paste0('taxon="', tx, '_last"'), script_new[idx_morpho:idx_morpho_end]))) {
            next 
          }
        }
        
        newLine <- gsub('value="([^"]+)"', 
                        paste0('value="', morpho_sub, '"'), 
                        script_new[idx_subs[i]])
        script_new[idx_subs[i]] <- newLine

      }
      
      writeLines(script_new, gsub('canidae_run', paste0(forl[forl != fl], 
                                                        '_canidae_run'),
                                  script_files[j]))
    }
  }
 
}

remove_extra_ranges <- function(script_files, out_files, species_all, 
                                ftrd_species) {
  for (j in 1:length(script_files)) {
    script <- readLines(paste0(script_files[j]))
    
    if (grepl("ftrd", script_files[j])) species <- ftrd_species else species <- species_all
    
    for (i in 1:length(species)) {
      sp <- species[i]
      
      idx_seqs <- grep(paste0('taxon="', sp, "_first", '"'), script)
      
      idx_taxon <- grep(paste0('id="', sp, "_first", '"'), script)
      
      idx_sampling <- grep(paste0('taxon="@', sp, "_first"), script)
      
      idx_age <- grep(paste0(sp, "_first", " = ([0-9.]+),"), script)
      if (length(idx_age) == 0) {
        idx_age <- grep(paste0(sp, "_first", " = ([0-9.]+)"), script)
        
        idx_last_comma <- idx_age - 1
        
        script[idx_last_comma] <- sub(',', '', script[idx_last_comma])
      }
      
      script <- script[-c(idx_seqs, idx_taxon, idx_age, idx_sampling)]
      

      idx_srange <- grep(paste0('firstOccurrence="@', sp, "_first", '"'), script)
        
      script[idx_srange] <- sub("_first", "_last", script[idx_srange])
    
      script[grep(sp, script)] <- gsub("_last", "_first", script[grep(sp, script)])
    }
    
    idx_chainlen <- grep("mcmc", script)
    script[idx_chainlen] <- sub('chainLength="([0-9]+)"', 'chainLength="2000000000"', 
                                script[idx_chainlen])
    
    writeLines(script, out_files[j])
  }
  
}


ranges <- read.delim(paste0(scripts_dir, "../data/ranges.tsv"), sep = "\t")
ftrd_ranges <- read.delim(paste0(scripts_dir, "../data/ftrd_ranges.tsv"), sep = "\t")
species_all <- ranges$taxon[ranges$fa_max == ranges$la_max &
                          ranges$fa_min == ranges$la_min]
ftrd_species <- ftrd_ranges$taxon[ftrd_ranges$fa_max == ftrd_ranges$la_max &
                               ftrd_ranges$fa_min == ftrd_ranges$la_min]
script_files <- c(paste0(scripts_dir, 1:4, "_canidae_run.xml"),
                  paste0(scripts_dir, 1:4, "_ftrd_canidae_run.xml"))
out_files <- gsub("/scripts/", "/new_scripts/new_", script_files)
remove_extra_ranges(script_files, out_files, species_all, ftrd_species)

script_files <- out_files
change_morpho_data(script_files)


remove_dna <- function(script_files) {
  for (i in 1:length(script_files)) {
    script <- readLines(script_files[i])
    
    idx_mol <- min(grep("<data", script))
    idx_end_mol <- min(grep("</data", script))
    
    idx_sequences <- grep("<sequence", script)
    idx_sequences <- idx_sequences[idx_sequences < idx_end_mol]
    
    mol_len <- length(strsplit(sub('.*value="([^ ]+)".*', '\\1',
                                   script[idx_sequences])[1], "")[[1]])
    mol_sub <- paste(rep("?", mol_len), collapse = "")
    
    for (j in 1:length(idx_sequences)) {
      taxon <- sub('.*taxon="([^ ]+)".*', '\\1', script[idx_sequences[j]])
      
      if (!grepl("_ext", taxon)) {
        newLine <- gsub('value="([^"]+)"', 
                        paste0('value="', mol_sub, '"'), 
                        script[idx_sequences[j]])
        script[idx_sequences[j]] <- newLine
      }
    }
    
    idx_chainlen <- grep("mcmc", script)
    script[idx_chainlen] <- sub('chainLength="([0-9]+)"', 'chainLength="2000000000"', 
                                script[idx_chainlen])
    
    writeLines(script, script_files[i])
  }
}

script_files <- paste0(fbds_scripts_dir,
                       paste(sort(paste(rep(1:4, 6), 
                                        c(rep(c("_first", "_last", "_both"), 4),
                  rep(c("_ftrd_first", "_ftrd_last", "_ftrd_both"), 4)), sep = "")),
                  ".xml", sep = ""))
remove_dna(script_files)

add_dna <- function(script_files, out_files) {
  for (i in 1:length(script_files)) {
    script <- readLines(script_files[i])
    
    idx_mol <- min(grep("<data", script))
    idx_end_mol <- min(grep("</data", script))
    
    idx_sequences <- grep("<sequence", script)
    idx_sequences <- idx_sequences[idx_sequences < idx_end_mol]
    
    for (j in 1:length(idx_sequences)) {
      specimen <- sub('.*taxon="([^ ]+)".*', '\\1', script[idx_sequences[j]])
      taxon <- sub('_last|_first|_ext', '', specimen)
      
      if (grepl("_last|_first", specimen) &&
          any(grepl(paste0(taxon, "_ext"), script[idx_sequences]))) {
        idx_dna <- grep(paste0('taxon="', taxon, "_ext"), script[idx_sequences])
        
        dna <- sub('.*value="([^ ]+)".*', '\\1', script[idx_sequences[idx_dna]])
        
        newLine <- gsub('value="([^"]+)"', 
                        paste0('value="', dna, '"'),
                        script[idx_sequences[j]])
        script[idx_sequences[j]] <- newLine
      }
    }
    
    idx_chainlen <- grep("mcmc", script)
    script[idx_chainlen] <- sub('chainLength="([0-9]+)"', 'chainLength="2000000000"', 
                                script[idx_chainlen])
    
    writeLines(script, out_files[i])
  }
}

script_files <- c(script_files,
                  paste0(fbds_scripts_dir, 
                         paste(c(rep("smorpho_", 8), rep("stm_", 8)),
                               paste(rep(1:4, 2),
                                     rep(c(rep("_both.xml", 4), 
                                           rep("_first.xml", 4))), 
                                     sep = ""), sep = "")))
script_files <- script_files[!grepl("_last|_ftrd", script_files)]
out_files <- sub("/scripts/", "/new_scripts/", script_files)
add_dna(script_files, out_files)

mol <- read.nexus.data(paste0(fbds_scripts_dir, "../data/mol_both.nex"))
for (i in 1:length(script_files)) {
  script <- readLines(script_files[i]) 
  
  idx_mol <- min(grep("<data", script))
  idx_end_mol <- min(grep("</data", script))
  
  idx_dusi <- ifelse(grepl("first", script_files[i]),
                     min(grep('taxon="Dusicyon_australis_first"', script)),
                     min(grep('taxon="Dusicyon_australis_last"', script)))
  idx_aeno <- ifelse(grepl("first", script_files[i]),
                     min(grep('taxon="Aenocyon_dirus_first"', script)),
                     min(grep('taxon="Aenocyon_dirus_first"', script)))
  
  dna_dusi <- paste(mol$Dusicyon_australis_last, collapse = "")
  dna_aeno <- paste(mol$Aenocyon_dirus_last, collapse = "")
  
  newLine_dusi <- gsub('value="([^"]+)"', 
                       paste0('value="', dna_dusi, '"'),
                       script[idx_dusi])
  newLine_aeno <- gsub('value="([^"]+)"', 
                       paste0('value="', dna_aeno, '"'),
                       script[idx_aeno])
  
  script[idx_dusi] <- newLine_dusi
  script[idx_aeno] <- newLine_aeno
  

  for (j in 1:length(extant_species)) {
    sp <- extant_species[j]
    
    idx_seqs <- grep(paste0('taxon="', sp, "_last"), script)
    
    if (length(idx_seqs) == 0) next
    
    idx_trait <- grep(paste0(sp, '_last=([0-9.]+),'), script)
    script[idx_trait] <- sub(paste0(sp, '_last=([0-9.]+),'), '', script[idx_trait])
    
    idx_prior <- grep(paste0('id="', sp, "_lastSet.prior"), script)
    idx_prior_end <- idx_prior + 5
    
    idx_operator <- grep(paste0('id="tipDatesSampler.', sp, '_lastSet"'), script)
    
    idx_log <- grep(paste0('idref="', sp, '_lastSet.prior'), script)
    
    script <- script[-c(idx_seqs, idx_prior:idx_prior_end, idx_operator, idx_log)]
  }
  
  idx_chainlen <- grep("mcmc", script)
  script[idx_chainlen] <- sub('chainLength="([0-9]+)"', 'chainLength="2000000000"', 
                              script[idx_chainlen])
  
  writeLines(script, out_files[i])
}
