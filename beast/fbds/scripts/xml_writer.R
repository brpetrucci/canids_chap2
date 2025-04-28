#######################################
##            Chapter 2              ##
##      Canids complete tree         ##
##        Edit xml scripts           ##
##     Bruno do Rosario Petrucci     ##
#######################################

###
# load libraries

# stringr
library(stringr)

###
# read scripts and data

# base directory
fbds_base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/fbds/"
extant_base_dir <- "/Users/petrucci/Documents/research/canids_chap2/beast/extant/"

# data and scripts directory
fbds_data_dir <- paste0(fbds_base_dir, "data/")
fbds_scripts_dir <- paste0(fbds_base_dir, "scripts/")
extant_data_dir <- paste0(extant_base_dir, "data/")
extant_scripts_dir <- paste0(extant_base_dir, "scripts/")

###
# auxiliary functions

# function to add spaces easily
add_spaces <- function(n) {
  paste(rep(' ', n), collapse = '')
}

# drawing from a dirichlet distribution
rdirichlet <- function (n, alpha) {
  # length of parameter vector
  l <- length(alpha)
  
  # make the draws from the gamma
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  
  # normalize the draws
  sm <- x %*% rep(1, l)
  
  # return normalized vector
  x / as.vector(sm)
}

###
# write function for reading and editing the xml
xml_editor <- function(template_file, out_dir, out_file, taxa_file = NULL, 
                       fossils = FALSE, n_scripts = 4) {
  # read template script
  temp <- readLines(template_file)
  
  # start xml as the same as template
  xml <- temp
  
  # if there are fossils on this datasets
  if (fossils) {
    # get index to write fossil age uncertainty distributions
    idx_prior <- max(grep("</prior>", temp))
    
    # index for last operator
    idx_operator <- max(grep("</operator>", temp))
    
    # index for idref line (2 lines before screenlog)
    idx_idref <- grep("screenlog", temp) - 2
    
    # get taxa file
    taxa <- read.delim(taxa_file)
    
    # taxa without extant singletons
    taxa_fossils <- taxa[taxa$min_age != taxa$max_age, ]
    
    # iterate through each taxon
    for (i in 1:nrow(taxa_fossils)) {
      # get this taxon
      taxon <- taxa_fossils[i, ]
      
      # start prior lines
      prior_lines <- vector("character", 6)
      
      # add distribution id line
      prior_lines[1] <- paste0(add_spaces(16), '<distribution id="',
                               taxon$taxon,
                               'Set.prior"',
                               ' spec="sa.math.distributions.SAMRCAPrior"',
                               ' tipsonly="true" tree="@Tree.t:tree">')
      
      # add taxon set lines
      prior_lines[2] <- paste0(add_spaces(20), '<taxonset id="',
                               taxon$taxon,
                               'Set" spec="TaxonSet">')
      prior_lines[3] <- paste0(add_spaces(24), '<taxon id="',
                               taxon$taxon,
                               '" spec="Taxon"/>')
      prior_lines[4] <- paste0(add_spaces(20), '</taxonset>')
      
      # add uniform distribution
      prior_lines[5] <- paste0(add_spaces(20), '<Uniform id="FossilAges.', i,
                               '" lower="', taxon$min_age,
                               '" name="distr" upper="', taxon$max_age,
                               '"/>')
      prior_lines[6] <- paste0(add_spaces(16), '</distribution>')
      
      # add to xml
      xml <- c(xml[1:idx_prior], prior_lines, xml[(idx_prior + 1):length(xml)])
      
      # add to prior index
      lines_added <- length(prior_lines)
      idx_prior <- idx_prior + lines_added
      
      # write operator line
      operator <- paste0(add_spaces(8), '<operator id="tipDatesSampler.',
                         taxon$taxon, 'Set" ',
                         'spec="sa.evolution.operators.',
                         'SampledNodeDateRandomWalker" ',
                         'taxonset="@', taxon$taxon, 'Set" ',
                         'tree="@Tree.t:tree" weight="1.0" windowSize="1.0"/>')
      
      # add to xml
      xml <- c(xml[1:(idx_operator + lines_added)], operator,
               xml[(idx_operator + lines_added + 1):length(xml)])
      
      # add to lines_added and index
      lines_added <- lines_added + 1
      idx_operator <- idx_operator + lines_added
      
      # write idref line
      idref <- paste0(add_spaces(12), '<log idref="', taxon$taxon, 'Set.prior"/>')
      
      # add to xml
      xml <- c(xml[1:(idx_idref + lines_added)], idref,
               xml[(idx_idref + lines_added + 1):length(xml)])
      
      # add to lines_added
      lines_added <- lines_added + 1
      idx_idref <- idx_idref + lines_added
    }
    
    # get index where starting ages are supposed to go
    idx_ages <- grep("<tree id=", xml)
  }
  
  # iterate through number of scripts
  for (i in 1:n_scripts) {
    # start final script
    script <- xml
    
    # if there are fossils 
    if (fossils) {
      # start dates lines
      dates <- paste0('<trait id="dateTrait.t:tree" ',
                      'spec="beast.base.evolution.tree.TraitSet" ',
                      'traitname="date-backward" value="')
      
      # iterate through taxa
      for (j in 1:nrow(taxa)) {
        # get taxon
        taxon <- taxa[j, ]
        
        # get age
        age <- runif(1, taxon$min_age, taxon$max_age)
        
        # append to dates
        dates <- paste0(dates, taxon$taxon, '=', age,
                        ifelse(j == nrow(taxa), '">', ','))
      }
      
      # add dates to script
      script <- c(script[1:idx_ages], paste0(add_spaces(16), dates),
                  script[(idx_ages + 1):(idx_ages + 5)],
                  paste0(add_spaces(16), "</trait>"),
                  script[(idx_ages + 6):length(script)])
    }
    
    # get index for </state>, i.e. right after the last parameter declaration
    idx_state <- grep("</state>", script)
    
    # index for all lines starting with <parameter that are before idx_state
    idx_params <- grep("<parameter", script[1:idx_state])
    
    # get all parameter lines
    params <- script[idx_params]

    # iterate through parameters
    for (j in 1:length(params)) {
      # parameter line
      param_line <- params[j]
      
      # get parameter name
      param <- sub('.*id="([^"]+)".*', '\\1', param_line)
      
      # prior on what?
      prior_target <- sub('.*:([^:]+).*', '\\1', param)
      
      # get parameter starting value
      param_start <- as.numeric(sub('.*>([0-9.]+)<.*', '\\1', param_line))
      
      # find whether there is a line with a prior for it
      idx_prior_param <- grep(paste0('x="@', param), script)
      
      # if not, next
      if (length(idx_prior_param) == 0 || param == "originFBD.t:tree") next
      
      # find how many lines after the prior ends
      prior_length <- 
        min(grep("</prior>", script[idx_prior_param:length(script)]))
      
      # get entire prior
      prior_param_lines <- 
        script[idx_prior_param:(idx_prior_param + prior_length - 1)]
      
      # get prior distribution name and id
      prior_dist <- sub('.*<([^ ]+) i.*', '\\1', prior_param_lines[2])
      prior_id <- sub('.*id="([^"]+)".*', '\\1', prior_param_lines[2])
      
      # get parameters of the prior
      prior_hyper <- prior_param_lines[grep("<parameter", prior_param_lines)]
      if (length(prior_hyper) == 0) {
        prior_hyper <- prior_param_lines[grep("<mean", prior_param_lines)]
      }
      
      # attempt to get lower and upper
      param_lower <- sub('.*lower="([0-9.]+)".*', '\\1', param_line)
      if (param_lower == param_line) param_lower <- -Inf
      param_upper <- sub('.*upper="([0-9.]+)".*', '\\1', param_line)
      if (param_upper == param_line) param_upper <- Inf
      param_lower <- as.numeric(param_lower)
      param_upper <- as.numeric(param_upper)
      
      # initialize condition for while loop
      cond <- FALSE
      
      # ensure that new_start is not too far from param_start
      while(!cond) {
        # switch for the case of each prior
        switch(prior_dist,
               "Exponential" = {
                 # extract mean
                 prior_mean <- sub('.*>([0-9.]+)<.*', '\\1', prior_hyper)
                 if (prior_mean == prior_hyper) {
                   prior_mean <- sub('.*value="([0-9.]+).*', '\\1', prior_hyper)
                 }
                 
                 # make it numeric
                 prior_mean <- as.numeric(prior_mean)
                 
                 # get the starting value
                 new_start <- rexp(1, 1/prior_mean)
               },
               "LogNormal" = {
                 # all the lognormals have the same mean and sd,
                 # so we can just draw it
                 new_start <- rlnorm(1, log(1), 0.2)
               },
               "Uniform" = {
                 # just between lower and upper
                 new_start <- runif(1, param_lower, param_upper)
               },
               "Beta" = {
                 # only one beta, so easy
                 new_start <- rbeta(1, 2, 2)
               },
               "Gamma" = {
                 # get the parameters
                 prior_rates <- as.numeric(sub('.*>([0-9.]+)<.*', '\\1', 
                                               prior_hyper))
                 
                 # draw a value
                 new_start <- rgamma(1, shape = prior_rates[1], 
                                     scale = prior_rates[2])
               },
               "distr" = {
                 # again all dirichlet have the same parameters, so 
                 new_start <- rdirichlet(1, c(4, 4, 4, 4))
                 
                 # round frequencies to 6 digits to avoid truncation later
                 new_start <- round(new_start, digits = 6)
                 
                 # ensure it still sums to 1
                 new_start[1] <- new_start[1] + (1 - sum(new_start))
                 # since it's just on the order of 10^-6, 
                 # hopefully won't change much
               },
               stop("Missing a distribution"))
        
        # check condition
        cond <- (mean(new_start) <= param_upper && 
                   mean(new_start) >= param_lower)
      }
      
      # if new_start is a higher length, make it one string
      if (length(new_start) > 1) {
        new_start <- paste(new_start, collapse = " ")
      }
        
      # get index of this particular parameter
      idx_param_sub <- grep(param_line, script)
      
      # change the starting value
      script[idx_param_sub] <- gsub(">([0-9.]+)<", paste0(">", new_start, "<"),
                                    param_line)
      
      # if parameter is from the tree prior, we want to add an init line
      if (prior_target == "tree") {
        # find init line
        idx_init <- max(grep("</init>", script))
        
        # make init lines vector
        init_lines <- vector("character", 3)
        
        # add first line
        init_lines[1] <- paste0(add_spaces(8),
                                '<init spec="feast.parameter.',
                                'RandomRealParameter" initial="@',
                                param, '">')
        
        # add line with distribution id
        init_lines[2] <- paste0(add_spaces(12),
                                '<distr idref="',
                                prior_id,
                                '"/>')
        init_lines[3] <- paste0(add_spaces(8), "</init>")
        
        # add to script
        script <- c(script[1:idx_init],
                    init_lines,
                    script[(idx_init + 1):length(script)])
      }
    }
    
    # save script
    writeLines(script, paste0(out_dir, i, '_', out_file))
  }
}

###
# collect files and run function

# out directory for fbds scripts
out_dir <- fbds_scripts_dir

# last ages
template_file <- paste0(fbds_scripts_dir, "both_template.xml")
taxa_file <- paste0(fbds_data_dir, "taxa_both.tsv")
out_file <- "both.xml"
xml_editor(template_file, out_dir, out_file, taxa_file, fossils = TRUE)

# first ages
template_file <- paste0(fbds_scripts_dir, "first_template.xml")
taxa_file <- paste0(fbds_data_dir, "taxa_first.tsv")
out_file <- "first.xml"
xml_editor(template_file, out_dir, out_file, taxa_file, fossils = TRUE)

# last ages
template_file <- paste0(fbds_scripts_dir, "last_template.xml")
taxa_file <- paste0(fbds_data_dir, "taxa_last.tsv")
out_file <- "last.xml"
xml_editor(template_file, out_dir, out_file, taxa_file, fossils = TRUE)

# out directory for fbds ftrd scripts
out_dir <- fbds_scripts_dir

# last ages
template_file <- paste0(fbds_scripts_dir, "ftrd_both_template.xml")
taxa_file <- paste0(fbds_data_dir, "ftrd_taxa_both.tsv")
out_file <- "ftrd_both.xml"
xml_editor(template_file, out_dir, out_file, taxa_file, fossils = TRUE)

# first ages
template_file <- paste0(fbds_scripts_dir, "ftrd_first_template.xml")
taxa_file <- paste0(fbds_data_dir, "ftrd_taxa_first.tsv")
out_file <- "ftrd_first.xml"
xml_editor(template_file, out_dir, out_file, taxa_file, fossils = TRUE)

# last ages
template_file <- paste0(fbds_scripts_dir, "ftrd_last_template.xml")
taxa_file <- paste0(fbds_data_dir, "ftrd_taxa_last.tsv")
out_file <- "ftrd_last.xml"
xml_editor(template_file, out_dir, out_file, taxa_file, fossils = TRUE)

# out directory for extant scripts
out_dir <- extant_scripts_dir

# molecular data
template_file <- paste0(extant_scripts_dir, "mol_recent_template.xml")
out_file <- "mol_recent.xml"
xml_editor(template_file, out_dir, out_file)

# morphological data
template_file <- paste0(extant_scripts_dir, "morpho_recent_template.xml")
out_file <- "morpho_recent.xml"
xml_editor(template_file, out_dir, out_file)

# combined data
template_file <- paste0(extant_scripts_dir, "comb_recent_template.xml")
out_file <- "comb_recent.xml"
xml_editor(template_file, out_dir, out_file)
