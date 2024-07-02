rm(list = ls())
library(stringr)

base_path <- "/Users/petrucci/Documents/research/canids_chap2/beast/"
template_path <- paste0(base_path, "scripts/canidae_template.xml")
ranges_path <- paste0(base_path, "data/canidae_ranges.tsv")
run_path <- paste0(base_path, "scripts/canidae_run.xml")

file.remove(run_path)

ranges <- read.delim(ranges_path)
# tmp_id <- which(ranges$la_max>ranges$fa_min)
# ranges$la_max[tmp_id]<-ranges$fa_min[tmp_id]

uncertain_fa <- which(ranges$fa_max-ranges$fa_min!=0)
uncertain_la <- which(ranges$la_max-ranges$la_min!=0)
certain_fa <- which(ranges$fa_max-ranges$fa_min==0)
certain_la <- which(ranges$la_max-ranges$la_min==0)

xml <- readLines(template_path)
idx  <- grep('totalcount="5"', xml)
morph <- xml[idx]
morph_first <- gsub('" spec', '_first" spec', morph)
morph_first <- gsub('" totalcount', '_first" totalcount', morph_first)

morph_last <- gsub('" spec', '_last" spec', morph)
morph_last <- gsub('" totalcount', '_last" totalcount', morph_last)

idx_dna <- grep('totalcount="4"', xml)
dna <- xml[idx_dna]
dna_first <- gsub('" spec', '_first" spec', dna)
dna_first <- gsub('" totalcount', '_first" totalcount', dna_first)
dna_first <- gsub('value=".*"', paste0('value="',strrep('-', 16570),'"'), dna_first)

dna_last <- gsub('" spec', '_last" spec', dna)
dna_last <- gsub('" totalcount', '_last" totalcount', dna_last)

# copy_xml <- paste(xml[1:(idx[1]-1)], morph_first, morph_last,
#               xml[(idx[length(idx)]+1):(idx_dna[1]-1)], dna_first, dna_last,
#               xml[(idx_dna[length(idx_dna)]+1):length(xml)], collapse = '\n')

write(xml[1:(idx[1]-1)],file=run_path,append=TRUE)
write(morph_first,file=run_path,append=TRUE)
write(morph_last,file=run_path,append=TRUE)
write(xml[(idx[length(idx)]+1):(idx_dna[1]-1)],file=run_path,append=TRUE)
write(dna_first,file=run_path,append=TRUE)
write(dna_last,file=run_path,append=TRUE)
write(xml[(idx_dna[length(idx_dna)]+1):length(xml)],file=run_path,append=TRUE)

rm(xml)
copy_xml <- readLines(run_path)




ranges_str <-character(nrow(ranges))
ranges_idref_str <-character(nrow(ranges))
dates_first <- character(nrow(ranges))
dates_last <- character(nrow(ranges))
for (i in 1:nrow(ranges)){
  ranges_str[i] <- paste0('<stratigraphicRange id="r',i,
  '" spec="StratigraphicRange" firstOccurrence="@',ranges$taxon[i],
  '_first" lastOccurrence="@',ranges$taxon[i],'_last"/>')
  ranges_idref_str[i] <- paste0('<stratigraphicRange idref="r',i,'"/>')
  
  f <- runif(1, ranges$fa_min[i], ranges$fa_max[i])
  l <- runif(1, ranges$la_min[i], ranges$la_max[i])
  while (f<l){
    f <- runif(1, ranges$fa_min[i], ranges$fa_max[i])
    l <- runif(1, ranges$la_min[i], ranges$la_max[i])
  }
  dates_first[i] <- paste0(ranges$taxon[i],'_first = ',f)
  dates_last[i] <- paste0(ranges$taxon[i],'_last = ',l)
}

ranges_str <-str_c(ranges_str, collapse = '\n')
ranges_idref_str <-str_c(ranges_idref_str, collapse = '\n')
dates_str <-str_c(c(dates_first,dates_last), collapse = ',\n')

copy_xml<-gsub('<insertDates/>', dates_str, copy_xml)
copy_xml<-gsub('<insertRanges/>', ranges_str, copy_xml)
copy_xml<-gsub('<insertRangeIds/>', ranges_idref_str, copy_xml)


certain_taxon<- character(length(certain_fa)+length(certain_la))
for (i in 1:length(certain_fa)){
  certain_taxon[i] <- paste0('<taxon id="',ranges$taxon[certain_fa[i]],'_first" spec="Taxon"/>')
}  
for (i in (length(certain_fa)+1):(length(certain_fa)+length(certain_la))){
  certain_taxon[i] <- paste0('<taxon id="',ranges$taxon[certain_la[i-length(certain_fa)]],'_last" spec="Taxon"/>')
}
certain_taxon <-str_c(certain_taxon, collapse = '\n')
copy_xml<-gsub('<insertCertainTaxon/>', certain_taxon, copy_xml)

uncertain_taxon<- character(length(uncertain_fa)+length(uncertain_la))
uncertainty_ranges <- character(length(uncertain_fa)+length(uncertain_la))
for (i in 1:length(uncertain_fa)){
  uncertain_taxon[i] <- paste0('<taxon id="',ranges$taxon[uncertain_fa[i]],'_first" spec="Taxon"/>')
  uncertainty_ranges[i] <- paste0('<samplingDates id="samplingDate',i,
                                  '" spec="sa.evolution.tree.SamplingDate" taxon="@',
                                  ranges$taxon[uncertain_fa[i]],'_first" lower="',
                                  ranges$fa_min[uncertain_fa[i]],'" upper="',
                                  ranges$fa_max[uncertain_fa[i]],'"/>')
}

for (i in (length(uncertain_fa)+1):(length(uncertain_fa)+length(uncertain_la))){
  uncertain_taxon[i] <- paste0('<taxon id="',ranges$taxon[uncertain_la[i-length(uncertain_fa)]],'_last" spec="Taxon"/>')
  uncertainty_ranges[i] <- paste0('<samplingDates id="samplingDate',i,
                                  '" spec="sa.evolution.tree.SamplingDate" taxon="@',
                                  ranges$taxon[uncertain_la[i-length(uncertain_fa)]],'_last" lower="',
                                  ranges$la_min[uncertain_la[i-length(uncertain_fa)]],'" upper="',
                                  ranges$la_max[uncertain_la[i-length(uncertain_fa)]],'"/>')
}

uncertain_taxon <-str_c(uncertain_taxon, collapse = '\n')
uncertainty_ranges <-str_c(uncertainty_ranges, collapse = '\n')

copy_xml <- gsub('<insertTaxonForRandomWalk/>', uncertain_taxon, copy_xml)
copy_xml <- gsub('<insertRandomWalkRangeBounds/>', uncertainty_ranges, copy_xml)


writeLines(copy_xml, run_path)

