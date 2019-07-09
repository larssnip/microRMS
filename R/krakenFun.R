

readKrakenReport <- function(krak.file, aggregate.to = "species"){
  patt <- switch(aggregate.to, "strain" = "S1", "species" = "S", "genus" = "G")
  suppressMessages(read_delim(file.path(krk.dir, fnam), delim = "\t", trim_ws = T,
                              col_names = c("Percentage", "Clade.count", "Tax.count", "Rank", "Taxid", "Name"))) %>% 
    filter(Rank == patt) %>% 
    select(Taxon = Name, Taxid, Rank, Clade.count) -> tbl
  return(tbl)
}