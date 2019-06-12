getGeneByGOTerms <- function(term, orgDb){
  go2egdb <- sub(".db", "GO2ALLEGS", orgDb)
  e <- mget(term, get(go2egdb))
  eg <- lapply(e, unique)
  egs <- unlist(eg)
  egs <- unique(egs)
  sym <- mget(egs, get(sub("GO2ALLEGS", "SYMBOL", go2egdb)), ifnotfound = NA)
  sym <- unique(unlist(sym))
  sym <- sym[!is.na(sym)]
  return(sym)
}