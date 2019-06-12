pkg <- installed.packages()
pkg <- rownames(pkg)
if(!"biomaRt" %in% pkg){
  if(!"BiocManager" %in% pkg){
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(c("biomaRt", "Biostrings"))
}

library(biomaRt)
filename <- "genelist.txt"
template <- "example.txt"
outfolder <- "output"
primer3path <- "primer3_core"
speciesdb <- c("human"="hsapiens", "mouse"="mmusculus", "zebrafish"="drerio")
ids <- c("human"="hgnc_symbol", "mouse"="mgi_symbol", "zebrafish"="zfin_id_symbol")
species <- "human"
barcode="ACGT"
UMI="ACGT"
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  message("Use default filename: genelist.txt; species: human; template: example.txt;",
          "outfolder: output; primer3path: primer3_core; barcode: ACGT; UMI: ACGT")
}else{
  for(i in seq_along(args)){
    print(args[[i]])
    eval(parse(text=args[[i]]))
  }
}

## create out folder
dir.create(outfolder, recursive = TRUE)
## read template
param <- paste(readLines(template), collapse = "\n")
## read gene list
genes <- readLines(filename)
## use mart
mart <- useMart("ensembl", dataset=paste0(speciesdb[species], "_gene_ensembl"))
## get cDNA sequence
seq <- getSequence(id=genes, type=ids[species], seqType="cdna", mart=mart)
seq <- seq[order(seq[, 2]), ]
## set unique id
id <- rle(seq[, 2])
seq$id <- paste(seq[, 2], unlist(sapply(id$lengths, seq.int)), sep="_")

## prepare parameters
library(Biostrings)
primer3end <- DNAString(paste0("TTTTTTTAAGCAGTGGTATCAACGCAGAGTAC",
                               barcode,
                               UMI,
                               "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"))
primer3end.rc <- reverseComplement(primer3end)
seqlen <- nchar(seq$cdna)
seqlen0 <- seqlen -200
seqlen0[seqlen0<50] <- 50
seq$cdna <- paste0(seq$cdna, primer3end.rc)
output <- file.path(outfolder, "p3config.txt")
seq_id <- paste0("SEQUENCE_ID=", seq$id)
seq_tmp <- paste0("SEQUENCE_TEMPLATE=", seq$cdna)
primer3 <- paste0("PRIMER_MUST_MATCH_THREE_PRIME=nnn", primer3end)
product_size <- paste0("PRIMER_PRODUCT_SIZE_RANGE=", 
                       ifelse(seqlen>1000, 
                              "800-2000", paste0(seqlen0, "-", seqlen)))
outStr <- rbind(seq_id, seq_tmp, primer3, product_size, rep(param, length(seq_id)))
outStr <- as.character(outStr)
writeLines(outStr, output)

## run primer3

out <- system(paste0(primer3path, " --output=", file.path(outfolder, "out.txt"),
                     "  --error=", file.path(outfolder, "err.txt"), " ", output))

if(out!=0){
  stop("error code: ", out)
}

## parse the results
primers <- readLines(file.path(outfolder, "out.txt"))
## split the primer for each cDNA
breaker <- which(primers=="=")
breaker <- diff(c(0, breaker))
bf <- rep(seq_along(breaker), breaker)
stopifnot(length(bf)==length(primers))
primers <- split(primers, bf)
## make as a table
primers <- lapply(primers, function(.ele){
  .ele <- do.call(rbind, strsplit(.ele, "="))
  .out <- .ele[, 2]
  names(.out) <- .ele[, 1]
  .out
})

coln <- unique(unlist(lapply(primers, names)))
rown <- unlist(lapply(primers, `[`, "SEQUENCE_ID"))
mat <- matrix(NA, nrow = length(primers), ncol = length(coln), 
              dimnames = list(rown, coln))
for(i in seq_along(rown)){
  for(j in seq_along(coln)){
    if(coln[j] %in% names(primers[[i]])) mat[i, j] <- primers[[i]][coln[j]]
  }
}

write.csv(mat, file.path(outfolder, "primers.csv"))
