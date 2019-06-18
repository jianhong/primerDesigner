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
golist <- FALSE
template <- "example.txt"
outfolder <- "output"
primer3path <- "primer3_core"
speciesdb <- c("human"="hsapiens", "mouse"="mmusculus", "zebrafish"="drerio")
ids <- c("human"="hgnc_symbol", "mouse"="mgi_symbol", "zebrafish"="zfin_id_symbol")
orgDbs <- c("human"="org.Hs.eg.db", "mouse"="org.Mm.eg.db", "zebrafish"="org.Dr.eg.db")
species <- "human"
barcode="ACGT"
UMI="ACGT"
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  message("Use default filename: genelist.txt; species: human; template: example.txt;",
          "outfolder: output; primer3path: primer3_core; barcode: ACGT; UMI: ACGT;",
          "golist: FALSE")
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
if(golist){
  ## check org db
  orgDb <- orgDbs[species]
  if(!orgDb %in% pkg){
    BiocManager::install(orgDb)
  }
  library(orgDb, character.only = TRUE)
  source("getGeneByGoTerms.R")
  genes <- getGeneByGOTerms(genes, orgDb)
}
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
primer3 <- "AAGCAGTGGTATCAACGCAGAGT"
primer3 <- paste0("SEQUENCE_PRIMER_REVCOMP=", primer3)
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

gene_name <- sub("_\\d+$", "", mat[, "SEQUENCE_ID"])
out.mat <- list()
k <- 1
for(i in seq_along(gene_name)){
  J <- as.numeric(mat[i, "PRIMER_LEFT_NUM_RETURNED"])
  if(J>0){
    for(j in seq.int(J)){
      sequence <- mat[i, paste0("PRIMER_LEFT_", j-1, "_SEQUENCE")]
      primer_length <- nchar(sequence)
      product_length <- nchar(mat[i, "SEQUENCE_TEMPLATE"]) -
        as.numeric(strsplit(mat[i, paste0("PRIMER_LEFT_", j-1)], ",")[[1]][1]) + 1
      primer_TM_value <- mat[i, paste0("PRIMER_LEFT_", j-1, "_TM")]
      out.mat[[k]] <- c(mat[i, "SEQUENCE_ID"], gene_name[i], paste0(mat[i, "SEQUENCE_ID"], "_", j),
                        sequence, primer_length, product_length, primer_TM_value)
      k <- k+1
    }
  }
}
out.mat <- do.call(rbind, out.mat)
colnames(out.mat) <- c("sequence_id", "gene_name", "primer_id", "sequence", 
                       "primer_length", "product_length", "primer_TM_value")
write.csv(out.mat, file.path(outfolder, "summary.csv"), row.names = FALSE)