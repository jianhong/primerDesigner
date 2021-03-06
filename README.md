# primerDesigner
pipeline for primer3

# Howto

## requirements

 [primer3](http://primer3.org/manual.html#installLinux)
 
 [R](https://www.r-project.org/) and package [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html), [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

## install

```
git clone https://github.com/jianhong/primerDesigner.git
```

## template file

Template file is the input file by removing the SEQUENCE_ID and SEQUENCE_TEMPLATE. Default is example.txt.
Last line must be `=`. Reference: http://primer3.org/manual.html#globalTags

For example, if you want to change the number of primers, reset the line of PRIMER_NUM_RETURN=5.
If you want to change the primer TM value to 42, you can add three lines in the template file:
```
PRIMER_OPT_TM=42.0
PRIMER_MIN_TM=39.0
PRIMER_MAX_TM=45.0
```

## gene list file

The target gene symbol list. One symbol per line. Default is genelist.txt.
if GO term is supplied, set golist to TRUE.

## run scipt

```
R CMD BATCH --no-save --no-restore '--args species="human"' primerDesigner.R out.log.txt &
```

By GO term
```
R CMD BATCH --no-save --no-restore '--args species="human" filename="goterms.txt" golist=TRUE' primerDesigner.R out.log.txt &
```

Run it in slurm
```
srun --mem=4GB R CMD BATCH --no-save --no-restore '--args species="human" filename="goterms.txt" golist=TRUE outfolder="output2"' primerDesigner.R out.log.txt &
```

## arguments

species, could be human, mouse, zebrafish, default `human`.

filename, the gene list file name, default `genelist.txt`.

template, the template file name, default `example.txt`.

outfolder, output folder, default `output`.

primer3path, the path to primer3_core, default `primer3_core`.

## error and log

Log file is out.log.txt for primerDesigner.R and ${outfolder}/out.txt for primer3. 
Error file for primer3 is ${outfolder}/err.txt.
