---
title: "README.md"
author: "Yuan Zhen"
date: "2022-10-17"
output: html_document
---

#### an R package to fast find Genome-wide gene family

You can find the information you are interested in from the massive data of NCBI. I have made full use of the R packages: rentrez, tidyverse, biostrings, Peptides, and ggrepl to quickly find specific gene families in the species you submit. The following is my implementation process:

```R
library(rentrez)
library(tidyverse)
library(Biostrings)
library(Peptides)
library(ggrepel)
library(proEvlu)

### you can get the pfam ID or the taxonomy ID from NCBI
pfm=get_pfm("VQ")
my_tax=get_taxonomy("Ananas comosus")

###you shold find the correct pfam and taxonomy ID,then you can get the protein accession ID
xm=proEvlu::get_protein(pfam = "pfam05678",taxonomyID = 4615)

#### to get the chr and gene information
yz=proEvlu::gene_sturcture(yy,my_tax$taxonomy_id)

#### Visualize the location of genes on chromosomes
proEvlu::plot_genes(yz)

##### get the protein sequence
yzxm=proEvlu::protein_seq(xm)

####Information on the physicochemical properties of proteins is obtained
proEvlu::protein_properties(yzxm)
```

This is a basic usage. It can also be combined with other R packages (an enhanced version will be developed later) to carry out downstream analysis, such as gene collinearity analysis, evolutionary tree analysis, gene structure analysis, conservative structure analysis, etc
