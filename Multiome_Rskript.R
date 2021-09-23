#############
# Install Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install libraries
BiocManager::install(c("Signac","dplyr", "org.Mm.eg.db", "ggplot2", "GenomeInfoDb", "patchwork", "EnsDb.Mmusculus.v79"))
install.packages('Seurat')
install.packages('hdf5r') # Might not run on Mac; Error: 
install.packages("EnsDb.Mmusculus.v79")
install.packages("stringr")
install.packages("rtracklayer")

#Package which is only available in source form, and may need compilation of C/C++/Fortran: ‘hdf5r’
#Do you want to attempt to install these from sources? (Yes/no/cancel) ; Recquires installation of Xcode.

#install.packages("devtools")
#devtools::install_github("hhoeflin/hdf5r")

#install.packages()

# Libraries #

# used packages
# if not installed please install with either bioconductor BiocManager::install() (see useful links) or install.packages()



library(Seurat)
library(Signac)
library(dplyr)
library(org.Mm.eg.db)
library(ggplot2)
library(hdf5r)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(rtracklayer)
library(stringr)
SeqinfoForUCSCGenome("mm10")



#We’ll create a Seurat object based on the gene expression data, and then add in the 
#ATAC-seq data as a second assay. You can explore the Signac getting started vignette 
#for more information on the creation and processing of a ChromatinAssay object.

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("D:/Google Drive/ILCDEV/R/Multiome/Multiome_Seurat/filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

?standardChromosomes

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts, species = "Mus_musculus")
atac_counts <- atac_counts[as.vector(grange.use), ]

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style since the data was mapped to hg19

#seqlevelsStyle(annotations) <- 'UCSC' #changed because of bug of GenomeInfoDb https://github.com/Bioconductor/GenomeInfoDb/issues/27

ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotations) <- ucsc.levels
genome(annotations) <- 'mm10'

frag.file <- "D:/Google Drive/ILCDEV/R/Multiome/Multiome_Seurat/atac_fragments.tsv.gz"

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'BSgenome.Mmusculus.UCSC.mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

pbmc[["ATAC"]] <- chrom_assay






head(metadata)
#References
#https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
#https://satijalab.org/signac/articles/install.html
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html

