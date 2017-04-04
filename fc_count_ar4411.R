#run this script using : run_fc_script.sh
setwd("/scratch/ar4411/Homework02")


#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsubread",lib.loc="/scratch/ar4411/Homework02")
#Rsubreads does not work on RStudio, therefore run this script as an Rscript, using bash script.
library(edgeR)		# need this for rpkm 
library(Rsubread)	# main library for featureCounts
library(methods)   	# (also need to include this library when calling on command-line using Rscript)


# featureCounts currently does not have a function to generate Tpm, so we convert from fpkm to tpm


#TPM:TPM stands for TPM "Transcripts per Kilobase of exons per million reads mapped".
#TPM is applied to normalize the sequencing depth of samples by accounting for the number of transcripts produced, for the samples that are being compared.
#"These approaches rely on normalizing methods that are based on total or effective counts, and tend to perform poorly when samples have heterogeneous 
#transcript distributions, that is, when highly and differentially expressed features can skew the count distribution".(Conesa, p-8)

fpkmToTpm <- function(fpkm, dolog=T)
{
  tpm = log2(fpkm) - log2(sum(fpkm)) + log2(1e6)
  if (!dolog) 
  {
    tpm = (2^tpm)
  }
  tpm	
}


#get this from scratch/ar4411/Homework02
annof  =  "mm10_genes.gtf"

# Upload all 8 bam files as a list. The file extension is .sorted.bam. Not using the indexed files (.bai)
bamfiles= list.files(pattern="*sorted.bam$")

# featureCounts is the main function to give the counts of reads over each gene / exon / genomic feature in the gtf file.
# going to use this funciton for both gene and exon counts.

{
  
  fc_PE_gene <- featureCounts(bamfiles, 
                              annot.ext = annof,  # gtf file
                              GTF.attrType="gene_id", # the label in the gtf describing a gene 
                              nthreads=6,   # use multiple cores for speed
                              minMQS=20,    # exclude multimapped reads (low MAPQ score) 
                              countMultiMappingReads = F,
                              strandSpecific = 0,    # this is not strand-specific data, we are working with both strands(would set to 1 or 2 otherwise).
                              isGTFAnnotationFile=T,
                              GTF.featureType = "exon",
                              useMetaFeatures=TRUE, # This indicates we want expression over genes; for exon-level expression set this to FALSE 
                              primaryOnly=TRUE,  # only include primary mappings
                              isPairedEnd=TRUE) # this is pair-ended data
}
{
  fc_PE_exon <- featureCounts(bamfiles, 
                              annot.ext = annof,  # gtf file
                              GTF.attrType="gene_id", # the label in the gtf describing a gene 
                              nthreads=5,   # use multiple cores for speed
                              minMQS=20,    # exclude multimapped reads (low MAPQ score) 
                              countMultiMappingReads = F,
                              strandSpecific = 0,    # this is not strand-specific data (would set to 1 or 2 to take strand into account for Illumina data, depending on strand-specific library construction)
                              isGTFAnnotationFile=T,
                              GTF.featureType = "exon",
                              useMetaFeatures=FALSE, # This indicates we want expression over genes; for exon-level expression set this to FALSE 
                              primaryOnly=TRUE,  # only include primary mappings
                              isPairedEnd=TRUE) # this is pair-ended data
}

# fc_PE$counts are the raw read counts overlying the features
# Note that geneids are given as the row names.
head(rownames(fc_PE_gene$counts))
head(colnames(fc_PE_gene$counts))

head(rownames(fc_PE_exon$counts))
head(colnames(fc_PE_exon$counts))

# fc_PE$annotation stores extra information on gene location etc.

# tidy up the column names by removing the "sorted.bam" extension, so that exon and gene tables only show the replicate name as column headers.
colnames(fc_PE_gene$counts) = gsub(".sorted.bam","", colnames(fc_PE_gene$counts), perl=T)
head(colnames(fc_PE_gene$counts))

colnames(fc_PE_exon$counts) = gsub("sorted.bam","", colnames(fc_PE_exon), perl = T)
head(colnames(fc_PE_exon$counts))

{
  
  # For later differential gene expression analysis, all we need is the counts
  # Here we calculate the expression in units of TPM to order genes by expression and compare expression of genes within one sample.
  # edgeR has the rpkm function to calculate the (obsolete) RPKM unit, so we use that and then convert to TPM.
  x_PE <- DGEList(counts=fc_PE_gene$counts, genes=fc_PE_gene$annotation)
  x_PE_rpkm = rpkm(x_PE,x_PE$genes$Length, log=T)
  x_PE_rpkm = data.frame(x_PE_rpkm)
  
  # convert from RPKM to TPM
  x_tpm  = fpkmToTpm(2^x_PE_rpkm, dolog=F)
  #end of gene level
  
  
  # save results in a binary RData file for later processing 
  save(fc_PE_gene, x_tpm, file="ar4411_gene_2.RData")
  
}
# exon-level
{
  # (Note: TPM is meaning less for exon-level counts)
  save(fc_PE_exon, file="ar4411_exon_2.RData")

  # we will use these RData files as input for the .Rmd files later on.
  
}