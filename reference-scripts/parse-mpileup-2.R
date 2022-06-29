library(tidyverse)
library(readxl)
library(optparse)

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-a", "--all_muts"), type="character", default=NULL,
              help="comma-separated file with all mutations in case", metavar="character"),
  make_option(c("-i", "--indel_parsed"), type="character", default=NULL,
              help="indel_parsed.tsv", metavar="character"),
  make_option(c("-m", "--mbc_parsed"), type="character", default=NULL,
              help="mbc_parsed.tsv", metavar="character"),
  make_option(c("-s", "--snv_parsed"), type="character", default=NULL,
              help="snv_parsed.tsv", metavar="character"),
  make_option(c("-c", "--CGID_sample"), type="character", default=NULL,
              help="CGID of sample", metavar="character"),
  make_option(c("-d", "--mutdata"), type="character", default=NULL,
              help="pgdx mutation data for sample", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #

fs::dir_create(opt$outdir)

all_muts <- read.csv(opt$all_muts) %>%
  tibble()
sample_id <- opt$CGID_sample
mut_data <- read.delim(opt$mutdata) %>%
  tibble()

if (ncol(mut_data) == 3) message("no read count info in mutdata")

if (colnames(mut_data)[1] == "gene") {
  colnames(mut_data)[1:2] <- c("Gene.Symbol", "Mutation.Position")
}

chr_order <- paste0("chr", c(1:22, "X", "Y"))

source("parse-mpileup-functions.R")

# ===================================================================================== #

sample_out <- all_muts %>%
  bind_cols(., tibble(CGID_sample = sample_id,
                      Distinct.Total.Reads = NA,
                      Distinct.Mutant.Reads = NA)) %>%
  arrange(chr)

# add distinct total and mutant reads from mut data to sample out 
#mutdata_sampleout_ind <- match(mut_data$Mutation.Position, sample_out$MutID)
#sample_out$Distinct.Total.Reads[mutdata_sampleout_ind] <- mut_data$Distinct.Total.Reads
#sample_out$Distinct.Mutant.Reads[mutdata_sampleout_ind] <- mut_data$Distinct.Mutant.Reads

missing_mut_inds <- which(is.na(sample_out$Distinct.Mutant.Reads))

if (length(missing_mut_inds) > 0) {
  # classify mutation type (snv, indel, mbc) for missing mutations
  mut_type_inds <- determineMutTypes(sample_out[missing_mut_inds, ])
  
  # SNVs
  if (file.exists(opt$snv_parsed)) {
    snv_mpileup <- read.delim(opt$snv,
                              stringsAsFactors = F)
    
    if (length(mut_type_inds$snv_ind) > 0) {
      missing_snv_ind <- missing_mut_inds[mut_type_inds$snv_ind]
      sample_out <- addSNVMpileupParsed(snv_mpileup, sample_out, missing_snv_ind)
    } else {
      message("no snvs to add")
    }
  }
  
  # Indels
  if (file.exists(opt$indel_parsed)) {
     if (length(mut_type_inds$indel_ind) > 0) {
    
    indel_mpileup <- read.delim(opt$indel_parsed,
                                stringsAsFactors = F,
                                colClasses = c("character",
                                               "numeric",
                                               "character",
                                               "integer",
                                               "integer",
                                               "integer",
                                               "integer",
                                               "integer",
                                               "integer",
                                               "character",
                                               "character",
                                               "character"))
    
    missing_indel_ind <- missing_mut_inds[mut_type_inds$indel_ind]
    sample_out <- addIndelMpileupParsed(indel_mpileup, sample_out, missing_indel_ind)
   } else {
      message("no indels to add")
   }
  }
  
  # Multi-base change muts
  if (file.exists(opt$mbc_parsed)) {

    if (length(mut_type_inds$mbc_ind) > 0) {
    mbc_mpileup <- read.delim(opt$mbc_parsed,
                              stringsAsFactors = F)
    
    missing_mbc_ind <- missing_mut_inds[mut_type_inds$mbc_ind]
    sample_out <- addMBCMpileupParsed(mbc_mpileup, sample_out, missing_mbc_ind)
    } else {
      message("no mbc to add")
    } 
  }
}


write.table(sample_out, 
            file.path(opt$outdir, paste0(sample_id, "_counts.tsv")),
            quote = F, sep = "\t", row.names = F, col.names = T)





