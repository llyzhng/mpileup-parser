library(tidyverse)
library(optparse)

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="full path of output file", metavar="character"),
  make_option(c("-a", "--all_muts"), type="character", default=NULL,
              help="tab-separated file with chr, position, ref, and alt alleles for all mutations", metavar="character"),
  make_option(c("-s", "--snv"), type="character", default=NULL,
              help="snv_mpileup.parsed.txt", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #
# ===================================================================================== #
# SNV parser

addSNVMpileupParsed <- function(snv_mpileup, sample_out) {
  # columns of parsed pileup are: chrom	pos	ref	cov	A	C	G	T	*	-	+ X
  snv_mpileup <- snv_mpileup %>%
    mutate(chr_start = paste(chrom, pos, sep = "_"))
  
  sample_out <- sample_out %>%
    mutate(chr_start = paste(CHROM, POS, sep = "_"))
  
  for (j in 1:nrow(sample_out)) {
    temp_ref <- sample_out$REF[j]
    temp_alt <- sample_out$ALT[j]
    temp_mpileup <- snv_mpileup[match(sample_out$chr_start[j], snv_mpileup$chr_start), ]
    
    if (nrow(temp_mpileup) > 1) stop("something is wrong?? matched to two mpileup results")
    
    mpileup_tot_count <- temp_mpileup$cov[1]
    mpileup_alt_count <- temp_mpileup[1, temp_alt]
    sample_out$variant_reads[j] <- mpileup_alt_count
    sample_out$total_reads[j] <- mpileup_tot_count
  }
    
  sample_out <- sample_out 
    
  sample_out$chr_start <- NULL
  return(sample_out)
}
#source("parse-mpileup-functions.R")

# ===================================================================================== #

all_muts <- read.delim(opt$all_muts) %>%
  tibble()

sample_out <- all_muts %>%
  bind_cols(., tibble(total_reads = NA,
                      variant_reads = NA))

snv_mpileup <- read.delim(opt$snv,
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

sample_out <- addSNVMpileupParsed(snv_mpileup, sample_out)


write.table(sample_out, 
            opt$outfile,
            quote = F, sep = "\t", row.names = F, col.names = T)





