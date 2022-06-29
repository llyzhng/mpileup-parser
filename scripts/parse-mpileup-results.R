library(optparse)
library(parallel)

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

# SNV parser

addSNVMpileupParsed <- function(snv_mpileup, sample_out) {
  # columns of parsed pileup are: chrom pos     ref     cov     A       C       G       T       *       -       + X
  #snv_mpileup <- snv_mpileup %>%
  #  mutate(chr_start = paste(chrom, pos, sep = "_"))
  
  #sample_out <- sample_out %>%
  #  mutate(chr_start = paste(CHROM, POS, sep = "_"))
  
  snv_mpileup$chr_start <- paste(snv_mpileup$chrom, snv_mpileup$pos, sep = "_")
  sample_out$chr_start <- paste(sample_out$CHROM, sample_out$POS, sep = "_")
  
  sample_out_list <- mclapply(1:nrow(sample_out), 
                         function(ind) updateRow(snv_mpileup, sample_out, ind),
                         mc.cores = detectCores())
  sample_out <- do.call("rbind", sample_out_list)
  
  sample_out$chr_start <- NULL
  return(sample_out)
}
#source("parse-mpileup-functions.R")

updateRow <- function(snv_mpileup, sample_out, j) {
  temp_ref <- sample_out$REF[j]
  temp_alt <- sample_out$ALT[j]
  temp_mpileup <- snv_mpileup[match(sample_out$chr_start[j], snv_mpileup$chr_start), ]
  
  if (nrow(temp_mpileup) > 1) stop("something is wrong?? matched to two mpileup results")
  
  mpileup_tot_count <- temp_mpileup$cov[1]
  mpileup_alt_count <- temp_mpileup[1, temp_alt]
  sample_out$variant_reads[j] <- mpileup_alt_count
  sample_out$total_reads[j] <- mpileup_tot_count
  return(sample_out[j, ])
}

# ===================================================================================== #

all_muts <- read.delim(opt$all_muts, stringsAsFactors = F)

sample_out <- all_muts
sample_out$total_reads <- NA
sample_out$variant_reads <- NA

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

sample_out_annot <- addSNVMpileupParsed(snv_mpileup, sample_out)


write.table(sample_out_annot, 
            opt$outfile,
            quote = F, sep = "\t", row.names = F, col.names = T)



