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
              help="snv_mpileup.parsed.txt", metavar="character"),
  make_option(c("c", "--cores"), type="numeric", default=NULL,
              help="number of cores", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #
# ===================================================================================== #

pickCol <- function(nt) {
  if (nt == "A") return (5)
  if (nt == "C") return (6)
  if (nt == "G") return (7)
  if (nt == "T") return (8)
}

# ===================================================================================== #


if (is.null(opt$cores)) {
  opt$cores <- detectCores()
}

sample_out <- read.delim(opt$all_muts, stringsAsFactors = F)


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

snv_mpileup$chr_start <- paste(snv_mpileup$chrom, snv_mpileup$pos, sep = "_")
sample_out$chr_start <- paste(sample_out$CHROM, sample_out$POS, sep = "_")

inds <- match(sample_out$chr_start, snv_mpileup$chr_start)
pileup_sub <- snv_mpileup[inds, ]

sample_out$total_reads <- pileup_sub$cov
alt.col.ind <- sapply(sample_out$ALT, pickCol)

alt.count <- mcmapply(function(i,j) pileup_sub[i,j],
                      i = 1:nrow(sample_out),
                      j = alt.col.ind,
                      mc.cores = opt$cores)
sample_out$variant_reads <- alt.count

write.table(sample_out, opt$outfile, 
            quote = F, sep = "\t", row.names = F, col.names = T)
