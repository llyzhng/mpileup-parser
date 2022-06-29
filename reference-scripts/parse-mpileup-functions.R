determineMutTypes <- function(sample_out) {
  # determine mutation type
  indel_ind <- which(nchar(sample_out$ref) != nchar(sample_out$alt) | 
                       is.na(sample_out$alt))
  # indices of multibase change mutations 
  multibasechange_ind <- which((nchar(sample_out$ref) > 1 | nchar(sample_out$alt) > 1)
                               & (nchar(sample_out$ref) == nchar(sample_out$alt)))
  ### SNVs only
  not_snv_ind <- c(indel_ind, multibasechange_ind)
  snv_ind <- which(!(1:nrow(sample_out) %in% not_snv_ind))
  
  return(list(indel_ind = indel_ind,
              mbc_ind = multibasechange_ind,
              snv_ind = snv_ind))
}

check_insertion <- function(ref) {
  if (is.na(ref) | ref == "") return(TRUE)
  return(FALSE)
}

addSNVMpileupParsed <- function(snv_mpileup, sample_out, snv_ind) {
  # columns of parsed pileup are: chrom	pos	ref	cov	A	C	G	T	*	-	+ X
  if (length(snv_ind) > 0) {
    snv_mpileup <- snv_mpileup %>%
      mutate(chr_start = paste(chrom, pos, sep = "_"))
    
    sample_out <- sample_out %>%
      mutate(chr_start = paste(chr, start, sep = "_"))
    
    for (j in snv_ind) {
      temp_ref <- sample_out$ref[j]
      temp_alt <- sample_out$alt[j]
      temp_mpileup <- snv_mpileup[match(sample_out$chr_start[j], snv_mpileup$chr_start), ]
      
      if (nrow(temp_mpileup) > 1) stop("something is wrong?? matched to two mpileup results")
      
      mpileup_tot_count <- temp_mpileup$cov[1]
      mpileup_alt_count <- temp_mpileup[1, temp_alt]
      sample_out$Distinct.Mutant.Reads[j] <- mpileup_alt_count
      sample_out$Distinct.Total.Reads[j] <- mpileup_tot_count
    }
    
    sample_out <- sample_out %>%
      select(Gene.Symbol, MutID, chr, start_stop, start, stop, ref, alt, CGID_sample, Distinct.Total.Reads, Distinct.Mutant.Reads)
    
  } else {
    message("no snv_ind supplied")
  }
  sample_out$chr_start <- NULL
  return(sample_out)
}

addIndelMpileupParsed <- function(indel_mpileup, sample_out, indel_ind){
  
  indel_mpileup <- indel_mpileup %>%
    mutate(chr_start = paste(chrom, pos, sep = "_"))
  
  sample_out <- sample_out %>%
    mutate(chr_start = paste(chr, start, sep = "_"))
  
  for (j in indel_ind) {
    # insertion or deletion? 
    is_ins <- check_insertion(sample_out$ref[j])
    if (is_ins) {
      sample_out <- addInsertion(sample_out, j, indel_mpileup)
    } else {
      sample_out <- addDeletion(sample_out, j, indel_mpileup)
    }
  }
  sample_out$chr_start <- NULL
  return(sample_out)
}

addDeletion <- function(sample_out, j, indel_mpileup) {
  # parsed mpileup column * (X.) is deletion of reference base
  # column - (X..1) is deletion 
  # when deletion is in "-" columns, deletion present at prev position in mpileup
  check_row_minus <- indel_mpileup %>%
    filter(pos == sample_out[j,]$start - 1)
  check_row_star <- indel_mpileup %>%
    filter(pos == sample_out[j,]$start)
  
  # check previous position first
  prev_pos_del <- check_row_minus[, "X..1"]
  if (prev_pos_del != "") {
    # deletion present in "-" column at previous position
    base_tb <- table(strsplit(prev_pos_del, ",")[[1]])
    if (nrow(base_tb) > 1) {
      message(paste0("j is ", j))
      message(paste(sample_out[j, ], collapse = ","))
      message("warning: more than one deleted string")
      message(print(base_tb))
      #stop("something weird... more than 1 deletion at this position?")
      
      del_strings <- names(base_tb)
      deleted_base <- sample_out[j,]$ref
      del_count <- unname(base_tb)[which(del_strings == deleted_base)]
      del_cov <- check_row_minus$cov
    } else {
      del_count <- unname(base_tb)
      del_bases <- names(base_tb)
      del_cov <- check_row_minus$cov
    }
    
    sample_out$Distinct.Total.Reads[j] <- del_cov
    sample_out$Distinct.Mutant.Reads[j] <- del_count
    return(sample_out)
  }
  
  # check deleted position
  pos_del <- check_row_star[, "X."]
  del_count <- pos_del
  del_cov <- check_row_star$cov
  sample_out$Distinct.Total.Reads[j] <- del_cov
  sample_out$Distinct.Mutant.Reads[j] <- del_count
  return(sample_out)
}



addDeletion2 <- function(sample_out, j, indel_mpileup) {
  # parsed mpileup column * (X.) is deletion of reference base
  # column - (X..1) is deletion 
  # when deletion is in "-" columns, deletion present at prev position in mpileup
  check_row_minus <- indel_mpileup %>%
    filter(pos == sample_out[j,]$start - 1)
  check_row_star <- indel_mpileup %>%
    filter(pos == sample_out[j,]$start)
  
  if (!is.na(check_row_star[, "X."]) &
      check_row_star[, "X."] != "") {
    
    if (check_row_star[, "X."] == 0) {
      # no deletion
      del_count <- 0
    } else {
      del_count <- check_row_star[, "X."]
      del_bases <- check_row_star[, "ref"]
      # double check that deleted base matches mutation
      if(sample_out$ref[j] != del_bases) {
        message("check this:")
        message(paste0("deleted base in mpileup: ", del_bases))
        message(paste0("deleted base in mut data: ", sample_out$ref[j]))
        #stop("something weird... deleted bases don't match mutation")
      }
    }
    # write read counts ----------------------------------------
    temp_cov <- check_row_star %>%
      pull(cov)
    sample_out$Distinct.Total.Reads[j] <- temp_cov
    sample_out$Distinct.Mutant.Reads[j] <- del_count
    
  } else if (!is.na(check_row_minus[, "X..1"]) & 
             check_row_minus[, "X..1"] != "") {
    # deletion present at prev position in mpileup
    base_tb <- table(strsplit(check_row_minus[, "X..1"], ",")[[1]])
    if (nrow(base_tb) > 1) {
      message(paste0("j is ", j))
      message(paste(sample_out[j, ], collapse = ","))
      stop("something weird... more than 1 deletion at this position?")
    }
    del_count <- unname(base_tb)
    del_bases <- names(base_tb)
    
    # double check that deleted base matches mutation
    if(sample_out$ref[j] != del_bases) stop("something weird... deleted bases don't match mutation")
    
    # write read counts ----------------------------------------
    # ** need to double check this -- how to count deletion reads? **
    temp_cov <- check_row_minus %>%
      pull(cov)
    sample_out$Distinct.Total.Reads[j] <- temp_cov
    sample_out$Distinct.Mutant.Reads[j] <- del_count
    
  } else {
    
    # deletion might not be present 
    if (!(sample_out[j,]$mutID %in% sample_mutIDs)) {
      # write coverage
      # write read count as 0
      temp_cov <- indel_mpileup %>%
        filter(pos == sample_out[j,]$start) %>%
        pull(cov)
      sample_out$Distinct.Total.Reads[j] <- temp_cov
      sample_out$Distinct.Mutant.Reads[j] <- 0
    } else {
      # deletion might be at individual positions 
      indel_mpileup %>%
        filter(pos == sample_out[j,]$start)
      message("need to write code for deletion")
      message(j)
    }
  }
  return(sample_out)
}

addInsertion <- function(sample_out, j, indel_mpileup) {
  # mutation is insertion ---------------------------------------------------------
  indel_col_ind <-11
  
  # check previous position
  check_row <- indel_mpileup %>%
    filter(pos == sample_out[j,]$start - 1)
  
  # if (is.na(check_row[, indel_col_ind])) {
  #   # no insertion? 
  #   alt_count <- 0
  #   # write read counts
  #   temp_cov <- indel_mpileup %>%
  #     filter(pos == sample_out[j,]$start) %>%
  #     pull(cov)
  #   sample_out$Distinct.Total.Reads[j] <- temp_cov
  #   sample_out$Distinct.Mutant.Reads[j] <- alt_count
  #   
  # } else 
    if (check_row[, indel_col_ind] != "") {
      base_tb <- table(strsplit(check_row[, indel_col_ind], ",")[[1]])
      if (nrow(base_tb) > 1) {
        #message(paste0("j is ", j))
        #message(paste(sample_out[j, ], collapse = ","))
        #stop("something weird... more than 1 insertion at this position?")
        
        # need to choose correct base
        correct_ind <- which(names(base_tb) == sample_out$alt[j])
        alt_count <- unname(base_tb[correct_ind])
        ins_bases <- names(base_tb[correct_ind])
        
      } else {
        alt_count <- unname(base_tb)
        ins_bases <- names(base_tb)
      }
      # make sure inserted bases match reported insertion
      if (sample_out$alt[j] != ins_bases) stop("something weird... inserted bases don't match mutation")
    
    # write read counts
    temp_cov <- indel_mpileup %>%
      filter(pos == sample_out[j,]$start) %>%
      pull(cov)
    sample_out$Distinct.Total.Reads[j] <- temp_cov
    sample_out$Distinct.Mutant.Reads[j] <- alt_count
    
  } else {
    
    # insertion not present?
    if (is.na(sample_out$Distinct.Total.Reads[j]) &
        is.na(sample_out$Distinct.Mutant.Reads[j])) {
      
      # insertion isn't present, grab coverage
      temp_cov <- indel_mpileup %>%
        filter(pos == sample_out[j,]$start) %>%
        pull(cov)
      sample_out$Distinct.Total.Reads[j] <- temp_cov
      sample_out$Distinct.Mutant.Reads[j] <- 0
      
    } else {
      stop("help i dunno why it's not in mpileup")
    }
  }
  
  return(sample_out)
}

addMBCMpileupParsed <- function(mbc_mpileup, sample_out, multibasechange_ind){
  
  mbc_mpileup <- mbc_mpileup %>%
    mutate(chr_start = paste(chrom, pos, sep = "_"))
  
  sample_out <- sample_out %>%
    mutate(chr_start = paste(chr, start, sep = "_"))
  
  for (j in multibasechange_ind) {
    pos_vec <- sample_out$start[j]:sample_out$stop[j]
    #temp_ref <- sample_out$ref[j]
    temp_alt <- sample_out$alt[j]
    
    #message(paste0("querying read counts for mbc: ", paste(sample_out[j, ], collapse = " ")))
    cov_vec <- c()
    alt_vec <- c()
    
    for (mut_base_ind in seq_len(length(pos_vec))) {
      temp_chr_start <- paste(sample_out$chr[j], pos_vec[mut_base_ind], sep = "_")
      #ref_base <- substr(temp_ref, mut_base_ind, mut_base_ind)
      alt_base <- substr(temp_alt, mut_base_ind, mut_base_ind)
      
      temp_mpileup <- mbc_mpileup[match(temp_chr_start, mbc_mpileup$chr_start), ]
      
      cov_vec <- c(cov_vec, temp_mpileup$cov[1])
      alt_vec <- c(alt_vec, temp_mpileup[1, alt_base])
    }
    #message(paste0("cov vec: ", paste(cov_vec, collapse = ", ")))
    #message(paste0("alt vec: ", paste(alt_vec, collapse = ", ")))
    cov_mean <- round(mean(cov_vec), digits = 0)
    alt_mean <- round(mean(alt_vec), digits = 0)
    
    sample_out$Distinct.Mutant.Reads[j] <- alt_mean
    sample_out$Distinct.Total.Reads[j] <- cov_mean
  }
  return(sample_out)
}

readParsedMpileup <- function(file) {
  tb <- read.delim(file,
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
  return(tb)
}
