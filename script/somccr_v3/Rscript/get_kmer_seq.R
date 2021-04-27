################################################
##### Function to get sequences of k-mers ######
################################################
## Load in packages
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(pbmcapply)

get_kmer_freq <- function(file, k) {
  ## Read in the data
  df <- fread(file)
  
  ## Convert data.table rows to individual list elements
  #df_list <- split(df, seq(nrow(df)))

  ##########################################################################################
  ##### Iterate through each variant to get the reference and alternate k-mer sequence #####
  ##########################################################################################
  kmer_df <- rbindlist(apply(df, 1, function(row, k) {
    #browser()
    ## Get information for seuqnece extraction
    chromosome <- row["chromosome"] %>% as.character()
    start <- as.numeric(row["stop"]) - ((k-1)/2) ## The stop coordinate is the 1-base position of the variant
    stop <- as.numeric(row["stop"]) + ((k-1)/2)
    
    ## Get reference kmer sequence
    ref_kmer_seq <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, chromosome, start, stop))
    
    ## Get the alternate allele 
    alt_allele <- row["alt"] %>% as.character()
    
    ## Get the alternate kmer sequence 
    pos <- ((k-1)/2) + 1
    alt_kmer_seq <- ref_kmer_seq
    substr(x = alt_kmer_seq, start = pos, stop = pos) <- alt_allele
    
    ## Make temp dataset of the kmers 
    kmer_df <- data.table("ref_kmer" = ref_kmer_seq, 
                          "alt_kmer" = alt_kmer_seq,
                          "kmer_id" = paste0(ref_kmer_seq, " > ", alt_kmer_seq))
    
    ## Combine the data
    final <- cbind(t(row), kmer_df)
    
    return(final)
    
  }, k=k))

  ##########################################
  ##### Get the frequency of each kmer #####
  ##########################################
  kmer_freq_table <- data.table(table(kmer_df$kmer_id)) %>% `colnames<-`(c("kmer_id", "count"))
  kmer_freq_table$total <- sum(kmer_freq_table$count)
  kmer_freq_table$freq <- round(kmer_freq_table$count/kmer_freq_table$total, digits = 3)
  
  ## Return dataset and kmer freq
  final <- list(kmer_df, kmer_freq_table)
  names(final) <- c("mut_df", "kmer_freq_df")
  return(final)
}

#############################################
##### Function to map k-mer frequencies #####
#############################################
map_kmer_freq <- function(df) {
  ## Load in the data
  mut_df <- df$mut_df
  kmer_df <- df$kmer_freq_df
  
  ## Get the unique kmer names
  uniq_kmer <- kmer_df$kmer_id %>% unique()
  
  ##############################################################
  ##### Iterate through each kmer and report its frequency #####
  ##############################################################
  kmer_df_freq <- rbindlist(lapply(uniq_kmer, function(kmer, mut_df, kmer_df){
    ## Subset the original dataset
    uniq_kmer_df <- mut_df[kmer_id == kmer]
    
    ## Get the frequency
    kmer_freq <- kmer_df$freq[which(kmer_df$kmer_id == kmer)]
    
    ## Add frequency to the original dataset
    uniq_kmer_df$freq <- kmer_freq
    
    return(uniq_kmer_df)
  }, kmer_df=kmer_df, mut_df=mut_df))
  
  ## Order the final dataset
  kmer_df_freq_sorted <- kmer_df_freq[order(chromosome, start, stop, gene, cancer_type, variant_class)]
  head(kmer_df_freq_sorted)

  return(kmer_df_freq_sorted)
}


file <- "~/Google Drive/Quinlan Lab - PhD/Projects/somatic_ccr/data/output/mc3_pcawg_v3/mc3_pcawg.sorted.filtered.snv.CDS.knownCDSlength.bed"
output_file <- "~/Google Drive/Quinlan Lab - PhD/Projects/somatic_ccr/data/output/mc3_pcawg_v3/mc3_pcawg.sorted.filtered.snv.CDS.knownCDSlength.kmer.bed"
k <- 3
df <- get_kmer_freq(file = file, k = k)
final <- map_kmer_freq(df = df)

## Write the final table
write.table(x = final, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)


