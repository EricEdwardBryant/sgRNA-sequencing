# Read FASTA or FASTQ files into a dataframe with columns 'id' and 'sequence'
read_sequences <- function(path) {
  seqs <- Biostrings::readDNAStringSet(path) %>% as.character
  data_frame(id = names(seqs), sequence = unname(seqs))
}
