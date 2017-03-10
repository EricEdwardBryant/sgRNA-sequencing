#' Read fastq file
#'
#' Reads sequence data from fastq files
#'
#' @param path String. A path to a fastq sequence file. See details.
#' @param trim_start Integer. Which 5' position in sequence to trim before
#' (e.g. 2 will trim all bases before the 2nd base).
#' Defaults to \code{1L} (i.e. no 5' trimming).
#' @param trim_end Integer. Which 3' position in sequence to trim after (e.g.
#' 30 will trim all bases after the 30th base). Defaults
#' to \code{-1L} (i.e. no 3' trimming).
#'
#' @details The fastq file can be either uncompressed (*.fq, *.fastq), or
#' compressed formats (e.g. *.fq.gz, *.fastq.gz). Each read is expected to have
#' 4 lines:
#'
#' \enumerate{
#'   \item Sequence ID
#'   \item Sequence
#'   \item (ignored)
#'   \item Phred33 encoded quality scores.
#' }
#'
#' @importFrom stringr str_sub str_extract
#' @importFrom readr read_lines
#' @importFrom tibble data_frame as_data_frame
#' @export

read_trimmed_sequences <- function(path, trim_start = 1L, trim_end = -1L) {
  # message('Reading: ', basename(path), ' ...')
  lines  <- read_lines(path, progress = F)
  class(lines) <- class_from_extension(path)
  as_data_frame(lines, trim_start, trim_end)
}

as_data_frame.fastq <- function(lines, trim_start, trim_end) {
  index  <- seq(1, length(lines), by = 4)
  data_frame(
    id       = lines[index],
    sequence = str_sub(lines[index + 1], start = trim_start, end = trim_end),
    quality  = str_sub(lines[index + 3], start = trim_start, end = trim_end)
  )
}

as_data_frame.fasta <- function(lines, trim_start, trim_end) {
  lines    <- str_trim(lines)
  id_index <- which(str_detect(lines, '^>'))
  data_frame(
    id       = lines[id_index],
    sequence = map2_chr(
      .x = id_index,
      .y = c(id_index[-1], length(lines) + 1),
      ~(lines[(.x + 1):(.y - 1)]) %>% str_c(collapse = '')
    )
  ) %>%
  mutate(sequence = str_sub(sequence, start = trim_start, end = trim_end))
}

#' Write fastq
write_sequence.data.frame <- function(sequence_df, file) {
  n <- nrow(reads)
  plus <- rep('+', n) # Create a vector of + to interleave at the third line
  interleave <- rep(1:n, each = 4) + seq(0, 4 - 1) * n  # 4 is number of vectors
  lines <- with(reads, c(id, sequence, plus, quality))[interleave]
  writeLines(lines, con = file)
}

add_quality_metrics <- function(reads) {
  reads %>%
    mutate(
      qs_avg   = quality %>% map_dbl(~mean(utf8ToInt(.) - 33)), # Illumina 1.8+ Phred+33,  raw reads typically (0, 41)
      length   = nchar(sequence)
    )
}

class_from_extension <- function(path) {
  extension <- str_extract(path, 'f(a|q|sa|sq|asta|astq)?(\\.gz)?$')
  switch(
    extension,
    fa       = 'fasta',
    fa.gz    = 'fasta',
    fsa      = 'fasta',
    fsa.gz   = 'fasta',
    fasta    = 'fasta',
    fasta.gz = 'fasta',
    fq       = 'fastq',
    fq.gz    = 'fastq',
    fsq      = 'fastq',
    fsq.gz   = 'fastq',
    fastq    = 'fastq',
    fastq.gz = 'fastq'
  )
}
