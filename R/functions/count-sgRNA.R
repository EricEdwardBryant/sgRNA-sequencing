count_sgRNA <- function(fastq, sgRNA, trim_before, trim_after, dir) {
  message('Counting sgRNA barcodes in ', nrow(fastq), ' fastq files\n',
          'Writing results to ', dir)
  if (!dir.exists(dir)) dir.create(dir)
  p <- progress_estimated(nrow(fastq) * 2)

  lapply(1:nrow(fastq), function(i) {
    reads <- read_trimmed_sequences(fastq$path[i], trim_before, trim_after)
    sample_i <- fastq$sample_id[i]

    p$tick()$print()
    counts <- count(reads, sequence)

    # Map reads to forward and reverse orientation
    mapped   <- sgRNA %>%
      left_join(counts %>% rename(n_fwd = n, fwd = sequence), by = 'fwd') %>%
      left_join(counts %>% rename(n_rev = n, rev = sequence), by = 'rev') %>%
      mutate(
        n_fwd = ifelse(is.na(n_fwd), 0, n_fwd),
        n_rev = ifelse(is.na(n_rev), 0, n_rev),
        n = n_fwd + n_rev,
        sample_id = sample_i
      ) %>%
      select(sample_id, id, sublibrary, sublibrary_name, gene, n, n_fwd, n_rev, sequence, fwd, rev)

    unmapped <- counts %>%
      filter(
        !(sequence %in% sgRNA$fwd),
        !(sequence %in% sgRNA$rev)
      ) %>%
      mutate(sample_id = sample_i)

    write_csv(mapped,   str_c(dir, '/', sample_i, '-mapped.csv'))
    write_csv(unmapped, str_c(dir, '/', sample_i, '-unmapped.csv'))

    p$tick()$print()
    return(invisible())
  })
  return(invisible(dir))
}
