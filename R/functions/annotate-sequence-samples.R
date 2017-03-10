
annotate_sequence_samples <- function(dir, name) {
  fastq_annotation <-
    data_frame(
      path = list.files(dir, 'fastq.gz', recursive = T, full.names = T),
      name = basename(path) %>% gsub('[.]fastq[.]gz$', '', .)
    ) %>%
    separate(name, into = c('name', 'sample', 'lane'), sep = '_') %>%
    mutate_at(vars(sample:lane), funs(as.integer(gsub('[^0-9]', '', .)))) %>%
    mutate(
      date = as.Date(file.info(path)$mtime),
      file_name = basename(path),
      size_mb = (file.info(path)$size / 1000000) %>% round(0),
      sample_id = str_c(date, '-', str_extract(file_name, '.*(?=(.fastq.gz))'))
    ) %>%
    select(sample_id, date, name:lane, size_mb, file_name, path)

  write_csv(fastq_annotation, name)
}
