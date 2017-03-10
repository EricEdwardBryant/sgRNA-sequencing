# This script takes Fasta formatted sgRNA reference libraries in forward and
# reverse orientation, combines them, parses the IDs for gene strand and
# coordinate. IDs are assumed to be in the format: "gene_[+-]_coordinate"
#
# It returns two tables:
# 1. ambiguously mapping sgRNAs (i.e. those that map to multiple genes)
# 2. Non-ambiguously mapping sgRNAs - if multiple strands and coordinates for
#    a given gene/sgRNA, they are collapsed and separated with " | ".
#    Thus only one unique sgRNA sequence per row should be represented in this
#    table.
convert_sgRNA_reference_to_CSV <- function(dir, vsn_sgRNA, fwd_fasta, rev_fasta, lib_tbl, sublibs) {

  barcodes <-
    list(
      fwd = read_sequences(fwd_fasta) %>% mutate(sequence = toupper(sequence)),
      rev = read_sequences(rev_fasta) %>% mutate(sequence = toupper(sequence))
    ) %>%
    bind_rows(.id = 'orientation') %>%
    spread(orientation, sequence) %>%
    mutate(id = str_replace(id, '^>', ''))

  library_full <-
    read_tsv(
      lib_tbl,
      col_types = 'ccccc',
      col_names = c('id', 'sublibrary', 'gene', 'transcripts', 'sequence'),
      skip = 1
    ) %>%
    left_join(barcodes, by = 'id') %>%
    mutate(
      sequence = toupper(sequence),
      sublibrary_name = sublibs[sublibrary]
    )

  library_sub <-
    library_full %>%
    filter(sublibrary %in% names(sublibs))

  library_sub_ambiguous <-
    library_sub %>%
    filter(duplicated(sequence) | duplicated(sequence, fromLast = T)) %>%
    group_by(sequence) %>%
    filter(length(unique(gene)) > 1L) %>%
    arrange(sequence) %>%
    select(id, sublibrary, sublibrary_name, gene, transcripts, sequence, fwd, rev)

  # Barcodes are allowed to map to more than one library/strand/coordinate so
  # long as they are for a single gene
  library_sub_distinct <-
    library_sub %>%
    anti_join(library_sub_ambiguous, by = 'sequence') %>%
    group_by(gene, sequence, fwd, rev) %>%
    summarise_all(funs(ifelse(length(unique(.)) > 1, str_c(., collapse = ' | '), .)))

  write_csv(library_full, str_c(dir, '/', 'CRISPRa-v2-human-29.csv'))
  write_csv(library_sub_ambiguous, str_c(dir, '/', vsn_sgRNA, '-ambiguous.csv'))
  write_csv(library_sub_distinct,  str_c(dir, '/', vsn_sgRNA, '.csv'))
  return(invisible(dir))
}
