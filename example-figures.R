library(tidyverse)

mapped   <- list.files('data/counts', '-mapped', full.names = T)   %>% map_df(read_csv, col_types = cols())
unmapped <- list.files('data/counts', '-unmapped', full.names = T) %>% map_df(read_csv, col_types = cols())

gg_unmapped_reads <-
  unmapped %>%
  group_by(sample_id) %>%
  summarise(reads = sum(n)) %>%
  ggplot(aes(x = sample_id, y = reads, fill = sample_id)) +
  geom_col(alpha = 0.5, color = 'black') +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = 'none') +
  coord_flip(ylim = c(0, 5.5e6)) +
  theme_bw() +
  labs(x = 'Sample ID', y = 'Unmapped reads') +
  theme(aspect.ratio = 1)
gg_unmapped_reads

ggsave('figures/unmapped-reads.pdf', gg_unmapped_reads, width = 5, height = 5)

gg_sublibrary_coverage <-
  mapped %>%
  group_by(sublibrary) %>%
  summarise(reads = sum(n)) %>%
  ggplot(aes(x = sublibrary, y = reads, fill = sublibrary)) +
  geom_col(alpha = 0.5, color = 'black') +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = 'none') +
  coord_flip(ylim = c(0, 5.5e6)) +
  theme_bw() +
  labs(x = 'Sub-library', y = 'Mapped reads') +
  theme(aspect.ratio = 1)
gg_sublibrary_coverage

ggsave('figures/sublibrary-coverage.pdf', gg_sublibrary_coverage, width = 5, height = 5)

gg_guide_coverage <-
  mapped %>%
  ggplot(aes(x = n)) +
  stat_ecdf() +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = 'Reads mapped per sgRNA', y = 'Cumulative distribution') +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(aspect.ratio = 1)

ggsave('figures/guide-coverage.pdf', gg_guide_coverage, width = 5, height = 5)

