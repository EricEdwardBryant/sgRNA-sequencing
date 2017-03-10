This is a sample analysis of a pooled CRISPRa NGS experiment. Our starting materials are:

1.  Forward and reverse sgRNA barcode reference files in fasta format (e.g. `CRISPRa_v2_human.trim_1_29_forward.fa`)
2.  A corresponding reference annotation file (e.g. `CRISPRa_v2_human`)
3.  Sequence files (`*.fastq.gz`)
4.  Code in this repository located in the `R` directory

If this is your first time running this analysis on your computer, you should first install a few items. Namely:

1.  Install [R](https://cloud.r-project.org)
2.  Install [Rstudio](https://www.rstudio.com/products/rstudio/download)
3.  Open Rstudio and run the following commands to install [Bioconductor](http://bioconductor.org) and the [tidyverse](http://tidyverse.org)

``` r
# Install the latest version of Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()

# Install the tidyverse
install.packages('tidyverse')
```

Phew! I'm exhausted, but we shall persist!

To begin our analysis, we will first load up all of the software requirements using the following commands:

``` r
library(tidyverse)
library(stringr)
requireNamespace('Biostrings')
# The quoted strings are file paths to some R scripts available in this GitHub 
# code repository that contain some helpful functions. Source them to reap their benefits!
source('R/functions/tidyseq.R')
source('R/functions/read-sequences.R')
source('R/functions/count-sgRNA.R')
source('R/functions/convert-sgRNA-reference-to-CSV.R')
source('R/functions/annotate-sequence-samples.R')
```

Now that we're all set up, the fun begins! I first convert the reference sequences to a CSV file. Why? I just find this tabular format easier to work with. To convert the reference sequences I use the handy dandy `conver_sgRNA_reference_to_CSV` function. It will take a few seconds because, well, these reference files are big-ish, and my code is slow.

``` r
convert_sgRNA_reference_to_CSV(
  # What directory should the resulting files be written to?
  dir = 'data/guide-reference', 
  # What will the name of this version of the sgRNA library be?
  vsn_sgRNA = 'CRISPRa-v2-human-29-top5-h1-h2-h3-h4-h5-h6-h7', 
  # What is the path to the forward oriented fasta reference file?
  fwd_fasta = 'data/guide-reference/CRISPRa_v2_human.trim_1_29_forward.fa',
  # What is the path to the reverse oriented fasta reference file?
  rev_fasta = 'data/guide-reference/CRISPRa_v2_human.trim_1_29_reverse.fa',
  # What is the path to the corresponding reference annotation file?
  lib_tbl   = 'data/guide-reference/CRISPRa_v2_human_librarytable.txt',
  # What sublibraries did you use to make your collection?
  sublibs   = c(
    'h1_top5'  = 'Kinases, phosphatases, and drug targets',
    'h2_top5'  = 'Cancer and apoptosis',
    'h3_top5'  = 'Stress proteostasis',
    'h4_top5'  = '?',
    'h5_top5'  = 'Gene expression',
    'h6_top5'  = '??',
    'h7_top5'  = 'Unassigned'#,
    #  'h1_supp5' = 'Kinases, phosphatases, and drug targets',
    #  'h2_supp5' = 'Cancer and apoptosis',
    #  'h3_supp5' = 'Stress proteostasis',
    #  'h4_supp5' = '?',
    #  'h5_supp5' = 'Gene expression',
    #  'h6_supp5' = '??',
    #  'h7_supp5' = 'Unassigned'
  )
)
```

If you have succeeded, you should have a custom reference table for your sgRNA library. Hooray! If you're curious what these files look like, you can read them in and check there contents with the following commands. Note that three files have been written. One has all possiple guides, one is your unambiguously mapped sub-collection, and another is a set of ambiguously mapped guides (problems to be delt with at a later time?).

``` r
sgRNA <- read_csv('data/guide-reference/CRISPRa-v2-human-29-top5-h1-h2-h3-h4-h5-h6-h7.csv')
sgRNA_ambiguous <- read_csv('data/guide-reference/CRISPRa-v2-human-29-top5-h1-h2-h3-h4-h5-h6-h7-ambiguous.csv')
# Let's have a look at each table. Also try clicking on the data set in the "Environment" panel in Rstudio
sgRNA
```

    ## # A tibble: 98,177 × 8
    ##     gene             sequence                          fwd
    ##    <chr>                <chr>                        <chr>
    ## 1   A1BG GAAGACAGGGAAGATGAAGC AAGACAGGGAAGATGAAGCGTTTAAGAG
    ## 2   A1BG GAGCAGCTCCAGGTAGAGTG AGCAGCTCCAGGTAGAGTGGTTTAAGAG
    ## 3   A1BG GCAGGGCCCCCATGGGGTCA CAGGGCCCCCATGGGGTCAGTTTAAGAG
    ## 4   A1BG GCGCGCCTGCGCCTCAGCCC CGCGCCTGCGCCTCAGCCCGTTTAAGAG
    ## 5   A1BG GCTGCCCCCGCCTTCCGGGA CTGCCCCCGCCTTCCGGGAGTTTAAGAG
    ## 6   A1BG GGGCAAAGGGTGAACTTCTG GGCAAAGGGTGAACTTCTGGTTTAAGAG
    ## 7   A1BG GGGGACACTCACGTGTGGCG GGGACACTCACGTGTGGCGGTTTAAGAG
    ## 8   A1BG GGTGCGGGGACACTCACGTG GTGCGGGGACACTCACGTGGTTTAAGAG
    ## 9   A1BG GTGGGCGCAGAGGGCTCCTC TGGGCGCAGAGGGCTCCTCGTTTAAGAG
    ## 10  A1BG GTGTCCTTCCCGGAAGGCGG TGTCCTTCCCGGAAGGCGGGTTTAAGAG
    ## # ... with 98,167 more rows, and 5 more variables: rev <chr>, id <chr>,
    ## #   sublibrary <chr>, transcripts <chr>, sublibrary_name <chr>

``` r
sgRNA_ambiguous
```

    ## # A tibble: 6,307 × 8
    ##                              id sublibrary
    ##                           <chr>      <chr>
    ## 1     NUSAP1_-_41624953.23-P1P2    h2_top5
    ## 2       OIP5_-_41624953.23-P1P2    h5_top5
    ## 3       DCAF8_+_160255023.23-P1    h3_top5
    ## 4     PEX19_+_160255023.23-P1P2    h3_top5
    ## 5     LRRC49_+_71184970.23-P1P2    h4_top5
    ## 6     THAP10_+_71184970.23-P1P2    h5_top5
    ## 7        ACIN1_+_23564543.23-P1    h3_top5
    ## 8  C14orf119_+_23564543.23-P1P2    h4_top5
    ## 9        PSMB9_+_32821773.23-P1    h3_top5
    ## 10      TAP1_+_32821773.23-P1P2    h1_top5
    ## # ... with 6,297 more rows, and 6 more variables: sublibrary_name <chr>,
    ## #   gene <chr>, transcripts <chr>, sequence <chr>, fwd <chr>, rev <chr>

Now we will annotate your sequence samples using the handy dandy `annotate_sequence_samples` function. Note that this function is a little delicate as it uses your file names and time stampts to figure some things out automagically, so if your file name format changes this function may need to be fixed. Of course, you could always just make the annotation file by hand!

``` r
annotate_sequence_samples(
  # Where are your sequence files?
  dir = 'data/fastq',
  # Where do you want to save this file and what should we name it?
  name = 'data/sample-annotation.csv'
)
```

Let's have a look at that file!

``` r
fastq <- read_csv('data/sample-annotation.csv')
fastq
```

    ## # A tibble: 8 × 8
    ##                   sample_id       date   name sample  lane size_mb
    ##                       <chr>     <date>  <chr>  <int> <int>   <dbl>
    ## 1 2017-02-24-hCRAv2_S12_001 2017-02-24 hCRAv2     12     1     198
    ## 2 2017-02-24-hCRAv2_S12_002 2017-02-24 hCRAv2     12     2     179
    ## 3 2017-02-24-hCRAv2_S12_003 2017-02-24 hCRAv2     12     3     203
    ## 4 2017-02-24-hCRAv2_S12_004 2017-02-24 hCRAv2     12     4     164
    ## 5 2017-02-24-hCRAv2_S13_001 2017-02-24 hCRAv2     13     1     209
    ## 6 2017-02-24-hCRAv2_S13_002 2017-02-24 hCRAv2     13     2     207
    ## 7 2017-02-24-hCRAv2_S13_003 2017-02-24 hCRAv2     13     3     228
    ## 8 2017-02-24-hCRAv2_S13_004 2017-02-24 hCRAv2     13     4     208
    ## # ... with 2 more variables: file_name <chr>, path <chr>

Time for the heavy lifting! We will count barcodes using the the `count_sgRNA` function. This will take a bit of time depending on how big your sequence files are. Be patient, grab a coffee!

``` r
# Trimming is key here as the counting is based on perfect matches to your sgRNA reference table
count_sgRNA(
  # This is your annotation table which contains paths and sample information for all
  # of your fastq sequence files
  fastq, 
  # This is your sgRNA library reference table that contains the expected sequences
  sgRNA,
  # The following settings will keep bases 37 to 64 of each read found in your fastq files
  # and match this sequence to your sgRNA reference table
  trim_before = 37, 
  trim_after = 64,
  # Where should the counts be saved (this directory should already exist)
  dir = 'data/counts'
)
```
