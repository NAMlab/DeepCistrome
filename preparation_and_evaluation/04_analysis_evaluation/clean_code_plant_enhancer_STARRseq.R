# by Dr. Tobias Jores; 18th November 2025

library(tidyverse)
library(ggsignif)


### function to display sample size on plot ###
give.n <- function(x, y = -Inf, vjust = -1, size = 8/.pt){
  return(c(y = y, vjust = vjust, size = size, label = length(x)))
}


### download missing data ###
## function to download file
download_file <- function(fname, baselink) {
  if (! file.exists(fname)) {
    download.file(
      paste0(baselink, fname),
      fname,
      mode = 'wb'
    )
  }
  
  return(NULL)
}

## list of required files
data_files <- tibble(
  file_name = c('plant-enhancer_variants.tsv', 'PEV_ld_data_main.Rdata', 'validation_refseq_all.fa', 'main_data_validation.Rdata'),
  set = c('enhancer', 'enhancer', 'core promoter', 'core promoter'),
  GitHub = c(
    'https://raw.github.com/tobjores/cooperativity-and-additivity-in-plant-enhancers/main/data/sequence_files/',
    'https://raw.github.com/tobjores/cooperativity-and-additivity-in-plant-enhancers/main/data/RData/',
    'https://raw.github.com/tobjores/Synthetic-Promoter-Designs-Enabled-by-a-Comprehensive-Analysis-of-Plant-Core-Promoters/main/subassembly/refseq/',
    'https://raw.github.com/tobjores/Synthetic-Promoter-Designs-Enabled-by-a-Comprehensive-Analysis-of-Plant-Core-Promoters/main/RData/'
  )
)

## download files
data_files |>
  rowwise() |>
  summarise(
    download_file(file_name, GitHub)
  ) |>
  invisible()


### load enhancer data ###
# sequences of enhancer variants
enhancer_sequences <- read_tsv('plant-enhancer_variants.tsv')

# enhancer strength data for variants
load('PEV_ld_data_main.Rdata') # imports object `data_mean_PEV_ld`


### load core promoter data ###
# sequences of promoter variants
promoter_sequences <- Biostrings::readDNAStringSet('validation_refseq_all.fa') |>
  as.character() |>
  enframe(
    name = 'id',
    value = 'sequence'
  ) |>
  separate_wider_delim(
    id,
    delim = '.',
    names = c('gene', 'motif', 'variant')
  )

# promoter strength data for variants
load('main_data_validation.Rdata') # imports objects `promoter.strength`, `enhancer.effect`, and `light.dependency`


### generate sequences for TFBS prediction ###
## enhancer variants
# define flanking sequences
flanks <- tibble(
  enhancer = rep(c('AB80', 'Cab-1', 'rbcS-E9'), each = 2),
  orientation = rep(c('fwd', 'rev'), times = 3),
  flank_5 = c(
    'GATAAGCTTGATATCGAATTCCACTCTAACTCGCCTCGATC', 'GATAAGCTTGATATCGAATTCCACTCAATGACCTGTCACGG',
    'GATAAGCTTGATATCGAATTCCACTCACATGGGGGATCATG', 'GATAAGCTTGATATCGAATTCCACTCCAGACTGGCTAATGC',
    'GATAAGCTTGATATCGAATTCCACTCGCTCCAGTTCCCAAC', 'GATAAGCTTGATATCGAATTCCACTCTCAGACTGGACGATC'
  ),
  flank_3 = c(
    'CCGTGACAGGTCATTAGGAGCTGGCAAGACCCTTCCTCTA', 'GATCGAGGCGAGTTAAGGAGCTGGCAAGACCCTTCCTCTA',
    'GCATTAGCCAGTCTGAGGAGCTGGCAAGACCCTTCCTCTA', 'CATGATCCCCCATGTAGGAGCTGGCAAGACCCTTCCTCTA',
    'GATCGTCCAGTCTGAAGGAGCTGGCAAGACCCTTCCTCTA', 'GTTGGGAACTGGAGCAGGAGCTGGCAAGACCCTTCCTCTA'
  )
)

# build sequences
TF_pred_seqs_enh <- enhancer_sequences |>
  filter(part == 'B' & type %in% c('substitution', 'WT')) |>
  inner_join(
    flanks,
    by = 'enhancer',
    relationship = 'many-to-many'
  ) |>
  mutate(
    sequence = if_else(
      orientation == 'rev',
      sequence |>
        Biostrings::DNAStringSet() |>
        Biostrings::reverseComplement() |>
        as.character(),
      sequence
    )
  ) |>
  unite(
    'sequence',
    flank_5,
    sequence,
    flank_3,
    sep = ''
  ) |>
  unite(
    'id',
    enhancer,
    part,
    variant,
    orientation
  ) |>
  select(id, sequence) |>
  deframe() |>
  Biostrings::DNAStringSet()


# save sequences to file
TF_pred_seqs_enh |>
  Biostrings::writeXStringSet('plant_enhancer_variants.fa')


## core promoter variants
# build sequences
TF_pred_seqs_pro <- promoter_sequences |>
  filter(motif %in% c('HSF', 'PIF', 'TCP(15)', 'TCP(22)', 'WRKY')) |>
  mutate(
    sequence = paste0('TCCGAGACCACAAGGCGCGCCTAGTGGTCTCCAGGAGCTG', sequence, 'TCTCCGAACTCCGAACCCCAGAACAGAGCAAAGCCTCCTC')
  ) |>
  unite(
    'id',
    gene,
    motif,
    variant
  ) |>
  select(id, sequence) |>
  deframe() |>
  Biostrings::DNAStringSet()

# save sequences to file
TF_pred_seqs_pro |>
  Biostrings::writeXStringSet('core_promoter_variants.fa')


### predict transcription factor binding for sequences in files `plant_enhancer_variants.fa` and `core_promoter_variants.fa` and continue here afterwards ###

### load TFBS predictions and combine with experimental data ###
# enhancers
TFBS_enh <- read_tsv('plant_enhancer_variants_TFBS.tsv')

data_enh <- TFBS_enh |>
  separate_wider_delim(
    SequenceHeader,
    delim = '_',
    names = c('enhancer', 'part', 'variant', 'orientation')
  ) |>
  inner_join(
    data_mean_PEV_ld,
    by = join_by(enhancer, part, variant, orientation)
  )

# promoters
TFBS_pro <- read_tsv('core_promoter_variants_TFBS.tsv')

data_pro <- TFBS_pro |>
  separate_wider_regex(
    SequenceHeader,
    patterns = c(gene = '.*', '_', motif = '.*', '_', variant = '.*')
  ) |>
  inner_join(
    promoter.strength.val,
    by = join_by(gene, motif, variant)
  )


### identify TFBS losses ###
## enhancers
# there are no predicted TFBS losses for AB80 and only few for Cab-1, so we focus on rbcS-E9 here
TFBS_loss_enh <- data_enh |>
  filter(enhancer == 'rbcS-E9') |>
  pivot_longer(
    -names(data_mean_PEV_ld),
    names_to = 'TF',
    values_to = 'TFBS'
  ) |>
  group_by(condition, enhancer, part, orientation, TF) |>
  filter(any(variant == 'WT')) |>
  mutate(
    TFBS_loss = TFBS[variant == 'WT'] > TFBS,
    rel_strength = enrichment - enrichment[variant == 'WT']
  ) |>
  group_by(condition, enhancer, TF) |>
  filter(sum(TFBS_loss) >= 15 & sum(! TFBS_loss) >= 15) |>
  ungroup() |>
  filter(variant != 'WT')


## core promoters
# to complement the tobacco data for the enhancers, we focus on maize data for the core promoters
TFBS_loss_pro <- data_pro |>
  filter(sys == 'proto' & ! enhancer) |>
  pivot_longer(
    -names(promoter.strength.val),
    names_to = 'TF',
    values_to = 'TFBS'
  ) |>
  group_by(gene, TF) |>
  filter(any(variant == 'WT')) |>
  mutate(
    TFBS_loss = TFBS[variant == 'WT'] > TFBS,
    rel_strength = enrichment - enrichment[variant == 'WT']
  ) |>
  group_by(TF) |>
  filter(sum(TFBS_loss) >= 15 & sum(! TFBS_loss) >= 15) |>
  ungroup() |>
  filter(variant != 'WT')


### identify significant differences between sequences with and without TFBS losses ###
## enhancers
signif_enh <- TFBS_loss_enh |>
  group_by(condition, TF) |>
  summarise(
    p_value = wilcox.test(rel_strength[TFBS_loss], rel_strength[! TFBS_loss])$p.value
  ) |>
  ungroup() |>
  mutate(
    p_adj = p_value * n_distinct(TF, condition),
    signif = symnum(
      p_adj,
      cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
      symbols = c('***', '**', '*', 'ns')
    )
  )


## core promoters
signif_pro <- TFBS_loss_pro |>
  group_by(TF) |>
  summarise(
    p_value = wilcox.test(rel_strength[TFBS_loss], rel_strength[! TFBS_loss])$p.value
  ) |>
  ungroup() |>
  mutate(
    p_adj = p_value * n_distinct(TF),
    signif = symnum(
      p_adj,
      cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
      symbols = c('***', '**', '*', 'ns')
    )
  )


### plot results ###
## enhancers
TFBS_loss_enh |>
  ggplot(aes(x = TF, y = rel_strength, fill = TFBS_loss)) +
  facet_grid(
    cols = vars(condition),
    labeller = labeller(
      condition = c(
        'dark' = 'tobacco leaves in dark',
        'light' = 'tobacco leaves in light'
      )
    )
  ) +
  geom_hline(
    yintercept = 0
  ) +
  geom_boxplot() +
  geom_signif(
    xmin = as.numeric(as.ordered(signif_enh$TF)) - 0.1875,
    xmax = as.numeric(as.ordered(signif_enh$TF)) + 0.1875,
    y_position = 1.2 * max(TFBS_loss_enh$rel_strength),
    tip_length = 0,
    annotation = signif_enh$signif,
    size = 1,
    textsize = 20/.pt
  ) +
  stat_summary(
    fun.data = give.n,
    geom = 'text',
    size = 20/.pt,
    position = position_dodge(0.75)
  ) +
  scale_x_discrete(
    name = 'tf family'
  ) +
  scale_y_continuous(
    name = expression(log[2]('enhancer strength, rel. to WT')),
    expand = expansion(mult = c(0.2, 0.1))
  ) +
  scale_fill_discrete(
    name = NULL,
    labels = c('TRUE' = 'loss of TFBS', 'FALSE' = 'no TFBS loss')
  ) +
  theme_bw(
    base_size = 24
  ) +
  theme(
    legend.position = 'top',
    axis.text = element_text(
      color = 'black'
    )
  )

ggsave('TFBS_loss_enhancers.pdf')


## core promters
TFBS_loss_pro |>
  ggplot(aes(x = TF, y = rel_strength, fill = TFBS_loss)) +
  facet_grid(
    cols = vars(sys),
    labeller = labeller(
      sys = c('proto' = 'maize protoplasts')
    )
  ) +
  geom_hline(
    yintercept = 0
  ) +
  geom_boxplot() +
  geom_signif(
    xmin = as.numeric(as.ordered(signif_pro$TF)) - 0.1875,
    xmax = as.numeric(as.ordered(signif_pro$TF)) + 0.1875,
    y_position = 1.2 * max(TFBS_loss_pro$rel_strength),
    tip_length = 0,
    annotation = signif_pro$signif,
    size = 1,
    textsize = 20/.pt
  ) +
  stat_summary(
    fun.data = give.n,
    geom = 'text',
    size = 20/.pt,
    position = position_dodge(0.75)
  ) +
  scale_x_discrete(
    name = 'tf family'
  ) +
  scale_y_continuous(
    name = expression(log[2]('promoter strength, rel. to WT')),
    expand = expansion(mult = c(0.2, 0.1))
  ) +
  scale_fill_discrete(
    name = NULL,
    labels = c('TRUE' = 'loss of TFBS', 'FALSE' = 'no TFBS loss')
  ) +
  theme_bw(
    base_size = 24
  ) +
  theme(
    legend.position = 'top',
    axis.text = element_text(
      color = 'black'
    )
  )

ggsave('TFBS_loss_promoters.pdf')
