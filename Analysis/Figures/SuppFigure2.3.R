
rm(list = ls(all = TRUE))
# SETUP ===========================================================================================
# *************************************************************************************************

library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Setup -------------------------------------------------------------------
### Output directory --------------------------------------------------------
outdir = "Analysis/Figures/"
out_plots="plots"
out_stats="stats"

input_dir = "Analysis/tables/"

### create directories ------------------------------------------------------
#create paths, change
out_plots = file.path(outdir, out_plots)
out_stats = file.path(outdir, out_stats)
#create directories, if do not exist
dir.create(out_plots, recursive = T)
dir.create(out_stats, recursive = T)


# Load intermediate -------------------------------------------------------


# 3'UTR to gene - unmerged to account for internal isoforms
utr3_to_genes = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithOppositeOrientationRepeat.sorted.csv", stringsAsFactors = F)

tandem_utr3_to_genes = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithSameOrientationRepeat.sorted.csv", stringsAsFactors = F)


# alu inverted regions 
alus = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RepeatsInOppositeOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

alus_merged = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RepeatsInOppositeOrientationAtRegions.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# alu tandem regions 
tandem_alus = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RepeatsInSameOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

tandem_alus_merged = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RepeatsInSameOrientationAtRegions.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# 3'UTR regions
utr3s = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithOppositeOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

tandem_utr3s = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithSameOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

# intersection results
intersection_res = fread("analysis/mm10.AluB1AndB2_3pUTR_minLen120_17022021.intersection_results.csv", stringsAsFactors = F) %>%
  mutate("RefSeqID" = str_remove(region_name, "_utr3_.*")) %>%
  distinct() %>%
  mutate(repeat_score = case_match(repeat_score, "Alu" ~ "B1", 
                                   .default = repeat_score))

# get only the relevant intersections
intersection_res_cytoindex_regions = intersection_res %>%
  filter(region %in% utr3s, repeat_cropped %in% alus) %>%
  inner_join(utr3_to_genes, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
  rename("UTR3Region" = region,
         "UTR3RegionLength" = region_length,
         "AluElement" = `repeat`,
         "Family" = repeat_score,
         "AluElementCroppedToUTR3" = repeat_cropped,
         "AluElementOrientation" = repeat_strand,
         "AluElementLength" = repeat_length,
         "GeneStrand" = Strand) %>%
  # inner_join(alu_families) %>%
  # inner_join(alu_divergence) %>%
  mutate("AluElementCroppedToUTR3Length" = as.numeric(AluElementCroppedToUTR3 %>% str_extract(pattern = "-[0-9]+") %>%
                                                        str_extract( pattern = "[0-9]+")) - 
           as.numeric(AluElementCroppedToUTR3 %>% 
                        str_extract(pattern = ":[0-9]+-") %>%
                        str_extract(pattern = "[0-9]+")),) %>%
  select(-contains("_"), -RegionLength)

# get only the relevant intersections
tandem_intersection_res_cytoindex_regions = intersection_res %>%
  filter(region %in% tandem_utr3s, repeat_cropped %in% tandem_alus) %>%
  inner_join(tandem_utr3_to_genes, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
  rename("UTR3Region" = region,
         "UTR3RegionLength" = region_length,
         "AluElement" = `repeat`,
         "Family" = repeat_score,
         "AluElementCroppedToUTR3" = repeat_cropped,
         "AluElementOrientation" = repeat_strand,
         "AluElementLength" = repeat_length,
         "GeneStrand" = Strand) %>%
  # inner_join(alu_families) %>%
  # inner_join(alu_divergence) %>%
  select(-contains("_"), -RegionLength)

# no refseq id and gene \ gene length - just different 3'UTRs per GeneSymbol
intersection_res_cytoindex_regions_unique = intersection_res_cytoindex_regions %>%
  select(-RefSeqID, -Gene, -GeneLength) %>%
  distinct()

tandem_intersection_res_cytoindex_regions_unique = tandem_intersection_res_cytoindex_regions %>%
  select(-RefSeqID, -Gene, -GeneLength, ) %>%
  distinct()


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")
chosen_colors_trio = c("white", index_colors)

# Supp. Fig. 3B -----------------------------------------------------------------

### count clusters ----
cluster_count_per_gene = intersection_res_cytoindex_regions_unique %>% 
  group_by(GeneSymbol, UTR3Region, Family) %>% 
  tally %>% 
  group_by(GeneSymbol, Family)%>%
  # take largest cluster per gene
  slice_max(order_by = n, with_ties = F) %>%
  mutate(Group = "Inverted") %>%
  bind_rows(tandem_intersection_res_cytoindex_regions_unique %>% 
              group_by(GeneSymbol, UTR3Region, Family) %>% 
              tally %>% 
              group_by(GeneSymbol, Family)%>%
              # take largest cluster per gene
              slice_max(order_by = n, with_ties = F) %>%
              mutate(Group = "Tandem")) %>%
  group_by(n, Group, Family) %>%
  tally(name = "Number of Genes") %>%
  ungroup()


suppFig3b = cluster_count_per_gene %>%
  filter(n > 1) %>%
  complete(Group, Family, n = 2:6, fill = list(`Number of Genes` = 0)) %>%
  # filter(n>1) %>%
  ggplot(aes(x = n, y = `Number of Genes`, fill = Group)) +
  # geom_bar(position = position_dodge()) +
  geom_col(position = position_dodge(), color= "black") +
  facet_grid(.~Family, scale = "free", space ="free") +
  theme_custom() +
  scale_x_continuous(breaks = c(2:6)) +
  # geom_text(stat='count', aes(label= after_stat(count)), position = position_dodge(width = 0.7), vjust=-1) +
  # geom_text(aes(label= `Number of Genes`,
  #               color = Group), position = position_dodge(width = 0.7), vjust=-1) +
  geom_text(aes(label= if_else(`Number of Genes` == 0, "", as.character(`Number of Genes`)),
                color = Group), position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5) +
  # scale_fill_gradient(low = case_control_colors[1], high = case_control_colors[2]) +
  scale_color_manual(values = index_colors) +
  scale_fill_manual(values = index_colors) +
  ggtitle("Inverted cluster size per gene",
          subtitle = "Largest cluster per gene (GeneSymbol)") +
  xlab("Number of elements in cluster per gene") +
  ylab("Number of genes") +
  guides(color = "none") +
  theme(legend.title = element_blank())

### expected inverted vs. tandem ratio -----
expected_inverted_to_tandem_ratio = data.frame(n = 2:6, 
                                               "Expected Ratio" = 2^(2:6 - 1) - 1,
                                               check.names = F) 

expected_inverted_relative_to_tandem = cluster_count_per_gene %>%
  filter(n > 1, Group == "Tandem") %>%
  complete(Group, Family, n = 2:6, fill = list(`Number of Genes` = 0)) %>%
  inner_join(expected_inverted_to_tandem_ratio) %>%
  mutate(Group = "Expected",
         `Number of Genes` = `Number of Genes` * `Expected Ratio`) %>%
  select(-`Expected Ratio`)


suppFig3b_inset_b1 = cluster_count_per_gene %>%
  filter(n > 1, Group == "Inverted") %>%
  complete(Group, n = 2:6, fill = list(`Number of Genes` = 0)) %>%
  mutate(Group = "Observed") %>%
  bind_rows(expected_inverted_relative_to_tandem) %>%
  filter(Family=="B1") %>%
  ggplot(aes(x = n, y = `Number of Genes`, fill = Group)) +
  # geom_bar(position = position_dodge()) +
  geom_col(position = position_dodge(), color= "black") +
  theme_custom(title_size = 12, text_size = 12) +
  scale_x_continuous(breaks = c(2:7), expand = c(0.02, 0.01)) +
  scale_y_continuous(breaks = seq(0, 550, 100)) +
  expand_limits(y = 150) +
  geom_text(aes(label= as.character(`Number of Genes`),
                color = Group), position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5, size = 3) +
  # scale_fill_gradient(low = case_control_colors[1], high = case_control_colors[2]) +
  scale_color_manual(values = c("black", index_colors[1])) +
  scale_fill_manual(values = c("white", index_colors[1])) +
  ggtitle("Expected vs. observed inverted genes",
          subtitle = "Relative to observed tandem genes") +
  xlab("# of elements in cluster per gene") +
  ylab("# of genes") +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA),
        legend.title = element_blank())


suppFig3b_inset_b2 = cluster_count_per_gene %>%
  filter(n > 1, Group == "Inverted") %>%
  complete(Group, n = 2:6, fill = list(`Number of Genes` = 0)) %>%
  mutate(Group = "Observed") %>%
  bind_rows(expected_inverted_relative_to_tandem) %>%
  filter(Family=="B2", n<=3) %>%
  ggplot(aes(x = n, y = `Number of Genes`, fill = Group)) +
  # geom_bar(position = position_dodge()) +
  geom_col(position = position_dodge(), color= "black") +
  theme_custom(title_size = 12, text_size = 12) +
  scale_x_continuous(breaks = c(2:7), expand = c(0.02, 0.01)) +
  scale_y_continuous(breaks = seq(0, 550, 100)) +
  expand_limits(y = 150) +
  geom_text(aes(label= as.character(`Number of Genes`),
                color = Group), position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5, size = 3) +
  # scale_fill_gradient(low = case_control_colors[1], high = case_control_colors[2]) +
  scale_color_manual(values = c("black", index_colors[1])) +
  scale_fill_manual(values = c("white", index_colors[1])) +
  ggtitle("Expected vs. observed inverted genes",
          subtitle = "Relative to observed tandem genes") +
  xlab("# of elements in cluster per gene") +
  ylab("# of genes") +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA),
        legend.title = element_blank())

# Supp. Fig. 3C -----------------------------------------------------------------

# # Chi squared
chiseq_pvals = cluster_count_per_gene %>%
  pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
  # compute expected tandem elements
  inner_join(expected_inverted_to_tandem_ratio) %>%
  mutate("Empiric Ratio" = Inverted/Tandem)%>%
  filter(!is.na(`Empiric Ratio`)) %>%
  rowwise() %>%
  mutate(
    total = Inverted + Tandem,
    p = list(c(`Expected Ratio`, 1) / (`Expected Ratio` + 1)),
    chisq = list(chisq.test(x = c(Inverted, Tandem), p = p)),
    p_value = chisq$p.value,
    stat = chisq$statistic
  ) %>%
  ungroup() %>%
  select(n, Family, Inverted, Tandem, `Empiric Ratio`, `Expected Ratio`, stat, p_value) %>%
  mutate(p.adjust = p.adjust(p_value, method = "fdr"))

suppFig3c = cluster_count_per_gene %>%
  pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
  # only take rows where both inverted and tandem exist
  filter(complete.cases(.)) %>%
  mutate("Empirical Ratio" = Inverted / Tandem) %>%
  # compute expected tandem elements
  inner_join(expected_inverted_to_tandem_ratio) %>%
  # only use cases where there are at least 10 tandem regions
  filter(Tandem > 10) %>%
  select(-Tandem, -Inverted) %>%
  mutate(RatioOfRatios = `Empirical Ratio` / `Expected Ratio`) %>%
  ggplot(aes(x = n, y = RatioOfRatios)) +
  geom_col(position = position_dodge(), color = "black", fill = index_colors[2]) +
  facet_grid(.~Family) +
  theme_custom() +
  scale_x_continuous(breaks = c(2:7)) +
  geom_text(aes(label= scales::label_percent()(signif(RatioOfRatios, digits = 3))), 
            position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5) +
  scale_y_continuous(labels = scales::label_percent()) +
  ggtitle("Depletion of inverted B1 and B2 clusters", subtitle = "") +
  xlab("Number of elements in cluster per gene") +
  ylab("% inverted clusters found empirically relative to expected") +
  ggpubr::stat_pvalue_manual(data = chiseq_pvals %>%
                               filter(Tandem > 10) %>%
                               mutate(group1 = "Empirical Ratio",
                                      group2 = "Expected Ratio",
                                      xmin = as.numeric(n),
                                      xmax = as.numeric(n),
                                      y.position = `Empiric Ratio`/`Expected Ratio` + 0.05,
                                      label = paste0("p = ", as.character(signif(p_value, digits = 3)))) %>%
                               filter(label!="ns"),
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             remove.bracket = TRUE,
                             inherit.aes = FALSE)

# Join --------------------------------------------------------------------
library(cowplot)
suppFig3 <- plot_grid(suppFig3b +
                      # inset
                      patchwork::inset_element(suppFig3b_inset_b1,
                                               left   = 0.15, right = 0.43,
                                               bottom = 0.52, top   = 0.995) +
                      # inset
                      patchwork::inset_element(suppFig3b_inset_b2,
                                               left   = 0.72, right = 0.87,
                                               bottom = 0.52, top   = 0.995) , 
                      suppFig3c, rel_widths = c(2, 1),
                      labels=c("A", "B"), ncol = 2, nrow = 1, align = 'hv', axis = "lbt", label_size=18)
save_plot(file.path(out_plots,"SuppFig2.3_MouseB1B2.pdf"), suppFig3, ncol = 1, nrow = 3, base_height = 3, base_width = 15, bg = "white")
save_plot(file.path(out_plots,"SuppFig2.3_MouseB1B2.png"), suppFig3, ncol = 1, nrow = 3, base_height = 3, base_width = 15, bg = "white")




