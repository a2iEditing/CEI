
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


### Human -------------------------------------------------------------------

# 3'UTR to gene - unmerged to account for internal isoforms
utr3_to_genes = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithOppositeOrientationRepeat.sorted.csv", stringsAsFactors = F)

tandem_utr3_to_genes = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithSameOrientationRepeat.sorted.csv", stringsAsFactors = F)


# alu inverted regions 
alus = fread("hg38.Alu3pUTR_minLen200_17022021.RepeatsInOppositeOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# alu tandem regions 
tandem_alus = fread("hg38.Alu3pUTR_minLen200_17022021.RepeatsInSameOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# 3'UTR regions
utr3s = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithOppositeOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

tandem_utr3s = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithSameOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

# intersection results
intersection_res = fread("analysis/hg38.Alu3pUTR_minLen200_17022021.intersection_results.csv", stringsAsFactors = F) %>%
  mutate("RefSeqID" = str_remove(region_name, "_utr3_.*")) %>%
  distinct()

# get only the relevant intersections
intersection_res_cytoindex_regions = intersection_res %>%
  filter(region %in% utr3s, repeat_cropped %in% alus) %>%
  inner_join(utr3_to_genes, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
  rename("UTR3Region" = region,
         "UTR3RegionLength" = region_length,
         "AluElement" = `repeat`,
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


### Mouse -------------------------------------------------------------------
utr3_to_genes_mouse = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithOppositeOrientationRepeat.sorted.csv", stringsAsFactors = F)

tandem_utr3_to_genes_mouse = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithSameOrientationRepeat.sorted.csv", stringsAsFactors = F)


# alu inverted regions 
alus_mouse = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RepeatsInOppositeOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# alu tandem regions 
tandem_alus_mouse = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RepeatsInSameOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# 3'UTR regions
utr3s_mouse = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithOppositeOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

tandem_utr3s_mouse = fread("mm10.AluB1AndB2_3pUTR_minLen120_17022021.RegionsWithSameOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

# intersection results
intersection_res_mouse = fread("analysis/mm10.AluB1AndB2_3pUTR_minLen120_17022021.intersection_results.csv", stringsAsFactors = F) %>%
  mutate("RefSeqID" = str_remove(region_name, "_utr3_.*")) %>%
  distinct() %>%
  mutate(repeat_score = case_match(repeat_score, "Alu" ~ "B1", 
                                   .default = repeat_score))

# get only the relevant intersections
intersection_res_cytoindex_regions_mouse = intersection_res_mouse %>%
  filter(region %in% utr3s_mouse, repeat_cropped %in% alus_mouse) %>%
  inner_join(utr3_to_genes_mouse, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
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
tandem_intersection_res_cytoindex_regions_mouse = intersection_res_mouse %>%
  filter(region %in% tandem_utr3s_mouse, repeat_cropped %in% tandem_alus_mouse) %>%
  inner_join(tandem_utr3_to_genes_mouse, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
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
intersection_res_cytoindex_regions_unique_mouse = intersection_res_cytoindex_regions_mouse %>%
  select(-RefSeqID, -Gene, -GeneLength) %>%
  distinct()

tandem_intersection_res_cytoindex_regions_unique_mouse = tandem_intersection_res_cytoindex_regions_mouse %>%
  select(-RefSeqID, -Gene, -GeneLength, ) %>%
  distinct()


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")
chosen_colors_trio = c("white", index_colors)

# Fig. 2B -----------------------------------------------------------------

### count clusters ----
cluster_count_per_gene = intersection_res_cytoindex_regions_unique %>% 
  group_by(GeneSymbol, UTR3Region) %>% 
  tally %>% 
  group_by(GeneSymbol)%>%
  # take largest cluster per gene
  slice_max(order_by = n, with_ties = F) %>%
  mutate(Group = "Inverted") %>%
  bind_rows(tandem_intersection_res_cytoindex_regions_unique %>% 
              group_by(GeneSymbol, UTR3Region) %>% 
              tally %>% 
              group_by(GeneSymbol)%>%
              # take largest cluster per gene
              slice_max(order_by = n, with_ties = F) %>%
              mutate(Group = "Tandem")) %>%
  group_by(n, Group) %>%
  tally(name = "Number of Genes") %>%
  ungroup()


fig2b = cluster_count_per_gene %>%
  filter(n > 1) %>%
  complete(Group, n = 2:21, fill = list(`Number of Genes` = 0)) %>%
  # filter(n>1) %>%
  ggplot(aes(x = n, y = `Number of Genes`, fill = Group)) +
  # geom_bar(position = position_dodge()) +
  geom_col(position = position_dodge(), color= "black") +
  theme_custom() +
  scale_x_continuous(breaks = c(2:21), expand = c(0.02, 0.01)) +
  scale_y_continuous(breaks = seq(0, 550, 100)) +
  geom_text(aes(label= if_else(`Number of Genes` == 0, "", as.character(`Number of Genes`)),
                color = Group), position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5) +
  # scale_fill_gradient(low = case_control_colors[1], high = case_control_colors[2]) +
  scale_color_manual(values = index_colors) +
  scale_fill_manual(values = index_colors) +
  ggtitle(expression(bold("Distribution of genes with " * bolditalic("Alu") * " clusters")),
          subtitle = "Largest cluster per gene (GeneSymbol)") +
  xlab("Number of elements in cluster per gene") +
  ylab("Number of genes") +
  guides(color = "none") +
  theme(legend.title = element_blank()) +
  geom_label(data = cluster_count_per_gene %>%
              filter(n >=2, n<=4) %>%
              pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
              mutate(EmpiricalRatio = Inverted/Tandem), 
            aes(x = n,
                # pick the taller of the two bars so the label clears both:
                y = 0,
                label = sprintf("%.2f", EmpiricalRatio)),
            position     = position_dodge(width = 0.9),
            vjust        = -0.5,
            hjust = 0.5,
            inherit.aes  = FALSE) 


### expected inverted vs. tandem ratio -----
expected_inverted_to_tandem_ratio = data.frame(n = 2:7, 
                                               "Expected Ratio" = 2^(2:7 - 1) - 1,
                                               check.names = F) 

expected_inverted_relative_to_tandem = cluster_count_per_gene %>%
  filter(n > 1, Group == "Tandem") %>%
  inner_join(expected_inverted_to_tandem_ratio) %>%
  mutate(Group = "Expected",
         `Number of Genes` = `Number of Genes` * `Expected Ratio`) %>%
  select(-`Expected Ratio`)
  

fig2b_inset = cluster_count_per_gene %>%
  filter(n > 1, n <= 7, Group == "Inverted") %>%
  complete(Group, n = 2:7, fill = list(`Number of Genes` = 0)) %>%
  mutate(Group = "Observed") %>%
  bind_rows(expected_inverted_relative_to_tandem) %>%
  ggplot(aes(x = n, y = `Number of Genes`, fill = Group)) +
  # geom_bar(position = position_dodge()) +
  geom_col(position = position_dodge(), color= "black") +
  theme_custom() +
  scale_x_continuous(breaks = c(2:7), expand = c(0.02, 0.01)) +
  scale_y_continuous(breaks = seq(0, 550, 100)) +
  expand_limits(y = 590) +
  geom_text(aes(label= if_else(`Number of Genes` == 0, "", as.character(`Number of Genes`)),
                color = Group), position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5, size = 3) +
  # scale_fill_gradient(low = case_control_colors[1], high = case_control_colors[2]) +
  scale_color_manual(values = c("black", index_colors[1])) +
  scale_fill_manual(values = c("white", index_colors[1])) +
  ggtitle("Expected vs. observed inverted genes",
          subtitle = "Relative to observed tandem genes") +
  xlab("# of elements in cluster per gene") +
  ylab("# of genes") +
  guides(color = "none") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA))


# Fig. 2C -----------------------------------------------------------------


# # Chi squared
chiseq_pvals = cluster_count_per_gene %>%
  pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
  # compute expected tandem elements
  inner_join(expected_inverted_to_tandem_ratio) %>%
  mutate("Empiric Ratio" = Inverted/Tandem)%>%
  rowwise() %>%
  mutate(
    total = Inverted + Tandem,
    p = list(c(`Expected Ratio`, 1) / (`Expected Ratio` + 1)),
    chisq = list(chisq.test(x = c(Inverted, Tandem), p = p)),
    p_value = chisq$p.value,
    stat = chisq$statistic
  ) %>%
  ungroup() %>%
  select(n, Inverted, Tandem, `Empiric Ratio`, `Expected Ratio`, stat, p_value) %>%
  mutate(p.adjust = p.adjust(p_value, method = "fdr"))

fig2c = cluster_count_per_gene %>%
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
  theme_custom() +
  scale_x_continuous(breaks = c(2:7)) +
  geom_text(aes(label= scales::label_percent()(signif(RatioOfRatios, digits = 3))), 
            position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5) +
  scale_y_continuous(labels = scales::label_percent()) +
  ggtitle(expression(bold("Depletion of inverted " * bolditalic("Alu") * " clusters"))) +
  xlab("Number of elements in cluster per gene") +
  ylab("% inverted clusters found empirically relative to expected") +
  ggpubr::stat_pvalue_manual(data = chiseq_pvals %>%
                               mutate(group1 = "Empirical Ratio",
                                      group2 = "Expected Ratio",
                                      xmin = as.numeric(n),
                                      xmax = as.numeric(n),
                                      y.position = 0.67,
                                      label = case_when(p_value < 0.05  ~ paste0("p = ", as.character(signif(p_value, digits = 3))),
                                        TRUE            ~ "ns")) %>%
                               filter(label!="ns"),
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             remove.bracket = TRUE,
                             inherit.aes = FALSE)

# Fig. 2D -----------------------------------------------------------------

### count elements in clusters ----
element_count_per_gene = intersection_res_cytoindex_regions_unique %>% 
  # add count of elements per 3'UTR
  group_by(GeneSymbol, UTR3Region) %>% 
  add_tally(name = "ClusterSize") %>% 
  group_by(AluElement)%>%
  # for each Alu element, put in the group of the largest cluster
  slice_max(order_by = ClusterSize, with_ties = F)%>%
  mutate(Group = "Inverted") %>%
  bind_rows(tandem_intersection_res_cytoindex_regions_unique %>% 
              # add count of elements per 3'UTR
              group_by(GeneSymbol, UTR3Region) %>% 
              add_tally(name = "ClusterSize") %>% 
              group_by(AluElement)%>%
              # for each Alu element, put in the group of the largest cluster
              slice_max(order_by = ClusterSize, with_ties = F) %>%
              mutate(Group = "Tandem"))  %>%
  mutate(`Cluster Size` = forcats::fct_other(as.factor(ClusterSize), keep = c("1", "2"), other_level = "3+") %>% 
           factor(levels = c("3+", "2", "1")))

fig2d = element_count_per_gene %>%
  ggplot(aes(x = Group, fill = `Cluster Size`))+
  geom_bar(position = position_stack(reverse = TRUE), color = "black") +
  theme_custom() +
  scale_fill_manual(values = rev(chosen_colors_trio)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_text(aes(label = after_stat(count)), position = position_stack(reverse = TRUE), stat='count', vjust = 2.5) +
  ggtitle("Elements in 3'UTRs",
          subtitle = "Human") +
  guides(fill = guide_legend(title = "Cluster\nSize")) +
  ylab(expression("Number of  " * italic("Alu") * " elements")) +
  theme(axis.title.x = element_blank()) +
  theme(plot.title = element_text(hjust = 1.3))

# Fig. 2E -----------------------------------------------------------------

### count elements in clusters ----
element_count_per_gene_mouse = intersection_res_cytoindex_regions_unique_mouse %>% 
  group_by(GeneSymbol, UTR3Region, Family) %>% 
  add_tally(name = "ClusterSize") %>% 
  group_by(AluElement, Family)%>%
  # take largest cluster per gene
  slice_max(order_by = ClusterSize, with_ties = F)%>%
  mutate(Group = "Inverted") %>%
  bind_rows(tandem_intersection_res_cytoindex_regions_unique_mouse %>% 
              group_by(GeneSymbol, UTR3Region, Family) %>% 
              add_tally(name = "ClusterSize") %>% 
              group_by(AluElement, Family)%>%
              # take largest cluster per gene
              slice_max(order_by = ClusterSize, with_ties = F) %>%
              mutate(Group = "Tandem"))  %>%
  mutate(`Cluster Size` = forcats::fct_other(as.factor(ClusterSize), keep = c("1", "2"), other_level = "3+") %>% 
           factor(levels = c("3+", "2", "1")))

fig2e = element_count_per_gene_mouse %>%
  ggplot(aes(x = Group, fill = `Cluster Size`))+
  facet_grid(.~ Family) +
  geom_bar(position = position_stack(reverse = TRUE), color = "black") +
  theme_custom() +
  scale_fill_manual(values = rev(chosen_colors_trio)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_text(aes(label = after_stat(count)), position = position_stack(reverse = TRUE), stat='count', vjust = 1.4) +
  ggtitle("Elements in 3'UTRs",
          subtitle = "Mouse") +
  ylab("Number of B1 and B2 elements") +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(title = "Cluster\nSize")) 


# Fig. 2F & 2G -----------------------------------------------------------------

fig2f = ggplot(data.frame(x = 1:1,
                              colors = c("white"),
                              text = "TODO: add\nDeepVenn Human"), 
                   aes(x, y = 0, fill = colors, label = text)) +
  geom_tile(width = .9, height = .9) + # make square tiles
  geom_text(color = "Red", size = 7) + # add white text in the middle
  scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
  # coord_fixed() + # make sure tiles are square
  theme_custom_void() +
  ggtitle(expression(bold("Genes with 3'UTR " * bolditalic("Alu") * " cluster"))) 

fig2g = ggplot(data.frame(x = 1:1,
                          colors = c("white"),
                          text = "TODO: add\nDeepVenn Mouse"), 
               aes(x, y = 0, fill = colors, label = text)) +
  geom_tile(width = .9, height = .9) + # make square tiles
  geom_text(color = "Red", size = 7) + # add white text in the middle
  scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
  # coord_fixed() + # make sure tiles are square
  theme_custom_void() +
  ggtitle("Genes with 3'UTR B1 and B2 clusters")



# Join --------------------------------------------------------------------
library(cowplot)
fig2 <- plot_grid(plot_grid(fig2b +
                              # inset
                              patchwork::inset_element(fig2b_inset,
                                                       left   = 0.075, right = 0.35,
                                                       bottom = 0.52, top   = 0.995), 
                            fig2c,
                            rel_widths = c(9, 3),
                            labels=c("B", "C"), ncol = 2, nrow = 1, align = 'hv', axis = "btl", label_size=18),
                  plot_grid(fig2d, fig2e, 
                            fig2f, fig2g, 
                            rel_widths = c(1.2, 2, 2.5, 4), 
                            labels=c("D", "E", "F", "G"), ncol = 4, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  rel_heights = c(1, 0.8), 
                  labels=c("", ""), ncol = 1, nrow = 2, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"Figure2B-E.pdf"), fig2, ncol = 1, nrow = 3, base_height = 6.7, base_width = 20)
save_plot(file.path(out_plots,"Figure2B-E.png"), fig2, ncol = 1, nrow = 3, base_height = 6.7, base_width = 20)




