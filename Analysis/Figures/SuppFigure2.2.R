
rm(list = ls(all = TRUE))


# PRE-RUN ===========================================================================================
# *************************************************************************************************

# # Create longest 3'UTR ----------------------------------------------------
# fread("ucscHg38_RefSeqCurated_3pUTR.17022021.bed") %>% 
#   mutate(ID = str_remove(V4, "_utr3_.*"),
#          region_length = V3-V2)%>% 
#   filter(!grepl(x = V1, "_")) %>%
#   inner_join(all_utr3_to_genes %>%distinct() , by = join_by("ID" == "name")) %>%
#   group_by(name2) %>%
#   slice_max(order_by = region_length, n = 1) %>%
#   select(-ID, -V4) %>%
#   distinct%>%
#   select(V1, V2, V3, name2, region_length, V6) %>%
#   fwrite(file = "ucscHg38_RefSeqCurated_3pUTR_LongestIsoform.17022021.bed", 
#          quote = F, sep = "\t", col.names = F)

# # then run command
###findOppositeOrientationRepeatsInRegions.py -r ucscHg38_repeatMasker_SINEAlu.17022021.sorted.bed.gz -i ucscHg38_RefSeqCurated_3pUTR_LongestIsoform.17022021.bed -o Validations --merge  --prefix "hg38.Alu3pUTR_17022021.test." --min_rep_len 200


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

### create directories ------------------------------------------------------
#create paths, change
out_plots = file.path(outdir, out_plots)
out_stats = file.path(outdir, out_stats)
#create directories, if do not exist
dir.create(out_plots, recursive = T)
dir.create(out_stats, recursive = T)


# Load intermediate -------------------------------------------------------

# 3'UTR to gene - unmerged to account for internal isoforms
utr3_to_genes = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RegionsWithOppositeOrientationRepeat.sorted.csv", stringsAsFactors = F)

tandem_utr3_to_genes = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RegionsWithSameOrientationRepeat.sorted.csv", stringsAsFactors = F)


# alu inverted regions 
alus = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RepeatsInOppositeOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

alus_merged = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RepeatsInOppositeOrientationAtRegions.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# alu tandem regions 
tandem_alus = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RepeatsInSameOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

tandem_alus_merged = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RepeatsInSameOrientationAtRegions.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# 3'UTR regions
utr3s = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RegionsWithOppositeOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

tandem_utr3s = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.RegionsWithSameOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

# intersection results
intersection_res = fread("minLen200_analysis/hg38.Alu3pUTR_17022021_minLen200.LongestIsoform.intersection_results.csv", stringsAsFactors = F) %>%
  mutate("RefSeqID" = str_remove(region_name, "_utr3_.*")) %>%
  distinct()

# get only the relevant intersections
intersection_res_cytoindex_regions = intersection_res %>%
  filter(region %in% utr3s, repeat_cropped %in% alus) %>%
  # inner_join(utr3_to_genes, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
  rename("UTR3Region" = region,
         "UTR3RegionLength" = region_length,
         # "GeneStrand" = Strand,
         GeneSymbol = RefSeqID,
         "AluElement" = `repeat`,
         "AluElementCroppedToUTR3" = repeat_cropped,
         "AluElementOrientation" = repeat_strand,
         "AluElementLength" = repeat_length) %>%
  mutate("AluElementCroppedToUTR3Length" = as.numeric(AluElementCroppedToUTR3 %>% str_extract(pattern = "-[0-9]+") %>%
                                                        str_extract( pattern = "[0-9]+")) - 
           as.numeric(AluElementCroppedToUTR3 %>% 
                        str_extract(pattern = ":[0-9]+-") %>%
                        str_extract(pattern = "[0-9]+")),) %>%
  select(-contains("_"))#, -RegionLength)

# get only the relevant intersections
tandem_intersection_res_cytoindex_regions = intersection_res %>%
  filter(region %in% tandem_utr3s, repeat_cropped %in% tandem_alus) %>%
  # inner_join(tandem_utr3_to_genes, by = join_by("region" == "Region", "RefSeqID"== "RefSeqID"))%>%
  rename("UTR3Region" = region,
         "UTR3RegionLength" = region_length,
         # "GeneStrand" = Strand,
         GeneSymbol = RefSeqID,
         "AluElement" = `repeat`,
         "AluElementCroppedToUTR3" = repeat_cropped,
         "AluElementOrientation" = repeat_strand,
         "AluElementLength" = repeat_length) %>%
  # inner_join(alu_families) %>%
  # inner_join(alu_divergence) %>%
  select(-contains("_"))#, -RegionLength)

# no refseq id and gene \ gene length - just different 3'UTRs per GeneSymbol
# intersection_res_cytoindex_regions_unique = intersection_res_cytoindex_regions %>%
#   select(-RefSeqID, -Gene, -GeneLength) %>%
#   distinct()

# tandem_intersection_res_cytoindex_regions_unique = tandem_intersection_res_cytoindex_regions %>%
#   select(-RefSeqID, -Gene, -GeneLength, ) %>%
#   distinct()


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")
chosen_colors_trio = c("white", index_colors)

# Supp. Fig. 2.2A -----------------------------------------------------------------

### count clusters ----
cluster_count_per_gene = intersection_res_cytoindex_regions %>% 
  group_by(GeneSymbol, UTR3Region) %>% 
  tally %>% 
  group_by(GeneSymbol)%>%
  # take largest cluster per gene
  slice_max(order_by = n, with_ties = F) %>%
  mutate(Group = "Inverted") %>%
  bind_rows(tandem_intersection_res_cytoindex_regions %>% 
              group_by(GeneSymbol, UTR3Region) %>% 
              tally %>% 
              group_by(GeneSymbol)%>%
              # take largest cluster per gene
              slice_max(order_by = n, with_ties = F) %>%
              mutate(Group = "Tandem")) %>%
  group_by(n, Group) %>%
  tally(name = "Number of Genes") %>%
  ungroup()


fig5b = cluster_count_per_gene %>%
  filter(n > 1) %>%
  complete(Group, n = 2:21, fill = list(`Number of Genes` = 0)) %>%
  # filter(n>1) %>%
  ggplot(aes(x = n, y = `Number of Genes`, fill = Group)) +
  # geom_bar(position = position_dodge()) +
  geom_col(position = position_dodge(), color= "black") +
  theme_custom() +
  scale_x_continuous(breaks = c(2:21)) +
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
  guides(color = "none")   +
  theme(legend.title = element_blank())

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


fig5b_inset = cluster_count_per_gene %>%
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



# Supp. Fig. 2.2B -----------------------------------------------------------------

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

fig5c = cluster_count_per_gene %>%
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
                                      y.position = 0.7,
                                      label = case_when(p_value < 0.05  ~ paste0("p = ", as.character(signif(p_value, digits = 3))),
                                        TRUE            ~ "ns")) %>%
                               filter(label!="ns"),
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             remove.bracket = TRUE,
                             inherit.aes = FALSE)

# Supp. Fig. 2.2C -----------------------------------------------------------------

# fig5d = tribble(~Group, ~Type,  ~n,
#                 # "Inverted Alu Cluster", "Merged Clusters",  length(alus_merged),
#                 "Inverted", "Alu elements",  length(alus),
#                 # "Tandem Alu Cluster", "Merged Clusters",  length(tandem_alus_merged),
#                 "Tandem", "Alu elements",  length(tandem_alus)) %>%
#   ggplot(aes(x = Group, y = n, fill = Group))+
#   geom_col(position = position_dodge(), color = "black") +
#   theme_custom() +
#   scale_fill_manual(values = index_colors) +
#   geom_text(aes(label = n, vjust = -0.5), position = position_dodge(width = 0.7)) +
#   ggtitle("Number of elements in orientation group") +
#   ylab("Number of Alu elements")


### count elements in clusters ----
element_count_per_gene = intersection_res_cytoindex_regions %>% 
  group_by(GeneSymbol, UTR3Region) %>% 
  add_tally(name = "ClusterSize") %>% 
  group_by(AluElement)%>%
  # take largest cluster per gene
  slice_max(order_by = ClusterSize, with_ties = F)%>%
  mutate(Group = "Inverted") %>%
  bind_rows(tandem_intersection_res_cytoindex_regions %>% 
              group_by(GeneSymbol, UTR3Region) %>% 
              add_tally(name = "ClusterSize") %>% 
              group_by(AluElement)%>%
              # take largest cluster per gene
              slice_max(order_by = ClusterSize, with_ties = F) %>%
              mutate(Group = "Tandem"))  %>%
  mutate(`Cluster Size` = forcats::fct_other(as.factor(ClusterSize), keep = c("1", "2"), other_level = "3+") %>% 
           factor(levels = c("3+", "2", "1")))

fig5d = element_count_per_gene %>%
  ggplot(aes(x = Group, fill = `Cluster Size`))+
  geom_bar(position = position_stack(reverse = TRUE), color = "black") +
  theme_custom() +
  scale_fill_manual(values = rev(chosen_colors_trio)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  geom_text(aes(label = after_stat(count)), position = position_stack(reverse = TRUE), stat='count', vjust = 2.5) +
  ggtitle("Elements in 3'UTRs") +
  ylab(expression("Number of  " * italic("Alu") * " elements"))  +
  theme(axis.title.x = element_blank())

# Supp. Fig. 2.2D -----------------------------------------------------------------

inverted_cluster_genes = intersection_res_cytoindex_regions$GeneSymbol %>% unique
tandem_cluster_genes = tandem_intersection_res_cytoindex_regions$GeneSymbol %>% unique

fig5e = tribble(~Group, ~`Cluster Types`,  ~n,
                "Inverted", "Inverted only",  length(setdiff(inverted_cluster_genes, tandem_cluster_genes)),
                # "Inverted", "Both",  length(intersect(inverted_cluster_genes, tandem_cluster_genes)),
                "Tandem", "Tandem only",  length(setdiff(tandem_cluster_genes, inverted_cluster_genes)),
                # "Tandem", "Both",  length(intersect(tandem_cluster_genes, inverted_cluster_genes))
) %>%
  ggplot(aes(x = Group, y = n, fill = `Cluster Types`))+
  geom_col(position = position_stack(), color = "black") +
  theme_custom(legend_position = "none") +
  scale_fill_manual(values = c(index_colors)) +
  # scale_fill_manual(values = c("white", index_colors)) +
  geom_text(aes(label = n), vjust = 1, position = position_stack(vjust = 0.85)) +
  ggtitle(expression(bold("Genes with 3'UTR " * bolditalic("Alu") * " cluster"))) +
  ylab("Number of genes") +
  theme(plot.title = element_text(hjust = 0.9)) +
  guides(fill = guide_legend(title = "Cluster\nTypes", ncol = 1))   +
  theme(axis.title.x = element_blank())


# Supp. Fig. 2.2E ---------------------------------------------------------
alu_elements = intersection_res_cytoindex_regions$AluElement%>%unique
tandem_alu_elements = tandem_intersection_res_cytoindex_regions$AluElement%>%unique
fig5f = tribble(~Group, ~`Cluster Types`,  ~n,
                "Inverted", "Inverted only",  length(setdiff(alu_elements, tandem_alu_elements)),
                "Inverted", "Both",  length(intersect(alu_elements, tandem_alu_elements)),
                "Tandem", "Tandem only",  length(setdiff(tandem_alu_elements, alu_elements)),
                "Tandem", "Both",  length(intersect(tandem_alu_elements, alu_elements))) %>%
  ggplot(aes(x = Group, y = n, fill = `Cluster Types`))+
  geom_col(position = position_stack(), color = "black") +
  theme_custom() +
  scale_fill_manual(values = c("white", index_colors)) +
  geom_text(aes(label = n), vjust = -0.5, position = position_stack(vjust = 0.5)) +
  ggtitle("Elements per orientation") +
  ylab(expression("Number of " * italic("Alu") * " elements")) +
  guides(fill = guide_legend(title = "Cluster Types", nrow = 1))   +
  theme(axis.title.x = element_blank(),
        legend.justification = c(1.3, 0))

# Join --------------------------------------------------------------------
library(cowplot)
fig5 <- plot_grid(plot_grid(fig5b +
                              # inset
                              patchwork::inset_element(fig5b_inset,
                                                       left   = 0.715, right = 0.99,
                                                       bottom = 0.52, top   = 0.995), 
                            labels=c("A"), ncol = 1, nrow = 1, align = 'hv', label_size=18),
                  plot_grid(fig5c, fig5d, fig5e, fig5f, rel_widths = c(3, 2, 2, 2), 
                            labels=c("B", "C", "D", "E"), ncol = 4, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  rel_heights = c(1, 0.8), 
                  labels=c("", "", ""), ncol = 1, nrow = 2, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"SuppFig2.2_LongestIsoformPreCalc.pdf"), fig5, ncol = 1, nrow = 3, base_height = 6, base_width = 15)
save_plot(file.path(out_plots,"SuppFig2.2_LongestIsoformPreCalc.png"), fig5, ncol = 1, nrow = 3, base_height = 6, base_width = 15)




