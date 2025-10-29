
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

### Element Analysis -----
# alu families
alu_families = fread("ucscHg38_repeatMasker_subset.Alu.17022021.sorted.bed", sep = "\t", header = F) %>%
  mutate(AluElement = paste0(V1, ":", V2, "-", V3),
         AluElementFamily = V4,
         AluElementFamilyGroup = case_when(str_detect(AluElementFamily, "^AluJ") ~ "AluJ",
                                           str_detect(AluElementFamily, "^AluS") ~ "AluS",
                                           str_detect(AluElementFamily, "^AluY") ~ "AluY",
                                           TRUE ~ "Other")) %>%
  select(-starts_with("V"))

# divergence
alu_divergence = fread("ucscHg38_repeatMasker_subset.Alu.17022021.sorted.bed", sep = "\t", header = F) %>%
  mutate(AluElement = paste0(V1, ":", V2, "-", V3),
         AluElementDivergencePcnt = V8 / 1000,
         AluElementDeletionsPcnt = V9 / 1000,
         AluElementInsertionsPcnt = V8 / 1000) %>%
  select(-starts_with("V")) 

# 3'UTR to gene - unmerged to account for internal isoforms
utr3_to_genes = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithOppositeOrientationRepeat.sorted.csv", stringsAsFactors = F)

tandem_utr3_to_genes = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithSameOrientationRepeat.sorted.csv", stringsAsFactors = F)


# alu inverted regions 
alus = fread("hg38.Alu3pUTR_minLen200_17022021.RepeatsInOppositeOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

alus_merged = fread("hg38.Alu3pUTR_minLen200_17022021.RepeatsInOppositeOrientationAtRegions.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# alu tandem regions 
tandem_alus = fread("hg38.Alu3pUTR_minLen200_17022021.RepeatsInSameOrientationAtRegions.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

tandem_alus_merged = fread("hg38.Alu3pUTR_minLen200_17022021.RepeatsInSameOrientationAtRegions.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
  mutate(AluRepeat = paste0(V1, ":", V2, "-", V3)) %>%
  pull(AluRepeat)

# 3'UTR regions
utr3s = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithOppositeOrientationRepeat.sorted.bed", stringsAsFactors = F, header = F) %>%
  mutate(UTR3Region = paste0(V1, ":", V2, "-", V3)) %>%
  pull(UTR3Region)

tandem_utr3s = fread("hg38.Alu3pUTR_minLen200_17022021.RegionsWithSameOrientationRepeat.sorted.merged.bed", stringsAsFactors = F, header = F) %>%
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
  inner_join(alu_families) %>%
  inner_join(alu_divergence) %>%
  # left_join(alu_conservation) %>%
  mutate(#AluElementConservationScore = replace_na(AluElementConservationScore, 0),
         "AluElementCroppedToUTR3Length" = as.numeric(AluElementCroppedToUTR3 %>% str_extract(pattern = "-[0-9]+") %>%
                                                        str_extract( pattern = "[0-9]+")) - 
           as.numeric(AluElementCroppedToUTR3 %>% 
                        str_extract(pattern = ":[0-9]+-") %>%
                        str_extract(pattern = "[0-9]+"))) %>%
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
  inner_join(alu_families) %>%
  inner_join(alu_divergence) %>%
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


# Supp. Fig. 1A -----------------------------------------------------------------


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

### expected inverted vs. tandem ratio -----
expected_tandem_to_inverted_ratio = data.frame(n = 2:7, 
                                               "Expected Ratio" = 2/(2^(2:7) - 2),
                                               check.names = F) 

# Chi squared
chiseq_pvals_suppFig1b = cluster_count_per_gene %>%
  pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
  # compute expected tandem elements
  inner_join(expected_tandem_to_inverted_ratio) %>%
  mutate("Empirical Ratio" = Tandem/Inverted)%>%
  rowwise() %>%
  mutate(
    total = Inverted + Tandem,
    p = list(c(`Expected Ratio`, 1) / (`Expected Ratio` + 1)),
    chisq = list(chisq.test(x = c(Tandem, Inverted), p = p)),
    p_value = chisq$p.value,
    stat = chisq$statistic
  ) %>%
  ungroup() %>%
  select(n, Inverted, Tandem, `Empirical Ratio`, `Expected Ratio`, stat, p_value) %>%
  mutate(p.adjust = p.adjust(p_value, method = "fdr"))

suppFig1a = cluster_count_per_gene %>%
  pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
  # only take rows where both inverted and tandem exist
  filter(complete.cases(.)) %>%
  mutate("Empirical Ratio" = Tandem / Inverted) %>%
  # compute expected tandem elements
  inner_join(expected_tandem_to_inverted_ratio) %>%
  select(-Tandem, -Inverted) %>%
  pivot_longer(cols = contains("Ratio"), names_to = "Group", values_to = "Ratio") %>%
  mutate(Ratio_label = round(Ratio, 3) %>%
           signif(digits = 2) %>%
           format %>%
           str_remove("0$")) %>%
  ggplot(aes(x = n, y = Ratio, fill = Group)) +
  geom_col(position = position_dodge(), color = "black") +
  theme_custom() +
  scale_x_continuous(breaks = c(2:7)) +
  geom_text(aes(label= Ratio_label), position = position_dodge(width = 0.9), vjust=-1, hjust = 0.5) +
  # ggpubr::stat_compare_means(label = "p.signif", symnum.args = list(cutpoints = c(0, Inf), symbols = c("****"))) +
  ggpubr::stat_pvalue_manual(data = chiseq_pvals_suppFig1b %>%
                                mutate(group1 = "Empirical Ratio",
                                       group2 = "Expected Ratio",
                                       xmin = as.numeric(n) - 0.25,
                                       xmax = as.numeric(n) + 0.25,
                                       y.position = 2.1,
                                       label = paste0("p = ", as.character(signif(p_value, digits = 2)))),
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             # bracket.size = 0,
                             tip.length = 0.02,
                             inherit.aes = FALSE) +
  ggpubr::stat_pvalue_manual(data = cluster_count_per_gene %>%
                               pivot_wider(names_from = Group, values_from = `Number of Genes`) %>%
                               # only take rows where both inverted and tandem exist
                               filter(complete.cases(.)) %>%
                               mutate("Empirical Ratio" = Inverted / Tandem) %>%
                               # compute expected tandem elements
                               inner_join(expected_tandem_to_inverted_ratio) %>%
                               mutate(`Expected Ratio` = 1/`Expected Ratio`) %>%
                               # only use cases where there are at least 10 tandem regions
                               filter(Tandem > 10) %>%
                               select(-Tandem, -Inverted) %>%
                               mutate(label = if_else(n<=4,
                                                      paste0(scales::label_percent()(`Expected Ratio` / `Empirical Ratio`),
                                                             "\nmore than\nexpected"),
                                                      ""),
                                      group1 = "Empirical Ratio",
                                      group2 = "Expected Ratio",
                                      xmin = as.numeric(n) - 0.2,
                                      xmax = as.numeric(n) + 0.2,
                                      y.position = 1.75),
                             label = "label",
                             xmin = "xmin",
                             xmax = "xmax",
                             y.position = "y.position",
                             bracket.size = 0,
                             tip.length = 0,
                             inherit.aes = FALSE) +
  scale_fill_manual(values = case_control_colors) +
  ggtitle("Comparison of expected and observed ratios of tandem vs. inverted clusters") +
  xlab("Number of elements in cluster per gene") +
  ylab("Number of tandem clusters / number of inverted clusters") +
  theme(legend.title = element_blank())

# Supp. Fig. 1B -----------------------------------------------------------------

divergence_score = bind_rows(intersection_res_cytoindex_regions_unique %>%
                               distinct(AluElement, AluElementDivergencePcnt, AluElementFamilyGroup) %>%
                               mutate(Orientation = "Inverted"),
                             tandem_intersection_res_cytoindex_regions_unique %>%
                               distinct(AluElement, AluElementDivergencePcnt, AluElementFamilyGroup) %>%
                               mutate(Orientation = "Tandem"))

# compute p-values
divergence_score_pvals_suppFig2.1b = divergence_score %>% 
  filter(AluElementFamilyGroup!="Other")%>%
  group_by(AluElementFamilyGroup) %>%
  summarise(p_value = wilcox.test(AluElementDivergencePcnt ~ Orientation)$p.value)


suppFig1b = divergence_score %>%
  # filter(Orientation!="Other") %>%
  ggplot(aes(x = Orientation, y = AluElementDivergencePcnt, fill = Orientation)) +
  geom_boxplot(outliers = F) + 
  facet_grid(.~AluElementFamilyGroup, space = "free", scales = "free") +
  # geom_violin(draw_quantiles = 0.5, trim = T) +
  theme_custom(legend_position = "none") +
  scale_y_continuous(labels = scales::label_percent(), breaks = seq(0, 0.25, 0.05)) +
  scale_fill_manual(values = index_colors) + #, labels = c("Inverted elements", "Tandem Alu elements")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab(expression(italic("Alu") * " element divergence %")) +
  ggtitle(expression(bold("Divergence percentage per major " * bolditalic("Alu") * " family"))) +
  # ggpubr::stat_compare_means(method = "wilcox.test", 
  #                            comparisons = list(c(1,2)))
  ggpubr::stat_pvalue_manual(data = divergence_score_pvals_suppFig2.1b %>%
                               mutate(group1 = "Inverted",
                                      group2 = "Tandem",
                                      y.position = 0.27,
                                      label = paste0("p = ", as.character(signif(p_value, digits = 2)))),
                             label = "label",
                             y.position = "y.position",
                             # bracket.size = 0,
                             tip.length = 0.02,
                             inherit.aes = FALSE)



# suppFig1b = intersection_res_cytoindex_regions_unique %>% 
#   group_by(GeneSymbol, UTR3Region) %>% 
#   add_tally(name = "ClusterSize") %>%
#   # only take inverted clusters of size 3
#   filter(ClusterSize == 3) %>%
#   group_by(GeneSymbol, UTR3Region, AluElementOrientation) %>%
#   add_tally(name = "NumElementsInThisDirection") %>%
#   group_by(NumElementsInThisDirection, AluElementFamilyGroup)%>% 
#   tally %>%
#   mutate(NumElementsInThisDirection = case_when(NumElementsInThisDirection==1 ~ "Minority",
#                                                 NumElementsInThisDirection==2 ~ "Majority",
#                                                 .default = NA)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = NumElementsInThisDirection, values_from = n) %>%
#   mutate("MajorityToMinorityRatio" = Majority/Minority) %>%
#   # plot
#   ggplot(aes(x = AluElementFamilyGroup, y = MajorityToMinorityRatio, fill = AluElementFamilyGroup)) +
#   geom_col(color = "black") +
#   theme_custom() +
#   geom_text(aes(label = signif(MajorityToMinorityRatio, digits = 3)), vjust = -1) +
#   scale_fill_manual(values = c("grey30", index_colors)) +
#   xlab("Major Family") +
#   ylab("Majority (2 elements) to minority (1 element) ratio") +
#   ggtitle("Cluster orientation ratios",
#           subtitle = "For inverted clusters of size 3") +
#   guides(fill = guide_legend(title = "Family"))


# Supp. Fig. 1C -----------------------------------------------------------------

# genome-wide families info
alu_families_count = alu_families %>%
  filter(AluElementFamilyGroup != "Other") %>%
  group_by(AluElementFamilyGroup) %>%
  tally() 

alu_family_prop = alu_families_count %>%
  mutate(GenomeWideProp = n/sum(n)) %>%
  select(-n) %>%
  { setNames(pull(., GenomeWideProp), pull(., AluElementFamilyGroup)) }

# get data
curr_data_suppFig1c = intersection_res_cytoindex_regions_unique %>% 
  group_by(GeneSymbol, UTR3Region) %>% 
  add_tally(name = "ClusterSize") %>%
  # only take inverted clusters of size 3
  filter(ClusterSize <= 6) %>%
  ungroup() %>%
  distinct(AluElement, ClusterSize, AluElementFamilyGroup) %>%
  group_by(ClusterSize, AluElementFamilyGroup) %>%
  tally %>%
  mutate(ClusterSize = as.factor(ClusterSize))

# compute chi squared
chiseq_pvals_suppFig1c = curr_data_suppFig1c %>%
  pivot_wider(names_from = AluElementFamilyGroup, values_from = n) %>%
  mutate(chisq = list(chisq.test(x = c_across(AluJ:AluY), p = alu_family_prop)),
         statistic = map_dbl(chisq, "statistic"),
         p_value   = map_dbl(chisq, "p.value")) %>%
  ungroup() %>%
  mutate(p.adjust = p.adjust(p_value, method = "fdr"))

# get plot
suppFig1c = curr_data_suppFig1c %>%
  # add distribution in general genomic population
  bind_rows(alu_families %>%
              filter(AluElementFamilyGroup != "Other") %>%
              group_by(AluElementFamilyGroup) %>%
              tally %>%
              mutate(ClusterSize = as.factor("Genome-wide"))) %>%
  mutate(ClusterSize = forcats::fct_relevel(ClusterSize, "Genome-wide"),
         total = sum(n),
         pcnt = n/total,
         label = if_else(ClusterSize == "Genome-wide", 
                         paste0(scales::label_percent(accuracy = 0.1)(pcnt), "\n(", n, ")"), 
                         paste0(scales::label_percent(accuracy = 1)(pcnt), "\n(", n, ")")),
         text_angle = if_else(ClusterSize == "Genome-wide", 90, 0),
         text_hjust = if_else(text_angle == 90, 1.2, 0.5),
         text_vjust = if_else(text_angle == 90, 0.5, 1.5)) %>%
  ggplot(aes(x = forcats::fct_reorder(AluElementFamilyGroup, -pcnt), y = pcnt, fill = AluElementFamilyGroup)) +
    facet_grid(.~ClusterSize) +
    geom_col(color = "black") +
    theme_custom(legend_position = "none") +
    geom_text(aes(label= label, angle = text_angle, hjust  = text_hjust,
                  vjust  = text_vjust), position = position_dodge(width = 0.7))  +
    # geom_text(aes(label = signif(MajorityToMinorityRatio, digits = 3)), vjust = -1) +
    scale_fill_manual(values = c("white", index_colors)) +
    scale_y_continuous(labels = scales::label_percent()) +
    xlab("Major Family") +
    ylab("% of elements ") +
    ggtitle(expression(bold("Distribution of major " * bolditalic("Alu") * " families according to cluster size"))) +
    guides(fill = guide_legend(title = "Family")) +
  geom_text(data = chiseq_pvals_suppFig1c %>%
              select(ClusterSize, p_value) %>%
              mutate(y.position = 0.72,
                     label = paste0("p = ", signif(p_value, 2))),
            aes(x = "AluJ",            # place at the right end of each panel
                y = y.position, 
                label = label),
            inherit.aes = FALSE)

# Supp. Fig. 1D -----------------------------------------------------------------
curr_data_suppFig1d = intersection_res_cytoindex_regions_unique %>% 
  group_by(GeneSymbol, UTR3Region) %>% 
  add_tally(name = "ClusterSize") %>%
  # only take inverted clusters of size 3
  filter(ClusterSize == 3)

# select longest 3'UTR
chosen_ut3regions = curr_data_suppFig1d %>%
  group_by(GeneSymbol) %>%
  slice_max(order_by = UTR3RegionLength, n = 1, with_ties = F) %>%
  pull(UTR3Region)

curr_data_suppFig1d_count = curr_data_suppFig1d %>%
  filter(UTR3Region %in% chosen_ut3regions) %>%
  group_by(GeneSymbol, UTR3Region, AluElementOrientation) %>%
  add_tally(name = "NumElementsInThisDirection") %>%
  group_by(NumElementsInThisDirection, AluElementFamilyGroup) %>% 
  tally %>%
  mutate(Group = case_when(NumElementsInThisDirection==1 ~ "Minority",
                           NumElementsInThisDirection==2 ~ "Majority",
                           .default = NA)) %>%
  group_by(Group) %>%
  mutate(total = sum(n),
         pcnt = n/total,
         label = paste0(scales::label_percent(accuracy = 1)(pcnt), "\n(", n, ")"))

# compute chi squared
chiseq_pvals_suppFig1d = curr_data_suppFig1d_count %>%
  pivot_wider(names_from = AluElementFamilyGroup, values_from = n, id_cols = c(Group)) %>%
  tibble::column_to_rownames("Group") %>%
  as.matrix() %>%
  chisq.test

suppFig1d = curr_data_suppFig1d_count %>%
  # plot
  ggplot(aes(x = forcats::fct_reorder(AluElementFamilyGroup, -pcnt), y = pcnt, fill = AluElementFamilyGroup)) +
  facet_grid(.~Group) +
  geom_col(color = "black") +
  theme_custom(legend_position = "none") +
  geom_text(aes(label= label), vjust=1.5, position = position_dodge(width = 0.7))  +
  # geom_text(aes(label = signif(MajorityToMinorityRatio, digits = 3)), vjust = -1) +
  scale_fill_manual(values = c("white", index_colors)) +
  scale_y_continuous(labels = scales::label_percent()) +
  xlab("Major Family") +
  ylab("% of elements ") +
  ggtitle("Distribution by orientation group", 
          subtitle = expression("For inverted " * italic("Alu") * " cluster of size 3")) +
  guides(fill = guide_legend(title = "Family")) 


# Join --------------------------------------------------------------------
library(cowplot)
suppFig1 <- plot_grid(plot_grid(suppFig1a, suppFig1b, rel_widths = c(2, 1),
                            labels=c("A", "B"), ncol = 2, nrow = 1, align = 'hv', axis = "lb", label_size=18),
                  plot_grid(suppFig1c, suppFig1d, rel_widths = c(3, 1),
                            labels=c("C", "D"), ncol = 2, nrow = 1, align = 'hv', axis = "lb", label_size=18),
                  labels=c("", ""), ncol = 1, nrow = 2, align = 'hv', label_size=18)

save_plot(file.path(out_plots,"SuppFig2.1_AluExtra.pdf"), suppFig1, ncol = 1, nrow = 2, base_height = 8, base_width = 15)
save_plot(file.path(out_plots,"SuppFig2.1_AluExtra.png"), suppFig1, ncol = 1, nrow = 2, base_height = 8, base_width = 15)




