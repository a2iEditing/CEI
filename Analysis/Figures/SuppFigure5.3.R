rm(list = ls(all = TRUE))
# SETUP ===========================================================================================
# *************************************************************************************************
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
# Setup -------------------------------------------------------------------
### Output directory --------------------------------------------------------
outdir = "Analysis/Figures/"
out_plots="plots"
out_stats="stats"

input_dir = "Analysis/tables_TCGA/"

### create directories ------------------------------------------------------
#create paths, change
out_plots = file.path(outdir, out_plots)
out_stats = file.path(outdir, out_stats)
#create directories, if do not exist
dir.create(out_plots, recursive = T)
dir.create(out_stats, recursive = T)


# Load Data ---------------------------------------

#metadata
## SRA table ---------------------------------------
metadata = "TCGA/Metadata/TCGA_all_samples_genome_sample_sheet.tsv" %>% 
  fread(stringsAsFactors = F) %>% 
  select("File ID", "Project ID", "Case ID", "Sample ID", "Sample Type")  %>%
  mutate(Tissue = str_remove(`Project ID`, pattern = "TCGA-")) %>%
  rename("Sample" = `File ID`)

## TissueColors ------------------------------------------------------------
tcga_colors = "TCGA/Metadata/tcga_hex_color_codes.csv" %>% 
  fread() %>%
  filter(cancer_type %in% unique(metadata$Tissue)) %>%
  arrange(cancer_type) %>%
  pull(cancer_type_color)

tcga_colors_primary_tumor = "TCGA/Metadata/tcga_hex_color_codes.csv" %>% 
  fread() %>%
  filter(cancer_type %in% unique(metadata %>% filter(`Sample Type` == "Primary Tumor") %>% pull(Tissue))) %>%
  arrange(cancer_type) %>%
  pull(cancer_type_color)

tcga_colors_normal = "TCGA/Metadata/tcga_hex_color_codes.csv" %>% 
  fread() %>%
  filter(cancer_type %in% unique(metadata %>% filter(`Sample Type` == "Solid Tissue Normal") %>% pull(Tissue))) %>%
  arrange(cancer_type) %>%
  pull(cancer_type_color)

## index ----

cytoindex_data = fread("TCGA/Summary/TCGA_3UTREditingIndex.csv", stringsAsFactors = F)
global_index_data = fread("TCGA/Summary/TCGA_filtered_AluEditingIndex.csv", stringsAsFactors = F)

l_EI_wMETA = cytoindex_data %>%
  select(Sample, ends_with("EditingIndex")) %>%
  pivot_longer(cols = ends_with("EditingIndex"), 
               values_to = "EditingIndex", 
               names_to = "Mismatch", 
               names_pattern = "([ACGT]2[ACGT])") %>%
  mutate("Editing Index Type" = "Cytoplasmic Index",
         EI_type = "IRAlu3UTR") %>%
  bind_rows(global_index_data %>%
              select(Sample, ends_with("EditingIndex")) %>%
              pivot_longer(cols = ends_with("EditingIndex"), 
                           values_to = "EditingIndex", 
                           names_to = "Mismatch", 
                           names_pattern = "([ACGT]2[ACGT])") %>%
              mutate("Editing Index Type" = "Global Index",
                     EI_type = "GenomeWideAlu")) %>%
  inner_join(metadata)


allData = inner_join(cytoindex_data %>%
                       rename_with(.cols = -Sample, ~paste0(.x, "IRAlu3UTR")), 
                     global_index_data %>%
                       select(any_of(intersect(colnames(global_index_data), colnames(cytoindex_data)))) %>%
                       rename_with(.cols = -Sample, ~paste0(.x, "GenomeWideAlu"))) %>%
  left_join(l_EI_wMETA %>%
              filter(Mismatch != "A2G") %>%
              group_by(Sample, EI_type) %>%
              summarise(NextBestMM = max(EditingIndex, na.rm = T)) %>% 
              mutate(EI_type = paste0("NextBestMMAfterA2G_", EI_type)) %>%
              pivot_wider(names_from = EI_type, values_from = NextBestMM),
            by = "Sample",
            suffix = c("", "_NextBestMMA2G")) %>%
  mutate(Signal2NoiseRatio_GenomeWideAlu = A2GEditingIndexGenomeWideAlu / NextBestMMAfterA2G_GenomeWideAlu,
         Signal2NoiseRatio_IRAlu3UTR = A2GEditingIndexIRAlu3UTR / NextBestMMAfterA2G_IRAlu3UTR,
         # contribution to index
         ContributionToIndex_IRAlu3UTR = IndexedMismatchesOfA2GIRAlu3UTR / sum(IndexedMismatchesOfA2GIRAlu3UTR + IndexedCanonicalOfA2GIRAlu3UTR)) %>%
  inner_join(metadata)

l_EI_wMETA_a2g =  l_EI_wMETA %>%
  filter(Mismatch == "A2G") %>%
  inner_join(l_EI_wMETA %>%
               filter(Mismatch != "A2G") %>%
               group_by(Sample, EI_type) %>%
               summarise(NextBestMM = max(EditingIndex, na.rm = T))) %>% 
  mutate(Signal2NoiseRatio = EditingIndex / NextBestMM)


# Filter ------------------------------------------------------------------

allData_filtered = allData %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)


# Remove replicates -------------------------------------------------------
allData_filtered_mergedReps = allData_filtered %>% 
  group_by(`Project ID`, `Case ID`, `Sample Type`, Tissue) %>%
  summarise(across(A2CEditingIndexIRAlu3UTR:ContributionToIndex_IRAlu3UTR, mean),
            NumReps = n()) %>%
  ungroup() 

l_EI_wMETA_a2g_filtered_mergedReps = allData_filtered_mergedReps %>%
  select(`Project ID`:Tissue, contains("A2GEditingIndex")) %>%
  pivot_longer(cols = contains("A2GEditingIndex"), names_to = "EI_type", 
               values_to = "EditingIndex", names_prefix = "A2GEditingIndex") %>%
  mutate("Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "GenomeWideAlu", 
                                                    "CEI" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("CEI"))

# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")



# Supp. Fig. 3A -----------------------------------------------------------
correlation_coef_suppFig5a = coef(lm(allData$A2GEditingIndexIRAlu3UTR ~ 0 + allData$A2GEditingIndexGenomeWideAlu))
suppFig5a = allData %>%
  ggplot(aes(x = A2GEditingIndexGenomeWideAlu, y = A2GEditingIndexIRAlu3UTR)) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  theme_custom() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_abline(slope = correlation_coef_suppFig5a, intercept = 0, color = "#8B1A1A") +
  scale_x_continuous(breaks = c(0:3)) +
  scale_y_continuous(breaks = c(0:6)) +
  ggtitle("Correlation of editing signal",
          subtitle = "Dashed line at x = y") +
  xlab("AEI") +
  ylab("CEI") +
  labs(color = "Sample density") +
  theme(legend.text  = element_text(angle = 45, hjust = 1)) +
  # annotate slope
  annotate("text",x=2,y=0.5, size = 5, hjust = 0, label=(paste0("Slope==",signif(correlation_coef_suppFig5a, digits = 3))),parse=TRUE, color = chosen_colors[2])+ 
  # # annotate line
  # annotate("text",x=0.5,y=12.5, size = 5, hjust = 0,label="Line: x==y",parse=TRUE) +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 6), expand = FALSE, default = FALSE, clip = "on")


# Supp. Fig. 3B -----------------------------------------------------------
correlation_coef_suppFig5b = coef(lm(allData$Signal2NoiseRatio_IRAlu3UTR ~ 0 + allData$Signal2NoiseRatio_GenomeWideAlu))
suppFig5b = allData %>%
  ggplot(aes(x=Signal2NoiseRatio_GenomeWideAlu,
             y=Signal2NoiseRatio_IRAlu3UTR))+
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  theme_custom() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_abline(slope = correlation_coef_suppFig5b, intercept = 0, color = "#8B1A1A") +
  scale_x_continuous(breaks = seq(0, 14, 2)) +
  scale_y_continuous(breaks = seq(0, 150, 25)) +
  ggtitle("Correlation of signal-to-noise ratios",
          subtitle = "Dashed line at x = y") +
  xlab("AEI") +
  ylab("CEI") +
  labs(color = "Sample density") +
  theme(legend.text  = element_text(angle = 45, hjust = 1)) +
  # annotate slope (forct intercept to be 0)
  annotate("text",x=1,y=125, size = 5, hjust = 0, label=(paste0("Slope==",signif(correlation_coef_suppFig5b, digits = 3))),parse=TRUE, color = chosen_colors[2])+
  # # annotate line
  # annotate("text",x=0,y=180, size = 5, hjust = 0,label="Line: x==y",parse=TRUE)
  coord_cartesian(xlim = c(0, 14), ylim = c(0, 150), expand = FALSE, default = FALSE, clip = "on")


# Supp. Fig. 3C -----------------------------------------------------------
# noise
suppFig5c = l_EI_wMETA_a2g %>% 
  mutate(`Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "Global Index", 
                                                    "CEI" = "Cytoplasmic Index") %>% 
           forcats::fct_relevel("CEI")) %>%
  ggplot(aes(x = NextBestMM, color = `Editing Index Type`)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  ggtitle("Cumulative distribution of noise",
          subtitle = "Dashed line at x = 0.3") +
  scale_x_log10() +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  theme_custom(legend_position = "inside") +
  scale_color_manual(values = index_colors) +
  xlab("Noise (maximum non-A2G mismatch)") +
  ylab("Cumulative proportion") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position.inside = c(0.1, 0.9),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent",   colour = NA))

suppFig5c_inset = l_EI_wMETA_a2g %>% 
  ggplot(aes(x = NextBestMM, fill = `Editing Index Type`)) +
  geom_histogram(binwidth = 0.01, alpha = 0.6, position = "identity") +
  ggtitle("Distribution of noise",
          subtitle = "Binwidth = 0.01") +
  scale_x_log10() +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  theme_custom(legend_position = "none") +
  scale_fill_manual(values = index_colors) +
  xlab("Noise") +
  ylab("Sample Count")

# Supp. Fig. 3D -----------------------------------------------------------

# compute comparisons
two_groups_comparison = function(all_data, case_group, control_group, tissue, variable, test = wilcox.test, min_samples = 3, paired = F) {
  print(tissue)
  
  total_samples = all_data%>% filter(Tissue == tissue, 
                                     (`Sample Type` == case_group | `Sample Type` == control_group) )%>% nrow
  
  if (paired) {
    # if this data is paired - calculate the p-value
    
    # check if both groups exist
    if (length(all_data%>% filter(Tissue == tissue, `Sample Type` == case_group ) %>% pull(variable)) > 0 &
        length(all_data%>% filter(Tissue == tissue, `Sample Type` == control_group ) %>% pull(variable)) > 0) {
      
      # Deal with samples that don't have a pair
      curr_data = all_data %>%
        filter(Tissue == tissue, 
               (`Sample Type` == case_group | `Sample Type` == control_group) ) %>%
        # choose relevant columns
        select(`Sample Type`, `Case ID`, all_of(c(variable))) %>%
        # make wide format to see who has a pair
        pivot_wider(values_from = variable, names_from = `Sample Type`, id_cols = `Case ID`) %>%
        # filter rows where either case or control are NA (incomplete rows)
        filter(if_all(everything(), ~!is.na(.x)))
      
      
      # get values for each group, which are in matched order as they are 2 columns in a dataframe
      case_values = curr_data %>% pull(case_group)
      control_values = curr_data %>% pull(control_group)
      
      # enough to check one group as they are equal in length
      if (length(case_values) >= min_samples) {
        pval = test(case_values, control_values, paired = T)$p.value } 
      else {
        pval = NA } 
    } else { 
      return(as.data.frame(list(Tissue = tissue,
                                TestedCase = case_group,
                                TestedControl = control_group,
                                Variable = variable,
                                MeanTestedCase = mean(all_data%>% filter(Tissue == tissue, `Sample Type` == case_group ) %>% pull(variable)),
                                MeanTestedControl = mean(all_data%>% filter(Tissue == tissue, `Sample Type` == control_group) %>% pull(variable)),
                                MedianTestedCase = median(all_data%>% filter(Tissue == tissue, `Sample Type` == case_group ) %>% pull(variable)),
                                MedianTestedControl = median(all_data%>% filter(Tissue == tissue, `Sample Type` == control_group ) %>% pull(variable)),
                                StdTestedCase = sd(all_data%>% filter(Tissue == tissue, `Sample Type` == case_group ) %>% pull(variable)),
                                StdTestedControl = sd(all_data%>% filter(Tissue == tissue, `Sample Type` == control_group ) %>% pull(variable)),
                                StdBothGroups = NA,
                                NumberOfSamplesCase = all_data%>% filter(Tissue == tissue, `Sample Type` == case_group ) %>% pull(variable) %>% length,
                                NumberOfSamplesControl = all_data%>% filter(Tissue == tissue, `Sample Type` == control_group ) %>% pull(variable) %>% length,
                                TotalSamples = total_samples,
                                PValue = NA))) }
  } else {
    # unpaired data - calculate the p-value
    
    # get values by indexes
    case_values = all_data%>% filter(Tissue == tissue, `Sample Type` == case_group ) %>% pull(variable)
    control_values = all_data %>% filter(Tissue == tissue, `Sample Type` == control_group) %>% pull(variable)
    
    # if observation number allows - calculate p-value
    if (length(case_values) >= min_samples & length(control_values) >= min_samples) {
      pval = test(case_values, control_values)$p.value } 
    else  {
      pval = NA }
  }
  
  return(as.data.frame(list(Tissue = tissue,
                            TestedCase = case_group,
                            TestedControl = control_group,
                            Variable = variable,
                            MeanTestedCase = mean(case_values, na.rm = T),
                            MeanTestedControl = mean(control_values, na.rm = T),
                            MedianTestedCase = median(case_values, na.rm = T),
                            MedianTestedControl = median(control_values, na.rm = T),
                            StdTestedCase = sd(case_values, na.rm = T),
                            StdTestedControl = sd(control_values, na.rm = T),
                            StdBothGroups = sd(c(case_values, control_values), na.rm = T),
                            NumberOfSamplesCase = length(case_values),
                            NumberOfSamplesControl = length(control_values),
                            TotalSamples = total_samples,
                            PValue = pval)))
}



### Primary vs. normal 
comparison_results = map_dfr(.x = allData$Tissue %>% unique,
                             ~ two_groups_comparison(all_data = allData_filtered_mergedReps,
                                                     case_group = "Primary Tumor",
                                                     control_group = "Solid Tissue Normal",
                                                     tissue = .x,
                                                     variable = "A2GEditingIndexGenomeWideAlu")) %>%
  rbind(map_dfr(.x = allData$Tissue %>% unique,
                ~ two_groups_comparison(all_data = allData_filtered_mergedReps,
                                        case_group = "Primary Tumor",
                                        control_group = "Solid Tissue Normal",
                                        tissue = .x,
                                        variable = "A2GEditingIndexIRAlu3UTR")) ) %>%
  # filter comparisons 
  # Rule 1: filter out comparisons where there are less than 4 samples for any of the groups 
  #         - implemented in function - p-value will be NA
  # Rule 2: filter out datasets where there are less than 10 samples in total
  # filter(TotalSamples >= 10, 
  #        !is.na(PValue),
  #        NumberOfSamplesCase > 3, NumberOfSamplesControl > 3) %>%
  mutate(CohensDEffectSize = (MeanTestedCase - MeanTestedControl) / StdBothGroups)


suppFig5d = comparison_results %>%
  # filter(TotalSamples >= 10,
  #        !is.na(PValue),
  #        NumberOfSamplesCase > 3, NumberOfSamplesControl > 3) %>%
  mutate(Significant = PValue < 0.05) %>%
  select(Tissue, CohensDEffectSize, Variable, Significant) %>%
  pivot_wider(values_from = c(CohensDEffectSize, Significant), names_from = Variable) %>%
  mutate(Significant = case_when(Significant_A2GEditingIndexGenomeWideAlu & Significant_A2GEditingIndexIRAlu3UTR ~ "Both Comparisons", 
                                 Significant_A2GEditingIndexGenomeWideAlu ~ "AEI Comparison",
                                 Significant_A2GEditingIndexIRAlu3UTR ~ "CEI Comparison",
                                 .default = "None") %>%
           factor(levels = c("Both Comparisons", "CEI Comparison", "AEI Comparison", "None"))) %>%
  ggplot(aes(x = CohensDEffectSize_A2GEditingIndexGenomeWideAlu, y = CohensDEffectSize_A2GEditingIndexIRAlu3UTR, color = Significant)) +
  geom_point(size = 4) +
  theme_custom() +
  ggtitle("Effect size of case-control comparison per editing index type across datasets",
          subtitle = "Green triangle marks positive (increasing) effect size greater in CEI") +
  # Cohen's D = (Mean editing index case - mean editing index control) / editing index sd
  scale_color_manual(values = c("#8B1A1A", "#e37e00", "#1a476f", "grey60")) +
  expand_limits(y = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize)),
                x = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize))) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept =  0 , linetype = "dashed", color = "grey") +
  geom_vline(xintercept =  0 , linetype = "dashed", color = "grey") +
  xlab("AEI effect size (Cohen's D)") +
  ylab("CEI effect size (Cohen's D)") +
  # draw triangles
  geom_segment(data = data.frame(x1 = 0, x2 = 0, y1 = 0, y2 = 2),
               mapping = aes(x = x1, y = y1, xend = x2, yend = y2), 
               colour = "#6E8B3D", linewidth = 1.5) +
  geom_segment(data = data.frame(x1 = 0, x2 = 2, y1 = 0, y2 = 2),
               mapping = aes(x = x1, y = y1, xend = x2, yend = y2), 
               colour = "#6E8B3D", linewidth = 1.5) +
  geom_segment(data = data.frame(x1 = 0, x2 = 2, y1 = 2, y2 = 2),
               mapping = aes(x = x1, y = y1, xend = x2, yend = y2), 
               colour = "#6E8B3D", linewidth = 1.5) +
  coord_cartesian(xlim = c(-2, 2.01), ylim = c(-2, 2.01), expand = FALSE, default = FALSE, clip = "on")


# Supp. Fig. 3E -----------------------------------------------------------
suppFig5e = comparison_results %>%
  mutate(Significant = PValue < 0.05) %>%
  select(Tissue, CohensDEffectSize, Variable, Significant) %>%
  pivot_wider(values_from = c(CohensDEffectSize, Significant), names_from = Variable) %>%
  filter(!is.na(CohensDEffectSize_A2GEditingIndexGenomeWideAlu)) %>%
  ggplot(aes(x = CohensDEffectSize_A2GEditingIndexGenomeWideAlu, y = CohensDEffectSize_A2GEditingIndexIRAlu3UTR, fill = Tissue)) +
  geom_point(size = 4, shape = 21) +
  theme_custom() +
  ggtitle("Effect size of case-control comparison per editing index type across datasets") +
  # Cohen's D = (Mean editing index case - mean editing index control) / editing index sd
  # scale_color_manual(values = c("#8B1A1A", "#e37e00", "#1a476f", "grey60")) +
  scale_fill_manual(values = tcga_colors_primary_tumor) +
  expand_limits(y = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize)),
                x = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize))) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept =  0 , linetype = "dashed", color = "grey") +
  geom_vline(xintercept =  0 , linetype = "dashed", color = "grey") +
  xlab("AEI effect size (Cohen's D)") +
  ylab("CEI effect size (Cohen's D)") +
  guides(fill = guide_legend(nrow = 3)) +
  # draw triangles
  # geom_segment(data = data.frame(x1 = 0, x2 = 0, y1 = 0, y2 = 2),
  #              mapping = aes(x = x1, y = y1, xend = x2, yend = y2), 
  #              colour = "#6E8B3D", linewidth = 1.5) +
  # geom_segment(data = data.frame(x1 = 0, x2 = 2, y1 = 0, y2 = 2),
  #              mapping = aes(x = x1, y = y1, xend = x2, yend = y2), 
  #              colour = "#6E8B3D", linewidth = 1.5) +
  # geom_segment(data = data.frame(x1 = 0, x2 = 2, y1 = 2, y2 = 2),
  #              mapping = aes(x = x1, y = y1, xend = x2, yend = y2), 
  #              colour = "#6E8B3D", linewidth = 1.5) +
  coord_cartesian(xlim = c(-2, 2.01), ylim = c(-2, 2.01), expand = FALSE, default = FALSE, clip = "on")+
  theme(legend.title = element_blank())


# Supp. Fig. 3F -------------------------------------------------------------

# only primary vs. normal
suppFig5f = l_EI_wMETA_a2g_filtered_mergedReps %>%
  filter(`Sample Type` == "Primary Tumor" | `Sample Type` == "Solid Tissue Normal") %>%
  mutate(`Sample Type` = factor(`Sample Type`, levels = c("Solid Tissue Normal", "Primary Tumor"))) %>%
  ggplot(aes(x = `Editing Index Type`, y = EditingIndex, fill = `Sample Type`)) +
  geom_boxplot()+
  # geom_violin(draw_quantiles = 0.5, trim = T) +
  facet_grid(.~Tissue) +
  ggtitle("A2G Editing Index in TCGA",
          subtitle = "Primary tumor and normal samples only") +
  scale_fill_manual(values = case_control_colors) +
  theme_custom() +
  theme(axis.text.x=element_text(angle=50, hjust=1, size = 10))  + 
  expand_limits(y = 0)+
  theme(legend.title = element_blank())


# Join --------------------------------------------------------------------
library(cowplot)
suppFig5 <- plot_grid(plot_grid(suppFig5a, 
                                suppFig5b, 
                                suppFig5c + 
                                  # inset
                                  patchwork::inset_element(suppFig5c_inset,
                                                           left   = 0.67, right = 0.98,
                                                           bottom = 0.05, top   = 0.6), 
                                rel_widths = c(1, 1, 2), 
                            labels=c("A", "B", "C"), ncol = 3, nrow = 1, align = 'hv', label_size=18),
                  plot_grid(suppFig5d, suppFig5e, 
                            labels=c("D", "E"), ncol = 2, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  plot_grid(suppFig5f, 
                            labels=c("F"), ncol = 1, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  rel_heights = c(0.8, 1, 0.7),
                  labels=c("", "", ""), ncol = 1, nrow = 3, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"SuppFig5.3_TCGA.pdf"), suppFig5, ncol = 1, nrow = 3, base_height = 10, base_width = 20)
save_plot(file.path(out_plots,"SuppFig5.3_TCGA.png"), suppFig5, ncol = 1, nrow = 3, base_height = 10, base_width = 20)

# ************************************************************