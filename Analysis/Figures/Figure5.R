
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


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")


# DATA: INFECTIOUS DISEASES -----------------------------------------------------


# Load intermediate -------------------------------------------------------
comparisons = fread("Analysis/tables/Classifications.Metadata.comparisons.new.v2.nonewlines.csv", stringsAsFactors = F)

pathogen_disease_classification = fread("Analysis/tables/pathogen_disease_classification_withBioProjects.csv", stringsAsFactors = F) %>%
  separate_longer_delim(cols = c(BioProject, Case), delim = ";")

l_EI_wMETA_infections = fread(file.path(input_dir, "EditingIndex.longFormat.wMetadata.csv"), stringsAsFactors = F) %>%
  mutate(disease_or_control = case_when(Disease_Status_manual_annotation == "Control" ~ "Control",
                                        Disease_Status_manual_annotation == "OtherDesign" ~ "OtherDesign",
                                        .default = "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case", "OtherDesign")),
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "Global Index" = "GenomeWideAlu", 
                                                    "Cytoplasmic Index" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("Cytoplasmic Index"))

allData_infections = fread(file.path(input_dir, "AllData.wMetadata.csv"), stringsAsFactors = F) %>%
  mutate(disease_or_control = case_when(Disease_Status_manual_annotation == "Control" ~ "Control",
                                        Disease_Status_manual_annotation == "OtherDesign" ~ "OtherDesign",
                                        .default = "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case", "OtherDesign"))) %>%
  # add next best mismatch
  left_join(l_EI_wMETA_infections %>%
              filter(Mismatch != "A2G") %>%
              group_by(Sample, EI_type) %>%
              summarise(NextBestMM = max(EditingIndex, na.rm = T)) %>% 
              mutate(EI_type = paste0("NextBestMMAfterA2G_", EI_type)) %>%
              pivot_wider(names_from = EI_type, values_from = NextBestMM),
            by = "Sample",
            suffix = c("", "_NextBestMMA2G")) %>%
  mutate(Signal2NoiseRatio_GenomeWideAlu = A2GEditingIndexGenomeWideAlu / NextBestMMAfterA2G_GenomeWideAlu,
         Signal2NoiseRatio_IRAlu3UTR = A2GEditingIndexIRAlu3UTR / NextBestMMAfterA2G_IRAlu3UTR,
         "Uniquely mapped Mbp" = (`Uniquely mapped reads number`/10^6)*`Average mapped length`)

l_EI_wMETA_infections_a2g =  l_EI_wMETA_infections %>%
  filter(Mismatch == "A2G") %>%
  inner_join(l_EI_wMETA_infections %>%
               filter(Mismatch != "A2G") %>%
               group_by(Sample, EI_type) %>%
               summarise(NextBestMM = max(EditingIndex, na.rm = T))) %>% 
  mutate(Signal2NoiseRatio = EditingIndex / NextBestMM)

# case-control
allData_infections_case_control = allData_infections %>%
  filter(StudyDesign=="Case-Control")

l_EI_wMETA_infections_case_control = l_EI_wMETA_infections %>%
  filter(StudyDesign=="Case-Control")

l_EI_wMETA_infections_a2g_case_control = l_EI_wMETA_infections_a2g %>%
  filter(StudyDesign=="Case-Control")


# Filter ------------------------------------------------------------------

allData_infections_filtered = allData_infections %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)

allData_infections_case_control_filtered = allData_infections_case_control %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)

l_EI_wMETA_infections_a2g_filtered = l_EI_wMETA_infections_a2g %>%
  filter(Sample %in% allData_infections_filtered$Sample)

l_EI_wMETA_infections_a2g_case_control_filtered = l_EI_wMETA_infections_a2g_case_control %>%
  filter(Sample %in% allData_infections_case_control_filtered$Sample)

# Fig. 5A ------------------------------------------------------------------
# compute comparisons
two_groups_comparison = function(all_data, comparison_line, variable, test = wilcox.test, min_samples = 3) {
  print(comparison_line[["BioProject"]])
  
  # get current data
  curr_data = all_data %>%
    filter(BioProject == comparison_line[["BioProject"]],
           Disease_Status_manual_annotation == comparison_line[["Case"]] | Disease_Status_manual_annotation == comparison_line[["Control"]])
  
  
  # if there are two read lengths - use recursion to do this for each
  if (curr_data$ReadLength %>% unique %>% length > 1) {
    readlengths = curr_data$ReadLength %>% unique
    d = map_dfr(.x = readlengths, 
                .f = ~ two_groups_comparison(all_data = curr_data %>% 
                                               filter(ReadLength == .x),
                                             comparison_line = comparison_line,
                                             variable = variable) %>%
                  mutate(ReadLength = .x))
    return(d)
  }
  
  
  # if this has different tissues - use recursion to recalculate this per-tissue
  # add tissues to case-control annotation to allow for comparison
  if (comparison_line[["MultipleTissuesData"]] == "Yes") {
    
    tissues = curr_data %>% pull(multiple_tissue_id) %>% unique()
    comparison_line[["MultipleTissuesData"]] = "No"
    d = map_dfr(.x = tissues, 
                .f = ~ two_groups_comparison(all_data = curr_data %>% 
                                               filter(multiple_tissue_id == .x),
                                             comparison_line = comparison_line,
                                             variable = variable) %>%
                  mutate(Tissue = .x))
    return(d)
  }
  
  
  # # remove replicates if required, group also by paired data just in case
  if (comparison_line[["TechicalReplicates"]] == "Yes") {
    curr_data = curr_data %>%
      # group and create mean per technical replicate
      group_by(Disease_Status_manual_annotation, technical_repeats_id, paired_sample_id) %>%
      summarise(across(all_of(c(variable)), mean)) %>%
      # release grouping
      ungroup()
  }
  
  
  if (comparison_line[["PairedData"]] == "Yes") {
    # if this data is paired - calculate the p-value
    
    # keep total number of samples before filtering
    total_samples = curr_data%>%nrow
    
    # Deal with samples that don't have a pair
    curr_data = curr_data %>%
      # choose relevant columns
      select(Disease_Status_manual_annotation, paired_sample_id, all_of(c(variable))) %>%
      # make wide format to see who has a pair
      pivot_wider(values_from = variable, names_from = Disease_Status_manual_annotation) %>%
      # filter rows where either case or control are NA (incomplete rows)
      filter(if_all(everything(), ~!is.na(.x)))
    
    # get values for each group, which are in matched order as they are 2 columns in a dataframe
    case_values = curr_data %>% pull(comparison_line[["Case"]])
    control_values = curr_data %>% pull(comparison_line[["Control"]])
    
    # enough to check one group as they are equal in length
    if (length(case_values) >= min_samples) {
      pval = test(case_values, control_values, paired = T)$p.value } 
    else {
      pval = NA }
    
  } else {
    # unpaired data - calculate the p-value
    
    # number of samples is as is
    total_samples = curr_data%>% nrow
    
    # get values by indexes
    case_values = curr_data %>% filter(Disease_Status_manual_annotation == comparison_line[["Case"]]) %>% pull(variable)
    control_values = curr_data %>% filter(Disease_Status_manual_annotation == comparison_line[["Control"]]) %>% pull(variable)
    
    # if observation number allows - calculate p-value
    if (length(case_values) >= min_samples & length(control_values) >= min_samples) {
      pval = test(case_values, control_values)$p.value } 
    else  {
      pval = NA }
  }
  
  return(as.data.frame(list(BioProject = comparison_line[["BioProject"]],
                            TestedCase = comparison_line[["Case"]],
                            TestedControl = comparison_line[["Control"]],
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

comparison_results = map_dfr(.x = c(1:nrow(comparisons)),
                             ~ two_groups_comparison(all_data = allData_infections_case_control_filtered,
                                                     comparison_line = comparisons[.x,] %>% as.vector,
                                                     variable = "A2GEditingIndexGenomeWideAlu")) %>%
  rbind(map_dfr(.x = c(1:nrow(comparisons)),
                ~ two_groups_comparison(all_data = allData_infections_case_control_filtered,
                                        comparison_line = comparisons[.x,] %>% as.vector,
                                        variable = "A2GEditingIndexIRAlu3UTR")) ) %>%
  mutate(Comparison = paste(TestedCase, replace_na(Tissue, ""), replace_na(as.character(ReadLength), ""), BioProject)) %>%
  # filter comparisons 
  # Rule 1: filter out comparisons where there are less than 4 samples for any of the groups 
  #         - implemented in function - p-value will be NA
  # Rule 2: filter out datasets where there are less than 10 samples in total
  filter(NumberOfSamplesCase+NumberOfSamplesControl >= 10, 
         !is.na(PValue),
         NumberOfSamplesCase > 3, NumberOfSamplesControl > 3) %>%
  mutate(CohensDEffectSize = (MeanTestedCase - MeanTestedControl) / StdBothGroups)


fig5a =  ggplot(data.frame(x = 1:1,
                           colors = c("white"),
                           text = "TODO: add\nfig5A from biorender"), 
                aes(x, y = 0, fill = colors, label = text)) +
  geom_tile(width = .9, height = .9) + # make square tiles
  geom_text(color = "Red", size = 7) + # add white text in the middle
  scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
  coord_fixed() + # make sure tiles are square
  theme_void() # remove any axis markings

# Fig. 5B ------------------------------------------------------------------

fig5b = comparison_results %>%
  mutate(Significant = PValue < 0.05) %>%
  select(BioProject, Comparison, TestedCase, TestedControl, CohensDEffectSize, Variable, Significant) %>%
  pivot_wider(values_from = c(CohensDEffectSize, Significant), names_from = Variable) %>%
  mutate(Significant = case_when(Significant_A2GEditingIndexGenomeWideAlu & Significant_A2GEditingIndexIRAlu3UTR ~ "Both Comparisons", 
                                 Significant_A2GEditingIndexGenomeWideAlu ~ "AEI Comparison",
                                 Significant_A2GEditingIndexIRAlu3UTR ~ "CEI Comparison",
                                 .default = "None") %>%
           factor(levels = c("Both Comparisons", "CEI Comparison", "AEI Comparison", "None"))) %>%
  ggplot(aes(x = CohensDEffectSize_A2GEditingIndexGenomeWideAlu, y = CohensDEffectSize_A2GEditingIndexIRAlu3UTR, color = Significant)) +
  geom_point(size = 3) +
  theme_custom() +
  ggtitle("Effect size of case-control comparison per editing index type across datasets",
          subtitle = "Green triangle marks positive (increasing) effect size greater in CEI") +
  theme(plot.title = element_text(vjust = 5),
        plot.subtitle = element_text(vjust = 4.5)) +
  # Cohen's D = (Mean editing index case - mean editing index control) / editing index sd
  scale_color_manual(values = c("#8B1A1A", "#e37e00", "#1a476f", "grey60")) +
  expand_limits(y = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize)),
                x = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize))) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept =  0 , linetype = "dashed", color = "grey") +
  geom_vline(xintercept =  0 , linetype = "dashed", color = "grey") +
  xlab("AEI effect size (Cohen's d)") +
  ylab("CEI effect size (Cohen's d)") +
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


# Fig. 5C ------------------------------------------------------------------

fig5c = allData_infections_filtered %>%
  filter(BioProject == "PRJNA352062") %>% 
  select(Sample, A2GEditingIndexGenomeWideAlu, A2GEditingIndexIRAlu3UTR, technical_repeats_id) %>% 
  inner_join(fread("Summary/Metadata/PRJNA352062.csv", stringsAsFactors = F) %>%
               select(Run, source_name, disease_state, timetonegativity, Time, treatmentresult),
             join_by("Sample" == "Run")) %>%
  # filter to wanted group
  filter(treatmentresult == "Definite Cure", disease_state == "TB Subjects") %>%
  mutate(Donor = str_extract(technical_repeats_id, pattern = "S[0-9]+"), 
         Timepoint = str_extract(technical_repeats_id, pattern = "_.*") %>% 
           str_remove("^_") %>% str_replace(pattern = "_", replacement = " ") %>% 
           str_to_title() %>%
           forcats::fct_recode(`Treatment Start` = "Dx") %>%
           forcats::fct_relevel(c("Treatment Start", "Day 7", "Week 4", "Week 24"))) %>%
  group_by(Donor, Timepoint, timetonegativity) %>%
  summarise(across(c(A2GEditingIndexGenomeWideAlu, A2GEditingIndexIRAlu3UTR), mean)) %>%
  pivot_longer(cols = c(A2GEditingIndexGenomeWideAlu, A2GEditingIndexIRAlu3UTR), names_to = "Editing Index Type", values_to = "Editing Index") %>%
  mutate(`Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "A2GEditingIndexGenomeWideAlu", 
                                                    "CEI" = "A2GEditingIndexIRAlu3UTR") %>% 
           forcats::fct_relevel("CEI")) %>%
  mutate(timetonegativity = forcats::fct_collapse(as.factor(timetonegativity), 
                                                  `Week 8-12` = c("Week08", "Week12"),
                                                  `Week 4` = "Week04",
                                                  `Week 24` = "Week24")) %>%
  # filter out cured at week 24
  filter(timetonegativity != "Week 24") %>%
  
  ggplot(aes(x = Timepoint, y = `Editing Index`, fill = Timepoint)) +
  geom_boxplot(outliers = F) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_grid(`Editing Index Type`~timetonegativity, scales = "free", space = "free") +
  # geom_jitter(width = 0.05, alpha = 0.5) +
  theme_custom(legend_position = "none") +
  scale_y_continuous(breaks = seq(0, 4, 0.5)) +
  ggtitle("Decrease in cytoplasmic editing in tuberculosis after treatment and recovery",
          subtitle = "Longitudinal study grouped by recovery time, PRJNA352062") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, vjust = 1)) +
  expand_limits(y = 0) +
  # expand_limits(y = 1.05) +
  scale_fill_manual(values = c(tail(chosen_colors, -1), head(chosen_colors, 1))) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format", comparisons = list(c(1,2),
                                                                                            c(2,3),
                                                                                            c(3,4)),
                             label.y = 0.5,
                             step.increase = -0.1,
                             tip.length = -0.02,
                             vjust = 3)


# DATA: TCGA --------------------------------------------------------------------

# Load Data ---------------------------------------

#metadata_tcga
## SRA table ---------------------------------------
metadata_tcga = "TCGA/Metadata/TCGA_all_samples_genome_sample_sheet.tsv" %>% 
  fread(stringsAsFactors = F) %>% 
  select("File ID", "Project ID", "Case ID", "Sample ID", "Sample Type")  %>%
  mutate(Tissue = str_remove(`Project ID`, pattern = "TCGA-")) %>%
  rename("Sample" = `File ID`)

## survival ---------------
# see TCGA_KaplanMeyer_analysis.bestPick.R
survival_analysis_tcga = fread("Analysis/tables_TCGA/KaplanMeier_Combinations_SurvivalAnalysis.PrimaryTumor.PositiveSurvival.csv") %>%
  rename(Tissue = tissue)

## TissueColors ------------------------------------------------------------
tcga_colors = "TCGA/Metadata/tcga_hex_color_codes.csv" %>% 
  fread() %>%
  { setNames(pull(., cancer_type_color), pull(., cancer_type)) }

tcga_colors_primary_tumor = "TCGA/Metadata/tcga_hex_color_codes.csv" %>% 
  fread() %>%
  filter(cancer_type %in% unique(metadata_tcga %>% filter(`Sample Type` == "Primary Tumor") %>% pull(Tissue))) %>%
  arrange(cancer_type) %>%
  pull(cancer_type_color)

tcga_colors_normal = "TCGA/Metadata/tcga_hex_color_codes.csv" %>% 
  fread() %>%
  filter(cancer_type %in% unique(metadata_tcga %>% filter(`Sample Type` == "Solid Tissue Normal") %>% pull(Tissue))) %>%
  arrange(cancer_type) %>%
  pull(cancer_type_color)

## index ----

cytoindex_data_tcga = fread("TCGA/Summary/TCGA_3UTREditingIndex.csv", stringsAsFactors = F)
global_index_data_tcga = fread("TCGA/Summary/TCGA_filtered_AluEditingIndex.csv", stringsAsFactors = F)

l_EI_wMETA_tcga = cytoindex_data_tcga %>%
  select(Sample, ends_with("EditingIndex")) %>%
  pivot_longer(cols = ends_with("EditingIndex"), 
               values_to = "EditingIndex", 
               names_to = "Mismatch", 
               names_pattern = "([ACGT]2[ACGT])") %>%
  mutate("Editing Index Type" = "Cytoplasmic Index",
         EI_type = "IRAlu3UTR") %>%
  bind_rows(global_index_data_tcga %>%
              select(Sample, ends_with("EditingIndex")) %>%
              pivot_longer(cols = ends_with("EditingIndex"), 
                           values_to = "EditingIndex", 
                           names_to = "Mismatch", 
                           names_pattern = "([ACGT]2[ACGT])") %>%
              mutate("Editing Index Type" = "Global Index",
                     EI_type = "GenomeWideAlu")) %>%
  inner_join(metadata_tcga)


allData_tcga = inner_join(cytoindex_data_tcga %>%
                            rename_with(.cols = -Sample, ~paste0(.x, "IRAlu3UTR")), 
                          global_index_data_tcga %>%
                            select(any_of(intersect(colnames(global_index_data_tcga), colnames(cytoindex_data_tcga)))) %>%
                            rename_with(.cols = -Sample, ~paste0(.x, "GenomeWideAlu"))) %>%
  left_join(l_EI_wMETA_tcga %>%
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
  inner_join(metadata_tcga)

l_EI_wMETA_tcga_a2g =  l_EI_wMETA_tcga %>%
  filter(Mismatch == "A2G") %>%
  inner_join(l_EI_wMETA_tcga %>%
               filter(Mismatch != "A2G") %>%
               group_by(Sample, EI_type) %>%
               summarise(NextBestMM = max(EditingIndex, na.rm = T))) %>% 
  mutate(Signal2NoiseRatio = EditingIndex / NextBestMM)


# Filter ------------------------------------------------------------------

allData_tcga_filtered = allData_tcga %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)


# Remove replicates -------------------------------------------------------
allData_tcga_filtered_mergedReps = allData_tcga_filtered %>% 
  group_by(`Project ID`, `Case ID`, `Sample Type`, Tissue) %>%
  summarise(across(A2CEditingIndexIRAlu3UTR:ContributionToIndex_IRAlu3UTR, mean),
            NumReps = n()) %>%
  ungroup() 

l_EI_wMETA_tcga_a2g_filtered_mergedReps = allData_tcga_filtered_mergedReps %>%
  select(`Project ID`:Tissue, contains("A2GEditingIndex")) %>%
  pivot_longer(cols = contains("A2GEditingIndex"), names_to = "EI_type", 
               values_to = "EditingIndex", names_prefix = "A2GEditingIndex") %>%
  mutate("Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "Global" = "GenomeWideAlu", 
                                                    "Cytoplasmic" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("Cytoplasmic"))

# Fig. 5D -----------------------------------------------------------------
fig5d = l_EI_wMETA_tcga_a2g_filtered_mergedReps %>%
  filter(`Sample Type` == "Primary Tumor" | `Sample Type` == "Solid Tissue Normal") %>%
  mutate(`Sample Type` = factor(`Sample Type`, levels = c("Solid Tissue Normal", "Primary Tumor"))) %>%
  group_by(`Editing Index Type`, `Sample Type`, Tissue) %>%
  summarise(across(c(EditingIndex), list(sd = sd,
                                         mean = mean))) %>% 
  mutate(`Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "Global", 
                                                    "CEI" = "Cytoplasmic") %>% 
           forcats::fct_relevel("CEI")) %>%
  mutate("Normalized standard deviation" = EditingIndex_sd / EditingIndex_mean) %>%
  ggplot(aes(x = `Editing Index Type`, y = `Normalized standard deviation`, fill = `Editing Index Type`)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(size = 2, width = 0.05) +
  theme_custom(legend_position = "none") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12)) +
  facet_grid(.~`Sample Type`)+
  ggtitle("Normalized standard deviation", subtitle = "Per tissue") +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c(index_colors, case_control_colors, "grey")) +
  ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2)), paired = T)

# Fig. 5E -----------------------------------------------------------------
# Survival 

# filter by minimum per tissue
n_threshold_tissue = 100
fig5e = survival_analysis_tcga %>% 
  # filter to include at least 100 of each tissue
  filter(n_total >= n_threshold_tissue)  %>%
  group_by(variable) %>%
  mutate(adjusted_pval = p.adjust(pval, "fdr")) %>%
  group_by(Tissue, variable) %>%
  # pick best p-value
  slice_min(order_by = pval, n = 1) %>%mutate("Editing Index Type" = variable %>% as.factor(),
                                              `Editing Index Type` = forcats::fct_recode(`Editing Index Type`,
                                                                                         "AEI" = "A2GEditingIndexGenomeWideAlu",
                                                                                         "CEI" = "A2GEditingIndexIRAlu3UTR")) %>%
  # mutate(adjusted_pval2 = -1*log10(adjusted_pval)) %>%
  pivot_wider(names_from = `Editing Index Type`, values_from = adjusted_pval, id_cols = Tissue) %>%
  mutate(significant = AEI<0.05 | CEI<0.05) %>%
  ggplot(aes(x = AEI, y = CEI, fill = Tissue, alpha = significant)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(size = 3, color = "black", shape = 21) +
  theme_custom() +
  ggtitle("Adjusted P-value of Kaplan-Meier for best possible data division", subtitle = "Primary tumors with at least 100 samples") +
  scale_fill_manual(values = tcga_colors) +
  scale_alpha_manual(values = c(0.4, 1), guide = "none") +
  # scale_x_continuous(breaks = c(0, 0.05, seq(0.1, 0.85, 0.1))) +
  # scale_y_continuous(breaks = c(0, 0.05, seq(0.1, 0.85, 0.1))) +
  scale_x_log10(
    breaks = c(1e-2, 5e-2, 1e-1, 5e-1, 1),
    labels = c("0.01", "0.05", "0.1", "0.5", "1")
  ) +
  scale_y_log10(
    breaks = c(1e-2, 5e-2, 1e-1, 5e-1, 1),
    labels = c("0.01", "0.05", "0.1", "0.5", "1")
  ) +
  # geom_abline(slope = 1, intercept = 0.03, linetype = "dashed", color = "grey") +
  # geom_abline(slope = 1, intercept = -0.03, linetype = "dashed", color = "grey") +
  geom_hline(yintercept =  0.05, linetype = "dashed", color = "grey") +
  geom_vline(xintercept =  0.05 , linetype = "dashed", color = "grey") +
  expand_limits(y = 0, x = 0) +
  # coord_cartesian(xlim = c(0, 0.75), ylim = c(0, 0.75), default = FALSE, clip = "on") +
  guides(fill = guide_legend(nrow = 2)) +
  theme(legend.title = element_blank()) 


# DATA: GTEX --------------------------------------------------------------

# Load intermediate -------------------------------------------------------
metadata_gtex = fread("GTEx/Metadata/SraRunTable_GTEX_RNAseq.filtered_noTumor_noStranded_paired.minimal.csv") %>%
  mutate(CohortOriginal = Cohort, 
         Cohort = if_else(Cohort == "Postmortem", Cohort, "Antemortem")) 

gtex_colors = fread("GTEx/Metadata/gtex_colors.csv") %>%
  { setNames(pull(., Color), pull(., Tissue)) }


allData_gtex = fread("GTEx/Summary/AluEditingIndex.csv", stringsAsFactors = F) %>%
  inner_join(fread("GTEx/Summary/UTR3EditingIndex.csv", stringsAsFactors = F),
             by = "Sample", suffix = c("GenomeWideAlu", "IRAlu3UTR")) 


l_EI_wMETA_gtex = allData_gtex %>%
  select(Sample, contains("EditingIndex")) %>%
  pivot_longer(cols = contains("EditingIndex"), values_to = "EditingIndex", names_to = c("Mismatch", "EI_type"), names_sep = "EditingIndex") %>%
  mutate(`Editing Index Type` = case_match(EI_type, 
                                           "GenomeWideAlu" ~ "Global Index",
                                           "IRAlu3UTR" ~ "Cytoplasmic Index") %>%
           forcats::fct_relevel("Cytoplasmic Index")) %>%
  inner_join(metadata_gtex, by = join_by(Sample == "Run"))

allData_gtex = allData_gtex %>%
  # add next best mismatch
  left_join(l_EI_wMETA_gtex %>%
              filter(Mismatch != "A2G") %>%
              group_by(Sample, EI_type) %>%
              summarise(NextBestMM = max(EditingIndex, na.rm = T)) %>% 
              mutate(EI_type = paste0("NextBestMMAfterA2G_", EI_type)) %>%
              pivot_wider(names_from = EI_type, values_from = NextBestMM),
            by = "Sample",
            suffix = c("", "_NextBestMMA2G")) %>%
  inner_join(metadata_gtex, by = join_by(Sample == "Run"))

l_EI_wMETA_gtex_a2g =  l_EI_wMETA_gtex %>%
  filter(Mismatch == "A2G") %>%
  inner_join(l_EI_wMETA_gtex %>%
               filter(Mismatch != "A2G") %>%
               group_by(Sample, EI_type) %>%
               summarise(NextBestMM = max(EditingIndex, na.rm = T))) %>% 
  mutate(Signal2NoiseRatio = EditingIndex / NextBestMM)



# Filter ------------------------------------------------------------------

allData_gtex_filtered = allData_gtex %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)

l_EI_wMETA_gtex_a2g_filtered = l_EI_wMETA_gtex_a2g %>%
  filter(Sample %in% allData_gtex_filtered$Sample)

allData_gtex_filtered_mergedReps = allData_gtex_filtered %>%
  group_by(Donor, Tissue, TissueHistologicalType, Sex, Cohort, CohortOriginal, Age, AgeGroup, BMI) %>%
  summarise(across(contains("EditingIndex"), mean),
            ReplicateNum = n()) 

l_EI_wMETA_gtex_a2g_filtered_mergedReps = l_EI_wMETA_gtex_a2g_filtered %>%
  group_by(EI_type, `Editing Index Type`, Donor, Tissue, TissueHistologicalType, Sex, Cohort, CohortOriginal, Age, AgeGroup, BMI) %>%
  summarise(EditingIndex = mean(EditingIndex),
            ReplicateNum = n()) 


# Linear regression -------------------------------------------------------

lm_per_cohort_not_through_zero = allData_gtex_filtered_mergedReps %>%
  # 1. group & bundle each tissue’s data
  group_by(Tissue, Cohort) %>%
  nest() %>%
  # 2. fit lm(...) for each tissue
  mutate(fit = map(data, ~ lm(A2GEditingIndexIRAlu3UTR ~ A2GEditingIndexGenomeWideAlu, 
                              data = .x)), 
         # 3. tidy() pulls out estimate & std.error
         tidy   = map(fit, broom::tidy),
         n = map_int(data, nrow)) %>%
  # 4. unnest & keep only the slope row
  unnest(tidy)

lm_slopes_per_cohort_not_through_zero = lm_per_cohort_not_through_zero %>%
  filter(term == "A2GEditingIndexGenomeWideAlu") %>%
  # 5. select/rename into a clean tibble
  select(Tissue, Cohort, slope = estimate, slope_SD = std.error, p_value = p.value, n)

# Fig. 5F -------------------------------------------------------------

estimates_vec_antemortem = lm_per_cohort_not_through_zero %>%
  filter(Cohort == "Antemortem",
         term == "A2GEditingIndexGenomeWideAlu") %>%
  { setNames(pull(., estimate), pull(., Tissue)) }


# antemortem slopes by order
fig5f = lm_slopes_per_cohort_not_through_zero %>%
  filter(Cohort == "Antemortem", n >= 10) %>%
  ggplot(aes(x = forcats::fct_reorder(Tissue, -slope), y = slope, color = Tissue)) +
  geom_point(size = 3) +
  # facet_grid(.~Cohort, scales = "free", space = "free") +
  geom_errorbar(aes(ymin = slope-slope_SD,ymax = slope+slope_SD), width = 0.5) + 
  theme_custom(legend_position = "none")+
  scale_color_manual(values = gtex_colors) +
  ggtitle("CEI/AEI ratio for different GTEx tissues (antemortem)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  expand_limits(y = 0) +
  scale_y_continuous(breaks = seq(0, 2, 0.25)) +
  xlab("Tissue") +
  ylab("Slope") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme(axis.title.x = element_blank()) 

fig5f_inset = allData_gtex_filtered_mergedReps %>%
  filter(Tissue == "Muscle - Skeletal" | Tissue == "Stomach", Cohort == "Antemortem") %>%
  ggplot(aes(x = A2GEditingIndexGenomeWideAlu, y = A2GEditingIndexIRAlu3UTR, color = Tissue)) +
  geom_point(size = 0.5) +
  theme_custom(legend_position = "right")+
  scale_color_manual(values = gtex_colors, 
                     labels = function(l) { paste0(l, " (slope = ", signif(estimates_vec_antemortem[l], digits = 3), ")")}) +
  # ggtitle("Tissues with highest and lowest slopes") +
  geom_abline(data = lm_per_cohort_not_through_zero %>%
                filter(Tissue == "Muscle - Skeletal" | Tissue == "Stomach", Cohort == "Antemortem") %>%
                pivot_wider(names_from = term, values_from = estimate, id_cols = c(Tissue, Cohort)),
              aes(slope = A2GEditingIndexGenomeWideAlu,intercept = `(Intercept)`, colour = Tissue)) +
  coord_cartesian(xlim = c(0, 2.2), ylim = c(0, 4.6), expand = FALSE, default = FALSE, clip = "on") +
  xlab("AEI") +
  ylab("CEI") +
  guides(color = guide_legend(title = "Tissues with top\nand bottom slopes")) +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA))+
  # annotate("text",x=1,y=0.5, size = 5, hjust = 0, label=(paste0("Slope==",signif(lm_per_cohort_not_through_zero %>%
  #                                                                                    filter(Tissue == "Muscle - Skeletal", 
  #                                                                                           Cohort == "Antemortem", 
  #                                                                                           term == "A2GEditingIndexGenomeWideAlu") %>% 
  #                                                                                    pull(estimate), digits = 3))),
  #          parse=TRUE, color = gtex_colors["Muscle - Skeletal"]) +
  # annotate("text",x=0.1,y=4, size = 5, hjust = 0, label=(paste0("Slope==",signif(lm_per_cohort_not_through_zero %>%
  #                                                                                  filter(Tissue == "Stomach", 
  #                                                                                         Cohort == "Antemortem", 
  #                                                                                         term == "A2GEditingIndexGenomeWideAlu") %>% 
  #                                                                                  pull(estimate), digits = 3))),
  #          parse=TRUE, color = gtex_colors["Stomach"]) +
  scale_x_continuous(breaks = seq(0, 2, 1))

# Fig. 5G -----------------------------------------------------------------
tissues_w2cohorts = allData_gtex_filtered_mergedReps %>%
  group_by(Tissue, Cohort) %>%
  tally %>%
  filter(n >=10) %>%
  group_by(Tissue) %>%
  tally%>%
  filter(n > 1) %>%
  pull(Tissue)

# compute std for ratios
curr_data_fig5g = l_EI_wMETA_gtex_a2g %>%
  filter(Tissue%in%tissues_w2cohorts)%>%
  group_by(Tissue, Cohort, `Editing Index Type`) %>%
  summarise(across(EditingIndex, list(mean = mean, sd = sd)),
            n = n()) %>%
  ungroup()

ratio_table_fig5g = curr_data_fig5g %>%
  pivot_wider(names_from = Cohort, values_from = EditingIndex_mean, id_cols = c(Tissue, `Editing Index Type`)) %>%
  mutate("Ratio" = Postmortem / Antemortem)

# Compute sd of ratios
stats_fig5g <- curr_data_fig5g %>%
  # compute normalized sd and normalized sd squared
  mutate(norm_sd_editing = EditingIndex_sd / EditingIndex_mean,
         norm_sd_editing_squared = norm_sd_editing^2,
         # divide by n-1
         final = norm_sd_editing_squared / n) %>%
  # group by same, NOT by treatment
  group_by(Tissue, `Editing Index Type`) %>%
  # compute the square root of the sum of normalized sd squared
  summarise(sqrt_sum_of_norm_sd_editing_squared = sqrt(sum(final))) %>%
  # add overall mean per group
  inner_join(ratio_table_fig5g %>%
               ungroup() %>%
               select(Tissue, `Editing Index Type`, Ratio)) %>%
  mutate(final_sd = sqrt_sum_of_norm_sd_editing_squared * Ratio)

tissue_order = ratio_table_fig5g %>%
  filter(`Editing Index Type`=="Global Index")%>%
  arrange(-Ratio) %>%
  pull(Tissue)

fig5g = ratio_table_fig5g %>%
  inner_join(stats_fig5g) %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order)) %>%
  mutate(`Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "Global Index", 
                                                    "CEI" = "Cytoplasmic Index") %>% 
           forcats::fct_relevel("CEI")) %>%
  ggplot(aes(x = Tissue, y = Ratio, fill = Tissue)) +
  geom_col(position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = Ratio-final_sd, ymax = Ratio+final_sd), position = position_dodge(width = 0.9), width = 0.2) +
  facet_grid(.~`Editing Index Type`) +
  theme_custom(legend_position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = gtex_colors) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Tissue") +
  ylab("Mean postmortem index / mean antemortem index") +
  ggtitle("Postmortem to antemortem mean A2G editing index ratio") +
  theme(axis.title.x = element_blank()) 

# inset
chosen_tissue = "Adipose - Visceral (Omentum)"

means_df = l_EI_wMETA_gtex_a2g_filtered_mergedReps %>% 
  filter(Tissue == chosen_tissue) %>%
  group_by(`Editing Index Type`, Cohort) %>%
  summarise(mean_editing = mean(EditingIndex))

means_vec_cyto = means_df %>%
  filter(`Editing Index Type` == "Cytoplasmic Index") %>%
  { setNames(pull(., mean_editing), pull(., Cohort)) }

fig5g_inset_cytoplsamic = l_EI_wMETA_gtex_a2g_filtered_mergedReps %>% 
  filter(Tissue == chosen_tissue, `Editing Index Type` == "Cytoplasmic Index") %>%
  ggplot(aes(x = EditingIndex, fill = Cohort)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  # ggtitle(paste("Editing distribution"),
  #         subtitle = "Binwidth = 0.1") +
  # facet_grid(.~`Editing Index Type`) +
  # scale_x_log10() +
  theme_custom(legend_position = "left") +
  scale_fill_manual(values = case_control_colors,
                    labels = function(l) { paste0(l, "\n(mean = ", signif(means_vec_cyto[l], digits = 3), ")")}) +
  xlab("CEI") +
  ylab("# Samples") +
  expand_limits(x = 0, y = 38) +
  geom_vline(data = means_df %>%
               filter(`Editing Index Type` == "Cytoplasmic Index"), aes(xintercept = mean_editing, color = Cohort), 
             linetype = "dashed") +
  scale_color_manual(values = colorspace::darken(case_control_colors, amount = 0.4)) +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA),
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.box.margin   = margin(t = 0, r = -5, b = 0, l = 0, unit = "pt"))+
  scale_y_continuous(expand = c(0.005, 1)) +
  guides(fill = guide_legend(title = chosen_tissue %>% 
                               str_remove(pattern = "\\(.*\\)") %>%
                               str_replace(pattern = " - ", replacement = " ")), 
         color = "none")

means_vec_global = means_df %>%
  filter(`Editing Index Type` == "Global Index") %>%
  { setNames(pull(., mean_editing), pull(., Cohort)) }

fig5g_inset_global = l_EI_wMETA_gtex_a2g_filtered_mergedReps %>% 
  filter(Tissue == chosen_tissue, `Editing Index Type` == "Global Index") %>%
  ggplot(aes(x = EditingIndex, fill = Cohort)) +
  geom_histogram(binwidth = 0.1, alpha = 0.6, position = "identity") +
  # ggtitle(paste("Editing distribution"),
  #         subtitle = "Binwidth = 0.1") +
  # facet_grid(.~`Editing Index Type`) +
  # scale_x_log10() +
  theme_custom(legend_position = "left") +
  scale_fill_manual(values = case_control_colors,
                    labels = function(l) { paste0(l, "\n(mean = ", signif(means_vec_global[l], digits = 3), ")")}) +
  xlab("AEI") +
  ylab("# Samples") +
  expand_limits(x = c(0, 5.641938)) +
  geom_vline(data = means_df %>%
               filter( `Editing Index Type` == "Global Index"), aes(xintercept = mean_editing, color = Cohort), 
             linetype = "dashed") +
  scale_color_manual(values = colorspace::darken(case_control_colors, amount = 0.4)) +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA),
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.box.margin   = margin(t = 0, r = -5, b = 0, l = 0, unit = "pt"))+
  scale_y_continuous(expand = c(0.005, 1)) +
  guides(fill = guide_legend(title = chosen_tissue %>% 
                               str_remove(pattern = "\\(.*\\)") %>%
                               str_replace(pattern = " - ", replacement = " ")), 
         color = "none")
    
  



## ***********************************************************************************************


# COMBINE -----------------------------------------------------------------


library(cowplot)
fig5 <- plot_grid(# infections
                  plot_grid(fig5a, fig5b, fig5c, rel_widths = c(1, 2, 2), 
                            labels=c("A", "B", "C"), ncol = 3, nrow = 1, align = 'hv', axis = "tlb", label_size=18),
                  # TCGA
                  plot_grid(fig5d, fig5e, rel_widths = c(1, 1.2),
                            labels=c("D", "E"), ncol = 2, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  # GTEx
                  plot_grid(fig5f +
                              # inset
                              patchwork::inset_element(fig5f_inset,
                                                       left   = 0.01, right = 0.7,
                                                       bottom = 0.01, top   = 0.5), 
                            fig5g +
                              # inset
                              patchwork::inset_element(fig5g_inset_global,
                                                       left   = 0.72, right = 1,
                                                       bottom = 0.66, top   = 1) +
                              # inset
                              patchwork::inset_element(fig5g_inset_cytoplsamic,
                                                       left   = 0.215, right = 0.495,
                                                       bottom = 0.66, top   = 1),
                            rel_widths = c(1, 1.8),
                            labels=c("F", "G"), ncol = 2, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  # rel_heights = c(1, 0.8, 1, 1),
                  labels=c("", "", "", ""), ncol = 1, nrow = 3, align = 'hv', axis = "ltb", label_size=18)

save_plot(file.path(out_plots,"Figure5.pdf"), fig5, ncol = 1, nrow = 3, base_height = 9, base_width = 20)
save_plot(file.path(out_plots,"Figure5.png"), fig5, ncol = 1, nrow = 3, base_height = 9, base_width = 20)




