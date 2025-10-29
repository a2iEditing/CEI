
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
comparisons = fread("Analysis/tables/Classifications.Metadata.comparisons.new.v2.nonewlines.csv", stringsAsFactors = F)

pathogen_disease_classification = fread("Analysis/tables/pathogen_disease_classification_withBioProjects.csv", stringsAsFactors = F) %>%
  separate_longer_delim(cols = c(BioProject, Case), delim = ";")

l_EI_wMETA = fread(file.path(input_dir, "EditingIndex.longFormat.wMetadata.csv"), stringsAsFactors = F) %>%
  mutate(disease_or_control = case_when(Disease_Status_manual_annotation == "Control" ~ "Control",
                                        Disease_Status_manual_annotation == "OtherDesign" ~ "OtherDesign",
                                        .default = "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case", "OtherDesign")),
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "Global Index" = "GenomeWideAlu", 
                                                    "Cytoplasmic Index" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("Cytoplasmic Index"))

allData = fread(file.path(input_dir, "AllData.wMetadata.csv"), stringsAsFactors = F) %>%
  mutate(disease_or_control = case_when(Disease_Status_manual_annotation == "Control" ~ "Control",
                                        Disease_Status_manual_annotation == "OtherDesign" ~ "OtherDesign",
                                        .default = "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case", "OtherDesign"))) %>%
  # add next best mismatch
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
         "Uniquely mapped Mbp" = (`Uniquely mapped reads number`/10^6)*`Average mapped length`)

l_EI_wMETA_a2g =  l_EI_wMETA %>%
  filter(Mismatch == "A2G") %>%
  inner_join(l_EI_wMETA %>%
               filter(Mismatch != "A2G") %>%
               group_by(Sample, EI_type) %>%
               summarise(NextBestMM = max(EditingIndex, na.rm = T))) %>% 
  mutate(Signal2NoiseRatio = EditingIndex / NextBestMM)

# case-control
allData_case_control = allData %>%
  filter(StudyDesign=="Case-Control")

l_EI_wMETA_case_control = l_EI_wMETA %>%
  filter(StudyDesign=="Case-Control")

l_EI_wMETA_a2g_case_control = l_EI_wMETA_a2g %>%
  filter(StudyDesign=="Case-Control")


# Filter ------------------------------------------------------------------

allData_filtered = allData %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)

allData_case_control_filtered = allData_case_control %>%
  filter(NextBestMMAfterA2G_GenomeWideAlu < 0.3, NextBestMMAfterA2G_IRAlu3UTR < 0.3)

l_EI_wMETA_a2g_filtered = l_EI_wMETA_a2g %>%
  filter(Sample %in% allData_filtered$Sample)

l_EI_wMETA_a2g_case_control_filtered = l_EI_wMETA_a2g_case_control %>%
  filter(Sample %in% allData_case_control_filtered$Sample)


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")


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
                             ~ two_groups_comparison(all_data = allData_case_control_filtered,
                                                     comparison_line = comparisons[.x,] %>% as.vector,
                                                     variable = "A2GEditingIndexGenomeWideAlu")) %>%
  rbind(map_dfr(.x = c(1:nrow(comparisons)),
                ~ two_groups_comparison(all_data = allData_case_control_filtered,
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



# Get comparison samples --------------------------------------------------


get_comparison_samples = function(all_data, comparison_line, variable) {
  print(comparison_line[["BioProject"]])
  
  # get current data
  curr_data = all_data %>%
    filter(BioProject == comparison_line[["BioProject"]],
           Disease_Status_manual_annotation == comparison_line[["Case"]] | Disease_Status_manual_annotation == comparison_line[["Control"]])
  
  
  # if there are two read lengths - use recursion to do this for each
  if (curr_data$ReadLength %>% unique %>% length > 1) {
    readlengths = curr_data$ReadLength %>% unique
    d = map_dfr(.x = readlengths, 
                .f = ~ get_comparison_samples(all_data = curr_data %>% 
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
                .f = ~ get_comparison_samples(all_data = curr_data %>% 
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
      group_by(Disease_Status_manual_annotation, BioProject, technical_repeats_id, paired_sample_id) %>%
      summarise(across(all_of(c(variable)), mean)) %>%
      # release grouping
      ungroup() %>%
      mutate("Sample" = technical_repeats_id)
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
      filter(if_all(everything(), ~!is.na(.x))) %>%
      pivot_longer(cols = c(comparison_line[["Case"]], comparison_line[["Control"]]), names_to = "Disease_Status_manual_annotation", values_to = variable) %>%
      rename("Sample" = "paired_sample_id") %>%
      mutate(BioProject = comparison_line[["BioProject"]])
    
  } 
  
  return(curr_data %>%
           pivot_longer(cols = variable, names_to = "EI_type", values_to = "EditingIndex") %>%
           select(Sample, Disease_Status_manual_annotation, BioProject, EditingIndex) %>%
           mutate(Variable = variable))
}

# get data points
comparison_samples_data = map_dfr(.x = c(1:nrow(comparisons)),
                                  ~ get_comparison_samples(all_data = allData_case_control_filtered,
                                                           comparison_line = comparisons[.x,] %>% as.vector,
                                                           variable = "A2GEditingIndexGenomeWideAlu")) %>%
  rbind(map_dfr(.x = c(1:nrow(comparisons)),
                ~ get_comparison_samples(all_data = allData_case_control_filtered,
                                         comparison_line = comparisons[.x,] %>% as.vector,
                                         variable = "A2GEditingIndexIRAlu3UTR")) ) %>%
  mutate(Tissue = replace_na(Tissue, ""),
         ReadLength = replace_na(as.character(ReadLength), ""),
         Data = paste(BioProject, Tissue, ReadLength)) %>%
  mutate(disease_or_control = if_else(Disease_Status_manual_annotation == "Control", "Control", "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case")),
         `Editing Index Type` = forcats::fct_recode(Variable, "AEI" = "A2GEditingIndexGenomeWideAlu", 
                                                    "CEI" = "A2GEditingIndexIRAlu3UTR") %>% 
           forcats::fct_relevel("CEI"))


# Supp. Fig. 4A -----------------------------------------------------------

data_order = c("PRJNA735648\n", 
               "PRJNA949617\nhiPSC-neurons", 
               "PRJNA949617\nhiPSC-macrophages", 
               "PRJNA789591\n", 
               "PRJNA638819\n", 
               "PRJNA660611\n", 
               "PRJNA600939\n", 
               "PRJNA285798\n", 
               "PRJNA551288\n", 
               "PRJNA882083\n", 
               "PRJNA285953\n", 
               "PRJNA227074\n", 
               "PRJNA723500\n", 
               "PRJNA918345\n", 
               "PRJNA692462\n", 
               "PRJNA352062\n", 
               "PRJNA695511\n", 
               "PRJNA750782\n", 
               "PRJNA562638\nHuman colon cancer cell line (HT-29)", 
               "PRJNA485140\n")
suppFig4a = comparison_samples_data %>%
  #filter significant comparisons on cytoindex only
  inner_join(comparison_results %>%
               mutate(Significant = PValue < 0.05) %>%
               select(BioProject, TestedCase, TestedControl, Variable, Significant, Tissue, ReadLength) %>%
               pivot_wider(values_from = Significant, names_from = Variable, names_prefix = "Significant_") %>%
               filter(!Significant_A2GEditingIndexGenomeWideAlu, Significant_A2GEditingIndexIRAlu3UTR) %>%
               select(BioProject, TestedControl, TestedCase, Tissue, ReadLength) %>%
               pivot_longer(cols = c(TestedControl, TestedCase), names_to = "X", values_to = "Disease_Status_manual_annotation") %>%
               select(-X) %>%
               distinct %>%
               mutate(Tissue = replace_na(Tissue, ""),
                      ReadLength = replace_na(as.character(ReadLength), ""))) %>%
  # make control always first
  mutate(Disease_Status_manual_annotation = forcats::fct_relevel(Disease_Status_manual_annotation, "Control"),
         # rank = rank(Disease_Status_manual_annotation, ties.method = "min"),
         Data = paste0(BioProject, "\n", replace_na(Tissue, "")),
         Data = factor(Data, levels = data_order, ordered = T)) %>% 
  ggplot(aes(x = Disease_Status_manual_annotation, y = EditingIndex)) +
  geom_boxplot(aes( fill = disease_or_control), outliers = F) +
  # geom_jitter(width = 0.05, alpha = 0.5) +
  theme_custom(legend_position = "none") +
  ggtitle("Infectious diseases comparisons that are significant only in CEI") +
  # Cohen's D = (Mean editing index case - mean editing index control) / editing index sd
  scale_fill_manual(values = chosen_colors) +
  facet_grid(`Editing Index Type`~Data, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        strip.text.x.top = element_text(size = 7)) +
  scale_y_continuous(breaks = c(0:7)) +
  ggpubr::stat_pvalue_manual(
    comparison_results %>%
      inner_join(comparison_results %>%
                   mutate(Significant = PValue < 0.05) %>%
                   select(BioProject, TestedCase, TestedControl, Variable, Significant, Tissue, ReadLength) %>%
                   pivot_wider(values_from = Significant, names_from = Variable, names_prefix = "Significant_") %>%
                   filter(!Significant_A2GEditingIndexGenomeWideAlu, Significant_A2GEditingIndexIRAlu3UTR) %>%
                   select(BioProject, TestedControl, TestedCase, Tissue, ReadLength)) %>%
      mutate(#Data = paste(BioProject, replace_na(Tissue, ""), replace_na(as.character(ReadLength), "")),
             Data = paste0(BioProject, "\n", replace_na(Tissue, "")),
             Data = factor(Data, levels = data_order, ordered = T),
             `Editing Index Type` = forcats::fct_recode(Variable, "AEI" = "A2GEditingIndexGenomeWideAlu", 
                                                        "CEI" = "A2GEditingIndexIRAlu3UTR"), 
             group1 = TestedControl,
             group2 = TestedCase,
             p = PValue, 
             y.position = if_else(Variable == "A2GEditingIndexIRAlu3UTR", 4.7, 1.8),
             # fix for multiple comparisons
             y.position = case_when(group2 == "Salmonella" | group2 == "Viral" ~ y.position - 0.25,
                                    group2 == "Shigella" ~ y.position - 0.5,
                                    .default = y.position),
             label = signif(x = p, digits = 2)) %>%
      select(Data, `Editing Index Type`, group1, group2, p, y.position, label),
    tip.length = 0.01,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position"
  ) +
  ylab("Editing Index")  +
  expand_limits(y=0) +
  expand_limits(y=1.9) +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank()) 



# Supp. Fig. 4B -----------------------------------------------------------

suppFig4b = comparison_results %>%
  mutate(Significant = PValue < 0.05) %>%
  select(BioProject, Comparison, TestedCase, TestedControl, CohensDEffectSize, Variable, Significant) %>%
  pivot_wider(values_from = c(CohensDEffectSize, Significant), names_from = Variable) %>%
  mutate(Significant = case_when(Significant_A2GEditingIndexGenomeWideAlu & Significant_A2GEditingIndexIRAlu3UTR ~ "Both Comparisons", 
                                 Significant_A2GEditingIndexGenomeWideAlu ~ "AEI Comparison",
                                 Significant_A2GEditingIndexIRAlu3UTR ~ "CEI Comparison",
                                 .default = "None") %>%
           factor(levels = c("Both Comparisons", "CEI Comparison", "AEI Comparison", "None"))) %>%
  left_join(pathogen_disease_classification, by = join_by("TestedCase"=="Case", "BioProject")) %>%
  mutate(`Pathogen Type` = replace_na(`Pathogen Type`, "Unknown pathogen")) %>%
  ggplot(aes(x = CohensDEffectSize_A2GEditingIndexGenomeWideAlu, y = CohensDEffectSize_A2GEditingIndexIRAlu3UTR, color = Significant)) +
  # draw lines
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept =  0 , linetype = "dashed", color = "grey") +
  geom_vline(xintercept =  0 , linetype = "dashed", color = "grey") +
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
  geom_point(size = 2) +
  theme_custom() +
  facet_grid(.~`Pathogen Type`)+
  ggtitle("Effect size of case-control comparison per editing index type across datasets and pathogen types",
          subtitle = "Green triangle marks positive (increasing) effect size greater in CEI")  +
  # Cohen's D = (Mean editing index case - mean editing index control) / editing index sd
  scale_color_manual(values = c("#8B1A1A", "#e37e00", "#1a476f", "grey60")) +
  expand_limits(y = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize)),
                x = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize))) +
  xlab("AEI effect size (Cohen's D)") +
  ylab("CEI effect size (Cohen's D)") +
  coord_cartesian(xlim = c(-2, 2.01), ylim = c(-2, 2.01), expand = FALSE, default = FALSE, clip = "on")


# Supp. Fig. 4C -----------------------------------------------------------

suppFig4c = comparison_results %>%
  mutate(Significant = PValue < 0.05) %>%
  select(BioProject, Comparison, TestedCase, TestedControl, CohensDEffectSize, Variable, Significant) %>%
  pivot_wider(values_from = c(CohensDEffectSize, Significant), names_from = Variable) %>%
  mutate(Significant = case_when(Significant_A2GEditingIndexGenomeWideAlu & Significant_A2GEditingIndexIRAlu3UTR ~ "Both Comparisons", 
                                 Significant_A2GEditingIndexGenomeWideAlu ~ "AEI Comparison",
                                 Significant_A2GEditingIndexIRAlu3UTR ~ "CEI Comparison",
                                 .default = "None") %>%
           factor(levels = c("Both Comparisons", "CEI Comparison", "AEI Comparison", "None"))) %>%
  left_join(pathogen_disease_classification, by = join_by("TestedCase"=="Case", "BioProject")) %>%
  mutate(`Clinical Manifestation` = replace_na(`Clinical Manifestation`, "Other"),
         `Clinical Manifestation` = forcats::fct_other(`Clinical Manifestation`, 
                                                       keep = c("Dermatological", "Febrile", "Gastrointestinal", "Respiratory"), 
                                                       other_level = "Other/Unknown")) %>%
  ggplot(aes(x = CohensDEffectSize_A2GEditingIndexGenomeWideAlu, y = CohensDEffectSize_A2GEditingIndexIRAlu3UTR, color = Significant)) +
  # draw lines
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept =  0 , linetype = "dashed", color = "grey") +
  geom_vline(xintercept =  0 , linetype = "dashed", color = "grey") +
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
  geom_point(size = 2) +
  theme_custom() +
  facet_grid(.~`Clinical Manifestation`)+
  ggtitle("Effect size of case-control comparison per editing index type across datasets and clinical manifestations",
          subtitle = "Green triangle marks positive (increasing) effect size greater in CEI")  +
  # Cohen's D = (Mean editing index case - mean editing index control) / editing index sd
  scale_color_manual(values = c("#8B1A1A", "#e37e00", "#1a476f", "grey60")) +
  expand_limits(y = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize)),
                x = c(max(comparison_results$CohensDEffectSize), min(comparison_results$CohensDEffectSize))) +
  xlab("AEI effect size (Cohen's D)") +
  ylab("CEI effect size (Cohen's D)") +
  coord_cartesian(xlim = c(-2, 2.01), ylim = c(-2, 2.01), expand = FALSE, default = FALSE, clip = "on")

# Combine -----------------------------------------------------------------


library(cowplot)
suppFig4 <- plot_grid(suppFig4a, suppFig4b, suppFig4c, rel_heights = c(1.5, 1, 1),
                      labels=c("AUTO"), ncol = 1, nrow = 3, align = 'hv', axis = "ltb", label_size=18)
save_plot(file.path(out_plots,"SuppFig5.1_Infections.pdf"), suppFig4, ncol = 1, nrow = 3, base_height = 7, base_width = 18)
save_plot(file.path(out_plots,"SuppFig5.1_Infections.png"), suppFig4, ncol = 1, nrow = 3, base_height = 7, base_width = 18)
