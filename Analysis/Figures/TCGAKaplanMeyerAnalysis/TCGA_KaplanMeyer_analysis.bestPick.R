rm(list = ls(all = TRUE))
# SETUP ===========================================================================================
# *************************************************************************************************

# RNA Editing Dimensionality Reduction & Covariate Analysis (TCGA, Cytoplasmic Regions)
# Input: editing_matrix.csv (rows = samples, columns = editing regions)
#        metadata.csv (columns = sample, tissue, condition)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
# library(umap)
# library(Rtsne)

### Output directory --------------------------------------------------------
outdir = "Analysis/"
out_plots="plots/TCGA"
out_boxplots = "boxplots"
out_distributions = "distributions"
out_scatterplots = "scatterplots"
out_correlations = "correlations"
out_barplots = "barplots"
out_survival = "survival"
out_stats="stats/TCGA"
out_tables = "tables_TCGA"


### create directories ------------------------------------------------------
#create paths, change
out_plots = file.path(outdir, out_plots)
out_boxplots = file.path(out_plots, out_boxplots)
out_distributions = file.path(out_plots, out_distributions)
out_scatterplots = file.path(out_plots, out_scatterplots)
out_correlations = file.path(out_plots, out_correlations)
out_barplots = file.path(out_plots, out_barplots)
out_survival = file.path(out_plots, out_survival)
out_stats = file.path(outdir, out_stats)
out_tables = file.path(outdir, out_tables)
#create directories, if do not exist
dir.create(out_boxplots, recursive = T)
dir.create(out_distributions, recursive = T)
dir.create(out_scatterplots, recursive = T)
dir.create(out_correlations, recursive = T)
dir.create(out_barplots, recursive = T)
dir.create(out_survival, recursive = T)
dir.create(out_stats, recursive = T)
dir.create(out_tables, recursive = T)


# Load Data ---------------------------------------

#metadata
## SRA table ---------------------------------------
metadata = "TCGA/Metadata/TCGA_all_samples_genome_sample_sheet.tsv" %>% 
  fread(stringsAsFactors = F) %>% 
  select("File ID", "Project ID", "Case ID", "Sample ID", "Sample Type")  %>%
  mutate(Tissue = str_remove(`Project ID`, pattern = "TCGA-")) %>%
  rename("Sample" = `File ID`)

## survival ---------------
clinical_data = fread("TCGA/Metadata/combined_study_clinical_data.tsv", stringsAsFactors = F) %>%
  select(`Patient ID`, `Diagnosis Age`, 
         `Overall Survival (Months)`, `Overall Survival Status`,
         Sex, `Prior Treatment`, `Race Category`, 
         `Patient's Vital Status`, `Year of Death`, `Year of Diagnosis`) %>%
  distinct() 

## data ---------------
cytoindex_data = fread("TCGA/Summary/TCGA_3UTREditingIndex.csv", stringsAsFactors = F)
global_index_data = fread("TCGA/Metadata/TCGA_filtered_AluEditingIndex.csv", stringsAsFactors = F)

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
  ungroup() %>%
  left_join(clinical_data, by = join_by("Case ID" == "Patient ID"))


l_EI_wMETA_a2g_filtered_mergedReps = allData_filtered_mergedReps %>%
  select(`Project ID`:Tissue, contains("A2GEditingIndex")) %>%
  pivot_longer(cols = contains("A2GEditingIndex"), names_to = "Editing Index Type", 
               values_to = "EditingIndex", names_prefix = "A2GEditingIndex") %>%
  left_join(clinical_data, by = join_by("Case ID" == "Patient ID"))

matched_samples = allData_filtered_mergedReps  %>% 
  filter(`Sample Type`=="Primary Tumor" | `Sample Type`=="Solid Tissue Normal")%>% 
  group_by(`Case ID`, Tissue, `Sample Type`) %>% 
  tally%>% 
  pivot_wider(names_from = `Sample Type`, values_from = n) %>%
  filter(!is.na(`Primary Tumor`), !is.na(`Solid Tissue Normal`)) %>%
  pull(`Case ID`)


# Survival Data -----------------------------------------------------------

allData_filtered_mergedReps_survival = allData_filtered_mergedReps %>%
  filter(`Overall Survival (Months)` > 0, `Sample Type`=="Primary Tumor") %>%
  # Step 5: Create high/low group and survival object 
  mutate(group_cytoplasmic = ifelse(A2GEditingIndexIRAlu3UTR >= median(A2GEditingIndexIRAlu3UTR, na.rm = TRUE), "High", "Low"),
         group_global = ifelse(A2GEditingIndexGenomeWideAlu >= median(A2GEditingIndexGenomeWideAlu, na.rm = TRUE), "High", "Low")) %>%
  separate(`Overall Survival Status`, sep = ":", into = c("status", "status_description"), convert = T)




# Survival analysis -------------------------------------------------------
# ALL TISSUES Q1 vs. Q4 

library(survival)
# library(survminer)
library(ggsurvfit)

# Prepare data
base_data <- allData_filtered_mergedReps %>%
  filter(`Overall Survival (Months)` > 0, `Sample Type` == "Primary Tumor") %>%
  separate(`Overall Survival Status`, sep = ":", into = c("status", "status_description"), convert = TRUE)


# Define a function for survival analysis
perform_survival_analysis <- function(tissue_name, variable_name, base_data, quantile_threshold = 0.5) {
  print(tissue_name)
  
  # Filter data for the specific tissue and prepare survival variables
  df_tissue <- base_data %>%
    filter(Tissue == tissue_name) %>%
    rename(
      survival_time = `Overall Survival (Months)`,
      survival_status = status
    )
  
  # Calculate the quantile threshold for the variable
  threshold <- quantile(df_tissue[[variable_name]], probs = quantile_threshold, na.rm = TRUE)
  
  # Divide into Low/High groups based on the threshold
  df_tissue <- df_tissue %>%
    mutate(
      group_index = if_else(!!rlang::sym(variable_name) >= threshold, 
                            "High", "Low"),
      group_index = factor(group_index, levels = c("Low", "High"), ordered = TRUE)
    )
  
  # Skip analysis if not enough data
  if (nrow(df_tissue) < 10 || length(unique(df_tissue$group_index)) < 2) {
    return(tibble(
      tissue = tissue_name,
      variable = as.character(variable_name),
      quantile = quantile_threshold,
      n_low = df_tissue%>%filter(group_index == "Low") %>%nrow,
      n_high = df_tissue%>%filter(group_index == "High") %>%nrow,
      n_total = nrow(df_tissue),
      pval = NA_real_
    ))
  }
  
  # Perform survival analysis
  fit <- survfit(Surv(survival_time, survival_status) ~ group_index, data = df_tissue)
  
  pval <- tryCatch({
    survdiff_res <- survdiff(Surv(survival_time, survival_status) ~ group_index, data = df_tissue)
    1 - pchisq(survdiff_res$chisq, df = length(survdiff_res$n) - 1)
  }, error = function(e) NA_real_)
  
  # Return results
  tibble(
    tissue = tissue_name,
    variable = as.character(variable_name),
    quantile = quantile_threshold,
    n_low = df_tissue%>%filter(group_index == "Low") %>%nrow,
    n_high = df_tissue%>%filter(group_index == "High") %>%nrow,
    n_total = nrow(df_tissue),
    pval = signif(pval, 4)
  )
}

# Create combinations of tissues, variables, and quantiles
# Use purrr to run the analysis for all combinations
analysis_combinations <- expand_grid(
  tissue_name = unique(allData_filtered_mergedReps_survival$Tissue),
  variable_name = c("A2GEditingIndexGenomeWideAlu", "A2GEditingIndexIRAlu3UTR"),
  quantile_threshold = seq(0.1, 0.9, 0.1)
)
 
# Run the analyses and collect results
results <- pmap_df(analysis_combinations, perform_survival_analysis, base_data = base_data) %>%
  group_by(variable) %>%
  mutate(adjusted_pval = p.adjust(pval, "fdr"))

## Adjust P-value & write table -----
fwrite(x = results, file.path(out_tables, "KaplanMeier_Combinations_SurvivalAnalysis.PrimaryTumor.PositiveSurvival.csv"), quote = F, row.names = F, scipen = 999)


## pick best -----
best_results = results %>%
  group_by(tissue, variable) %>%
  slice_min(order_by = pval, n = 1)

# best pick
fwrite(x = best_results, file.path(out_tables, "KaplanMeier_Combinations_SurvivalAnalysis.PrimaryTumor.PositiveSurvival.BestResults.csv"), quote = F, row.names = F, scipen = 999)



# ************************************************************