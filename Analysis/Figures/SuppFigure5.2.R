
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

### create directories ------------------------------------------------------
#create paths, change
out_plots = file.path(outdir, out_plots)
out_stats = file.path(outdir, out_stats)
#create directories, if do not exist
dir.create(out_plots, recursive = T)
dir.create(out_stats, recursive = T)


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
                                           "GenomeWideAlu" ~ "AEI",
                                           "IRAlu3UTR" ~ "CEI") %>%
           forcats::fct_relevel("CEI")) %>%
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


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")


# Supp. Fig5.2A ---------------------------------------------------------------
suppFig5.2a = ggplot(l_EI_wMETA_gtex_a2g_filtered_mergedReps, aes(x = Tissue, y = EditingIndex, fill = Tissue)) +
  # geom_violin(position = position_dodge(0.8)) +
  geom_boxplot() +
  facet_wrap(~ `Editing Index Type`) +
  scale_fill_manual(values = gtex_colors) +
  theme_custom(legend_position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("CEI and AEI across the GTEx database")+
  ylab("Editing Index") +
  theme(axis.title.x = element_blank())


# Supp. Fig5.2B ---------------------------------------------------------------
suppFig5.2b = l_EI_wMETA_gtex_a2g_filtered_mergedReps %>%
  group_by(Tissue) %>%
  filter(n_distinct(Cohort) == 2) %>%
  ggplot(aes(x = Tissue, y = EditingIndex, fill = Cohort)) +
  # geom_violin(position = position_dodge(0.8)) +
  geom_boxplot() +
  facet_wrap(~ `Editing Index Type`) +
  scale_fill_manual(values = chosen_colors) +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
  ylab("Editing Index")  +
  ggtitle("Differences in CEI and AEI across cohorts")+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

# Supp. Fig5.2C ---------------------------------------------------------------
tissues_w2cohorts_atLeast10samples = allData_gtex_filtered_mergedReps %>%
  group_by(Tissue, Cohort) %>%
  tally %>%
  filter(n >=10) %>%
  group_by(Tissue) %>%
  tally%>%
  filter(n > 1) %>%
  pull(Tissue)

# tissue order by diff
levels = allData_gtex_filtered_mergedReps %>%
  filter(Tissue%in%tissues_w2cohorts_atLeast10samples)%>%
  mutate("CEI/AEI Ratio" = A2GEditingIndexIRAlu3UTR / A2GEditingIndexGenomeWideAlu) %>%
  group_by(Tissue, Cohort) %>%
  summarise(mean = mean(`CEI/AEI Ratio`)) %>%
  pivot_wider(names_from = Cohort, values_from = mean, id_cols = Tissue) %>%
  mutate(Diff = Antemortem - Postmortem) %>%
  arrange(-Diff) %>%
  pull(Tissue)


# ratio between cohorts
suppFig5.2c = allData_gtex_filtered_mergedReps %>%
  filter(Tissue%in%tissues_w2cohorts_atLeast10samples)%>%
  mutate("CEI / AEI ratio" = A2GEditingIndexIRAlu3UTR / A2GEditingIndexGenomeWideAlu,
         Tissue = factor(Tissue, levels = levels)) %>%
  ggplot(aes(x = Tissue, y = `CEI / AEI ratio`, fill = Cohort)) +
  geom_boxplot() +
  theme_custom() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_manual(values = case_control_colors) +
  ggtitle("Ratio of CEI to AEI per sample")+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

# Supp. Fig5.2D ---------------------------------------------------------------
### slopes -------
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

suppFig5.2d = lm_slopes_per_cohort_not_through_zero %>%
  filter(Tissue %in% tissues_w2cohorts_atLeast10samples) %>%
  pivot_wider(names_from = Cohort, values_from = c(slope, slope_SD, n), id_cols = c(Tissue)) %>%
  mutate(# horizontal bar covers y=x if |Postmortem - Antemortem| ≤ SD_Postmortem
    h_cross = abs(slope_Postmortem - slope_Antemortem) <= slope_SD_Postmortem,
    # vertical bar covers y=x if |Postmortem - Antemortem| ≤ SD_Antemortem
    v_cross = abs(slope_Postmortem - slope_Antemortem) <= slope_SD_Antemortem,
    crosses_identity = h_cross | v_cross,
    Tissue = str_remove(string = Tissue, pattern = " \\(.*\\)")) %>%
  ggplot(aes(x = slope_Postmortem, y = slope_Antemortem, color = Tissue, alpha = crosses_identity)) +
  # facet_grid(.~Cohort, scales = "free", space = "free") +
  geom_errorbarh(aes(xmin = slope_Postmortem-slope_SD_Postmortem,xmax = slope_Postmortem+slope_SD_Postmortem),
                 height = 0.01, alpha = 0.5) +
  geom_errorbar(aes(ymin = slope_Antemortem-slope_SD_Antemortem,ymax = slope_Antemortem+slope_SD_Antemortem),
                width = 0.01, alpha = 0.5) +
  geom_point(size = 3) +
  theme_custom()+
  # expand_limits(y = 0, x = 0) +
  scale_color_manual(values = gtex_colors, 
                     labels = function(l) { str_replace(string = l, pattern = "Terminal", replacement = "\n  Terminal")}) +
  scale_alpha_manual(values = c(`TRUE`  = 0.3,   # error‐bars cross identity → semi‐transparent
                                `FALSE` = 1.0),  # no crossing      → fully opaque
                     guide  = FALSE) +
  ggtitle("Slopes of antemortem and postmortem cohorts (CEI ~ AEI)") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  xlab("Antemortem slope") +
  ylab("Postmortem slope") +
  guides(color = guide_legend(ncol = 3))  +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 9),
        legend.key.spacing.y = unit(0.05, "cm"),
        legend.box.margin   = margin(t = 0, r = 0, b = 0, l = -50, unit = "pt")) 


# Supp. Fig.5.2E ----------------------------------------------------------
estimates_vec_postmortem = lm_per_cohort_not_through_zero %>%
  filter(Cohort == "Postmortem",
         term == "A2GEditingIndexGenomeWideAlu") %>%
  { setNames(pull(., estimate), pull(., Tissue)) }


# Postmortem slopes by order
suppFig5.2e = lm_slopes_per_cohort_not_through_zero %>%
  filter(Cohort == "Postmortem", n >= 10) %>%
  ggplot(aes(x = forcats::fct_reorder(Tissue, -slope), y = slope, color = Tissue)) +
  geom_point(size = 3) +
  # facet_grid(.~Cohort, scales = "free", space = "free") +
  geom_errorbar(aes(ymin = slope-slope_SD,ymax = slope+slope_SD), width = 0.5) + 
  theme_custom(legend_position = "none")+
  scale_color_manual(values = gtex_colors) +
  ggtitle("CEI/AEI ratio for different GTEx tissues (postmortem)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(breaks = seq(0, 2, 0.25)) +
  xlab("Tissue") +
  ylab("Slope")
  geom_hline(yintercept = 1, linetype = "dashed") 

suppFig5.2e_inset = allData_gtex_filtered_mergedReps %>%
  filter(Tissue == "Liver" | Tissue == "Brain - Cortex", Cohort == "Postmortem") %>%
  ggplot(aes(x = A2GEditingIndexGenomeWideAlu, y = A2GEditingIndexIRAlu3UTR, color = Tissue)) +
  geom_point(size = 0.5) +
  theme_custom(legend_position = "right")+
  scale_color_manual(values = gtex_colors, 
                     labels = function(l) { paste0(l, " (slope = ", signif(estimates_vec_postmortem[l], digits = 3), ")")}) +
  # ggtitle("Tissues with highest and lowest slopes") +
  geom_abline(data = lm_per_cohort_not_through_zero %>%
                filter(Tissue == "Liver" | Tissue == "Brain - Cortex", Cohort == "Postmortem") %>%
                pivot_wider(names_from = term, values_from = estimate, id_cols = c(Tissue, Cohort)),
              aes(slope = A2GEditingIndexGenomeWideAlu,intercept = `(Intercept)`, colour = Tissue)) +
  coord_cartesian(xlim = c(0, 2.8), ylim = c(0, 6.7), expand = FALSE, default = FALSE, clip = "on") +
  xlab("AEI") +
  ylab("CEI") +
  guides(color = guide_legend(title = "Tissues with top\nand bottom slopes")) +
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent",   colour = NA),
        legend.background = element_rect(fill = "transparent",   colour = NA))+
  scale_x_continuous(breaks = seq(0, 2, 1))

# Combine -----------------------------------------------------------------


library(cowplot)

suppFig5.2 = plot_grid(plot_grid(suppFig5.2a,
                                 labels=c("A"), ncol = 1, nrow = 1, align = 'hv', axis = "ltb", label_size=18),
                       plot_grid(suppFig5.2b,
                                 labels=c("B"), ncol = 1, nrow = 1, align = 'hv', axis = "ltb", label_size=18),
                       plot_grid(suppFig5.2c,
                                 suppFig5.2e +
                                   # inset
                                   patchwork::inset_element(suppFig5.2e_inset,
                                                            left   = 0.01, right = 0.7,
                                                            bottom = 0.01, top   = 0.45),
                                 suppFig5.2d,
                                 rel_widths = c(1.1, 1.2, 1),
                                 labels=c("C", "D", "E"), ncol = 3, nrow = 1, align = 'hv', axis = "ltb", label_size=18),
                       labels=c("", "", ""), ncol = 1, nrow = 3, align = 'hv', axis = "ltb", label_size=18)
save_plot(file.path(out_plots,"SuppFig5.2_GTEx.pdf"), suppFig5.2, ncol = 1, nrow = 3, base_height = 10, base_width = 20)
save_plot(file.path(out_plots,"SuppFig5.2_GTEx.png"), suppFig5.2, ncol = 1, nrow = 3, base_height = 10, base_width = 20)




## ***********************************************************************************************
