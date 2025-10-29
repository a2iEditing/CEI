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
out_plots="plots"
out_stats="stats"

### create directories ------------------------------------------------------
#create directories, if do not exist
dir.create(out_plots, recursive = T)
dir.create(out_stats, recursive = T)


# Load intermediate -------------------------------------------------------
metadata = fread("Metadata/SraRunTable_GTEX_RNAseq.filtered_noTumor_noStranded_paired.minimal.csv") %>%
  mutate(CohortOriginal = Cohort, 
         Cohort = if_else(Cohort == "Postmortem", Cohort, "Antemortem")) 

gtex_colors = fread("GTEx/Metadata/gtex_colors.csv") %>%
  { setNames(pull(., Color), pull(., Tissue)) }

performance_analysis = fread("GTEx/Summary/EIRunStats.csv")

site_count = fread("EI_Resources/BaseCounts.hg38.Alu3pUTR_minLen200_17022021.InvertedRepeatsIn3pUTR.sorted.merged.csv", stringsAsFactors = F) %>%
  mutate(EI_type = "IRAlu3UTR") %>%
  bind_rows(fread("EI_Resources/ucscHg38Alu.csv", stringsAsFactors = F) %>%
              mutate(EI_type = "GenomeWideAlu")) %>%
  group_by(EI_type) %>%
  summarise(across(c("A", "C", "G", "T"), list(sum = sum))) %>%
  mutate(AT_site_count = A_sum + T_sum) %>%
  select(EI_type, AT_site_count)%>%
  pivot_wider(names_from = "EI_type", values_from = AT_site_count)

allData = fread("GTEx/Summary/AluEditingIndex.csv", stringsAsFactors = F) %>%
  inner_join(fread("GTEx/Summary/UTR3EditingIndex.csv", stringsAsFactors = F),
             by = "Sample", suffix = c("GenomeWideAlu", "IRAlu3UTR")) 


l_EI_wMETA = allData %>%
  select(Sample, contains("EditingIndex")) %>%
  pivot_longer(cols = contains("EditingIndex"), values_to = "EditingIndex", names_to = c("Mismatch", "EI_type"), names_sep = "EditingIndex") %>%
  mutate(`Editing Index Type` = case_match(EI_type, 
                                           "GenomeWideAlu" ~ "Global Index",
                                           "IRAlu3UTR" ~ "Cytoplasmic Index") %>%
           forcats::fct_relevel("Cytoplasmic Index")) %>%
  inner_join(metadata, by = join_by(Sample == "Run"))

allData = allData %>%
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
         
         # normalized coverage
         NormalizedCoveragePerSite_GenomeWideAlu = (NumOfAGenomeWideAlu + NumOfA2GMismatchesGenomeWideAlu + NumOfTGenomeWideAlu + NumOfT2CMismatchesGenomeWideAlu) / site_count$GenomeWideAlu,
         NormalizedCoveragePerSite_IRAlu3UTR = (NumOfAIRAlu3UTR + NumOfA2GMismatchesIRAlu3UTR + NumOfTIRAlu3UTR + NumOfT2CMismatchesIRAlu3UTR) / site_count$IRAlu3UTR) %>%
  inner_join(metadata, by = join_by(Sample == "Run"))

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

l_EI_wMETA_a2g_filtered = l_EI_wMETA_a2g %>%
  filter(Sample %in% allData_filtered$Sample)

allData_filtered_mergedReps = allData_filtered %>%
  group_by(Donor, Tissue, TissueHistologicalType, Sex, Cohort, CohortOriginal, Age, AgeGroup, BMI) %>%
  summarise(across(contains("EditingIndex"), mean),
            ReplicateNum = n()) 

l_EI_wMETA_a2g_filtered_mergedReps = l_EI_wMETA_a2g_filtered %>%
  group_by(EI_type, `Editing Index Type`, Donor, Tissue, TissueHistologicalType, Sex, Cohort, CohortOriginal, Age, AgeGroup, BMI) %>%
  summarise(EditingIndex = mean(EditingIndex),
            ReplicateNum = n()) 



# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")

# Supp. Fig. 3A -----------------------------------------------------------------

correlation_coef_fig3b = coef(lm(allData$A2GEditingIndexIRAlu3UTR ~ 0 + allData$A2GEditingIndexGenomeWideAlu))
fig3b = allData %>%
  ggplot(aes(x = A2GEditingIndexGenomeWideAlu, y = A2GEditingIndexIRAlu3UTR)) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  theme_custom() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_abline(slope = correlation_coef_fig3b, intercept = 0, color = "#8B1A1A") +
  scale_x_continuous(breaks = c(0:4)) +
  scale_y_continuous(breaks = c(0:8)) +
  ggtitle("Correlation of editing signal",
          subtitle = "Dashed line at x = y") +
  xlab("Global A2G editing index") +
  ylab("Cytoplasmic A2G editing index") +
  labs(color = "Sample density") +
  theme(legend.text  = element_text(angle = 45, hjust = 1)) +
  # annotate slope
  annotate("text",x=2,y=0.5, size = 5, hjust = 0, label=(paste0("Slope==",signif(correlation_coef_fig3b, digits = 5))),parse=TRUE, color = chosen_colors[2])+ 
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 8), expand = FALSE, default = FALSE, clip = "on")

# Supp. Fig. 3B -----------------------------------------------------------------
temp = allData %>%
  filter(Signal2NoiseRatio_IRAlu3UTR < Inf)
correlation_coef_fig3c = coef(lm(temp$Signal2NoiseRatio_IRAlu3UTR ~ 0 + temp$Signal2NoiseRatio_GenomeWideAlu))
fig3c = temp %>%
  ggplot(aes(x=Signal2NoiseRatio_GenomeWideAlu,
             y=Signal2NoiseRatio_IRAlu3UTR))+
  ggpointdensity::geom_pointdensity(adjust = 4) +
  viridis::scale_color_viridis() +
  theme_custom() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_abline(slope = correlation_coef_fig3c, intercept = 0, color = "#8B1A1A") +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  scale_y_continuous(breaks = seq(0, 400, 25)) +
  ggtitle("Correlation of signal-to-noise ratios",
          subtitle = "Dashed line at x = y") +
  xlab("Global index SNR") +
  ylab("Cytoplasmic index SNR") +
  labs(color = "Sample density") +
  theme(legend.text  = element_text(angle = 45, hjust = 1)) +
  # annotate slope (forct intercept to be 0)
  annotate("text",x=45,y=375, size = 5, hjust = 0, label=(paste0("Slope==",signif(correlation_coef_fig3c, digits = 5))),parse=TRUE, color = chosen_colors[2])+
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 400), expand = FALSE, default = FALSE, clip = "on")

# Supp. Fig. 3C -----------------------------------------------------------
correlation_coef_fig3d = coef(lm(allData$NormalizedCoveragePerSite_IRAlu3UTR ~ 0 + allData$NormalizedCoveragePerSite_GenomeWideAlu))
fig3d = allData %>%
  ggplot(aes(x = NormalizedCoveragePerSite_GenomeWideAlu, y = NormalizedCoveragePerSite_IRAlu3UTR)) +
  ggpointdensity::geom_pointdensity() +
  viridis::scale_color_viridis() +
  theme_custom() +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_abline(slope = correlation_coef_fig3d, intercept = 0, color = "#8B1A1A") +
  scale_x_continuous(breaks = seq(0, 6, 1)) +
  scale_y_continuous(breaks = seq(0, 55, 5)) +
  ggtitle("Correlation of avg. coverage per site",
          subtitle = "Dashed line at x = y") +
  xlab("Average coverage per site in global index") +
  ylab("Average coverage per site in cytoplasmic index") +
  labs(color = "Sample density")  +
  theme(legend.text  = element_text(angle = 45, hjust = 1)) +
  # annotate slope
  annotate("text",x=3.7,y=52, size = 5, hjust = 0, label=(paste0("Slope==",signif(correlation_coef_fig3d, digits = 3))),parse=TRUE, color = chosen_colors[2])+
  coord_cartesian(xlim = c(0, 6), ylim = c(0, 55), expand = FALSE, default = FALSE, clip = "on")




# Supp. Fig. 3D -----------------------------------------------------------------
fig3e = l_EI_wMETA_a2g_filtered %>%
  group_by(Tissue, `Editing Index Type`) %>%
  summarise(across(c(EditingIndex), list(sd = sd,
                                         mean = mean))) %>% 
  mutate(`Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "Global\nIndex" = "Global Index", 
                                                    "Cytoplasmic\nIndex" = "Cytoplasmic Index") %>% 
           forcats::fct_relevel("Cytoplasmic\nIndex")) %>%
  mutate("Normalized standard deviation" = EditingIndex_sd / EditingIndex_mean) %>%
  ggplot(aes(x = `Editing Index Type`, y = `Normalized standard deviation`, fill = `Editing Index Type`)) +
  geom_violin(draw_quantiles = 0.5) +
  geom_jitter(size = 2, width = 0.05) +
  theme_custom(legend_position = "none") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12)) +
  ggtitle("Normalized standard deviation", subtitle = "Per tissue") +
  ylab("Standard deviation / mean") +
  scale_fill_manual(values = c(index_colors, case_control_colors, "grey")) +
  ggpubr::stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2)), paired = T) 


# Supp. Fig. 3E -----------------------------------------------------------------

fig3f = l_EI_wMETA_a2g %>%
  ggplot(aes(x = NextBestMM, color = `Editing Index Type`)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  ggtitle("Cumulative distribution of noise",
          subtitle = "Dashed line at x = 0.3") +
  scale_x_log10() +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  theme_custom() +
  scale_color_manual(values = index_colors) +
  xlab("Noise (maximum non-A2G mismatch)") +
  ylab("Cumulative proportion") +
  scale_y_continuous(labels = scales::percent_format())

fig3f_inset = l_EI_wMETA_a2g %>% 
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

# Supp. Fig. 3F ---------------------------------------------------------------
# runtime
fig3g = performance_analysis %>%
  group_by(EI_type) %>%
  summarise(mean = mean(`Total elapsed time (seconds)`/60), sd = sd(`Total elapsed time (seconds)`/60)) %>%
  mutate(ymin = if_else(mean-sd > 0, mean-sd, 0),
         ymax = mean+sd,
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "Global\nIndex" = "GenomeWideAlu", 
                                                    "Cytoplasmic\nIndex" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("Cytoplasmic\nIndex")) %>%
  ggplot(aes(x = `Editing Index Type`, y = mean, fill = `Editing Index Type`)) +
  geom_col()+
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                position = position_dodge(width = 0.9), width = 0.3) +
  scale_fill_manual(values = index_colors) +
  theme_custom(legend_position = "none") +
  theme(axis.title.x = element_blank()) +
  ylab("Total elapsed time (minutes)") +
  scale_y_continuous(labels = scales::label_timespan(unit = "mins"), breaks = c(0:20)) +
  ggtitle("Runtime", subtitle = "Average across samples") 


# Supp. Fig. 3G -----------------------------------------------------------------
# memory
fig3h = performance_analysis %>%
  group_by(EI_type) %>%
  summarise(mean = mean(`Maximum resident set size (kbytes)`*1024), sd = sd(`Maximum resident set size (kbytes)`*1024)) %>%
  mutate(ymin = if_else(mean-sd > 0, mean-sd, 0),
         ymax = mean+sd,
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "Global\nIndex" = "GenomeWideAlu", 
                                                    "Cytoplasmic\nIndex" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("Cytoplasmic\nIndex")) %>%
  ggplot(aes(x = `Editing Index Type`, y = mean, fill = `Editing Index Type`)) +
  geom_col()+
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                position = position_dodge(width = 0.9), width = 0.3) +
  scale_fill_manual(values = index_colors) +
  theme_custom(legend_position = "none") +
  theme(axis.title.x = element_blank()) +
  ylab("Maximum resident set size (gigabytes)") +
  scale_y_continuous(labels = scales::label_bytes(units = "GB"), breaks = seq(0, 10000000000, 1000000000)) +
  ggtitle("Peak memory usage", subtitle = "Average across samples") 






# Join --------------------------------------------------------------------
library(cowplot) 

fig3 <- plot_grid(plot_grid(fig3b, fig3c, fig3d, fig3e, rel_widths = c(1, 1, 1, 1), 
                            labels=c("A", "B", "C", "D"), ncol = 4, nrow = 1, align = 'hv', label_size=18),
                  plot_grid(fig3f + 
                              # inset
                              patchwork::inset_element(fig3f_inset,
                                                       left   = 0.65, right = 0.98,
                                                       bottom = 0.05, top   = 0.6),
                            fig3g,  fig3h, rel_widths = c(2.2, 0.65, 0.65),
                            labels=c("E", "F", "G"), ncol = 3, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  rel_heights = c(1, 0.8), 
                  labels=c("", ""), ncol = 1, nrow = 2, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"SuppFig3.2_GTEx.pdf"), fig3, ncol = 1, nrow = 2, base_height = 8, base_width = 16)
save_plot(file.path(out_plots,"SuppFig3.2_GTEx.png"), fig3, ncol = 1, nrow = 2, base_height = 8, base_width = 16)

