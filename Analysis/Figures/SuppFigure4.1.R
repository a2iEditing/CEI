
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
l_EI_wMETA = fread(file.path(input_dir, "EditingIndex.longFormat.wMetadata.csv"), stringsAsFactors = F) %>%
  mutate(disease_or_control = case_when(Disease_Status_manual_annotation == "Control" ~ "Control",
                                        Disease_Status_manual_annotation == "OtherDesign" ~ "OtherDesign",
                                        .default = "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case", "OtherDesign")),
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "GenomeWideAlu", 
                                                    "CEI" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("CEI"))

allData = fread(file.path(input_dir, "AllData.wMetadata.csv"), stringsAsFactors = F) %>%
  mutate(disease_or_control = case_when(Disease_Status_manual_annotation == "Control" ~ "Control",
                                        Disease_Status_manual_annotation == "OtherDesign" ~ "OtherDesign",
                                        .default = "Case"),
         disease_or_control = factor(disease_or_control, levels = c("Control", "Case", "OtherDesign"))) %>%
  select(-Disease_Status_original_annotation, -Disease_Status_original_annotation_columns) %>%
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


# noise check
alu_index_noise_check = list.files("SummaryNonOptimalSubset/EI_results/", pattern = "_AluEditingIndex.csv",full.names = T) %>%
  map_dfr(fread, header = T, stringsAsFactors = F) %>%
  mutate(EI_type = "GenomeWideAlu") %>%
  # add cytoindex
  bind_rows(list.files("SummaryNonOptimalSubset/EI_results/", pattern = "_3UTREditingIndex.csv",full.names = T) %>%
              map_dfr(fread, header = T, stringsAsFactors = F) %>%
              mutate(EI_type = "IRAlu3UTR")) %>%
  mutate(Parameters = "Non-optimal Alignment") %>%
  pivot_longer(cols = contains("EditingIndex"),
               names_to = c("Mismatch", "IndexType"), names_sep = "EditingIndex", values_to = "EditingIndex") %>%
  select(-IndexType)
# merge with optimal aligment for sample subset
alu_index_noise_check = alu_index_noise_check %>% 
  select(all_of(intersect(colnames(alu_index_noise_check), 
                          colnames(l_EI_wMETA))), 
         Parameters) %>%
  bind_rows(l_EI_wMETA %>% 
              mutate(Parameters = "Optimal Alignment") %>%
              select(all_of(intersect(colnames(alu_index_noise_check), 
                                      colnames(l_EI_wMETA))), 
                     Parameters) %>%
              filter(Sample %in% alu_index_noise_check$Sample)) %>%
  mutate(`Editing Index Type` = forcats::fct_recode(EI_type, "AEI" = "GenomeWideAlu", 
                                                    "CEI" = "IRAlu3UTR") %>% 
           forcats::fct_relevel("CEI"))

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

# Supp. Fig. 2A -----------------------------------------------------------------

suppFig4a = l_EI_wMETA %>%
  mutate(index_category = if_else(Mismatch == "A2G", "Signal", "Noise")) %>%
  group_by(Mismatch, index_category, `Editing Index Type`) %>%
  summarise(MeanEditingIndex = mean(EditingIndex),
            SD = sd(EditingIndex)) %>%
  mutate(ymin = if_else(MeanEditingIndex-SD > 0, MeanEditingIndex-SD, 0),
         ymax = MeanEditingIndex+SD) %>%
  ggplot(aes(x = Mismatch, y = MeanEditingIndex, fill = `Editing Index Type`)) +
  geom_col(position = position_dodge(), color= "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), 
                position = position_dodge(width = 0.9), width = 0.3) +
  scale_fill_manual(values = index_colors) +
  scale_y_continuous(expand = c(0.005, 0.01)) +
  theme_custom() +
  # notice this is a trick to create the significance we saw in a previous plot
  ggpubr::stat_compare_means(label = "p.signif", symnum.args = list(cutpoints = c(0, Inf), symbols = c("****")),
                             label.y = -0.07) +
  ggtitle("Mean index per mismatch and index type") +
  ylab("Mean editing index") +
  # ggforce::facet_row(
  #   facets = vars(index_category),
  #   scales = "free",
  #   space  = "free") +
  facet_grid(.~index_category, space = "free", scales = "free") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_blank())


# Supp. Fig. 2B -----------------------------------------------------------------

suppFig4b = alu_index_noise_check %>%
  filter(Mismatch != "A2G") %>%
  group_by(Sample, Parameters, `Editing Index Type`) %>%
  slice_max(order_by = EditingIndex, n = 1, with_ties = F) %>%
  ggplot(aes(x = EditingIndex, color = Parameters)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  ggtitle("Empirical cumulative distribution of noise for different run parameters",
          subtitle = "Line at x = 0.3") + 
  scale_x_log10() +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
  theme_custom() +
  facet_grid(.~`Editing Index Type`) +
  scale_color_manual(values = rev(case_control_colors)) +
  xlab("Noise (maximal non-A2G mismatch)") +
  ylab("Cumulative proportion") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.title = element_blank())


# Supp. Fig. 2C -----------------------------------------------------------------

suppFig4c = l_EI_wMETA %>%
  group_by(Sample, `Editing Index Type`) %>%
  mutate(Rank = dense_rank(desc(EditingIndex))) %>%
  group_by(Mismatch, Rank, `Editing Index Type`) %>%
  tally() %>%
  ungroup() %>%
  complete(Mismatch, Rank, `Editing Index Type`, fill = list(n = 0)) %>%
  mutate(Rank = case_when(Rank == 1 ~ paste0(Rank, "st"),
                          Rank == 2 ~ paste0(Rank, "nd"),
                          Rank == 3 ~ paste0(Rank, "rd"),
                          .default = paste0(Rank, "th"))) %>%
  group_by(Rank, `Editing Index Type`) %>%
  mutate(total = sum(n),
         pcnt = n/total,
         label = if_else(pcnt > 0.05, paste0(n, "\n(", scales::label_percent(accuracy = 0.1)(pcnt), ")"), as.character(n))) %>%
  ggplot(aes(x = Mismatch, y = n, fill = Mismatch)) +
  geom_col(color = "black") +
  geom_text(aes(label = label), vjust=-0.3, size = 3) +
  theme_custom(legend_position = "none") +
  facet_grid(`Editing Index Type`~Rank) +
  ggtitle("Mismatch rank per sample") +
  ylab("Number of samples") +
  scale_fill_manual(values = chosen_colors_extended) +
  expand_limits(y = 17500) +
  theme(axis.title.x = element_blank())




# Join --------------------------------------------------------------------
library(cowplot) 

# Fig4 <- plot_grid(plot_grid(Fig4b, Fig4c,rel_widths = c(0.25, 0.75),
#                             labels=c("B", "C"), ncol = 2, nrow = 1, align = 'hv', label_size=18),
#                   plot_grid(Fig4d, Fig4e,
#                             labels=c("D", "E"), ncol = 2, nrow = 1, align = 'hv', label_size=18),
#                   labels=c("", ""), ncol = 1, nrow = 2, align = 'hv', label_size=18)
suppFig4 <- plot_grid(plot_grid(suppFig4a, #rel_widths = c(1.3, 2),
                                labels=c("A"), ncol = 1, nrow = 1, align = 'hv', label_size=18),
                      plot_grid(suppFig4b, suppFig4c, rel_heights = c(1.5, 2),
                                labels=c("B", "C"), ncol = 1, nrow = 2, align = 'hv', label_size=18, axis = "lb"),
                      # rel_heights = c(1.5, 2),
                      rel_widths = c(1, 3),
                      labels=c("", ""), ncol = 2, nrow = 1, align = 'hv', axis = "lb", label_size=18)
# Fig4 <- plot_grid(Fig4b, Fig4c,
#                   Fig4d, Fig4e,
#                   labels=c("B", "C", "D", "E"), ncol = 2, nrow = 2, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"SuppFig4.1.pdf"), suppFig4, ncol = 1, nrow = 2, base_height = 7, base_width = 20)
save_plot(file.path(out_plots,"SuppFig4.1.png"), suppFig4, ncol = 1, nrow = 2, base_height = 7, base_width = 20)


## ***********************************************************************************************
