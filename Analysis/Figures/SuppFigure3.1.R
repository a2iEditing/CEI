
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
global_editing_index = fread("ExtraDatasets/Summary/AluEditingIndex.csv", stringsAsFactors = F) %>%
  bind_rows(fread("MouseDatasets/Summary/AluEditingIndex.csv", stringsAsFactors = F)) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Global")) 

cytoplasmic_editing_index = fread("ExtraDatasets/Summary/UTR3EditingIndex.csv", stringsAsFactors = F) %>%
  bind_rows(fread("MouseDatasets/Summary/UTR3EditingIndex.csv", stringsAsFactors = F)) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Cytoplasmic"))

tandem_no_inverted_cytoindex = fread("ExtraDatasets/Summary/TandemNoInverted3UTREditingIndex.csv", stringsAsFactors = F) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_TandemNoInverted")) 

editing_data = global_editing_index %>% 
  inner_join(cytoplasmic_editing_index) %>%
  left_join(tandem_no_inverted_cytoindex)

editing_data_long = editing_data %>%
  pivot_longer(cols = contains("EditingIndex"),names_to = c("Mismatch", "EI_type"), names_sep = "_", values_to = "EditingIndex") %>%
  mutate(Mismatch = str_remove(Mismatch, "EditingIndex"),
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "Global", 
                                                    "CEI" = "Cytoplasmic",
                                                    "Tandem Index" = "TandemNoInverted") %>% 
           forcats::fct_relevel("CEI"))

editing_data_long_a2g = editing_data_long %>% 
  filter(Mismatch == "A2G")


# Classification ----------------------------------------------------------

BetaCells_Islets_wIFN = fread("ExtraDatasets/Metadata/BetaCells_Islets_wIFN.csv", stringsAsFactors = F) %>%
  filter(SourceName == "Islets of Langharns") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("2h", "8h", "18h", "48h")),
         Treatment = case_when(Dataset == "Cytokines" & Condition  == "Treatment" ~ "Cytokine",
                               Dataset == "IFN-Alpha" & Condition  == "Treatment" ~ "IFN-alpha treated",
                               .default = "Control"),
         Data = case_when(Dataset == "Cytokines" ~ "PRJNA427190",
                          Dataset == "IFN-Alpha" ~ "PRJNA622979") %>%
           factor(levels = c("PRJNA622979", "PRJNA427190"))) 

mouse_info = fread("MouseDatasets/Metadata/DRP006377_Adarp110KO_Adar2KO_PE_100_Unstranded.csv", stringsAsFactors = F) %>%
  rename(Sample = Run) %>%
  select(Sample, Sample_name) %>%
  separate(Sample_name, into = c("Genotype","Tissue"), sep = "_") %>%
  mutate(Name = "Adarp110KO vs. dKO (+ADAR2)",
         Data = "DRP006377", 
         Tissue = case_match(Tissue, 
                             "Br" ~ "Brain",
                             "Thy" ~ "Thymus")) %>%
  bind_rows(fread("MouseDatasets/Metadata/PRJNA780507_p110_p1150_transfection_PE_150_144Real_Stranded.csv", stringsAsFactors = F) %>%
              rename(Sample = Run,
                     Genotype = transfected_with_plasmid) %>%
              select(Sample, Genotype) %>%
              mutate(Genotype = str_remove(Genotype, pattern = ": "),
                     Data = "PRJNA780507",
                     Tissue = "MEFs")) %>%
  mutate(Genotype = Genotype %>%
           forcats::fct_relevel("KO")%>%
           forcats::fct_relevel("RFP")%>%
           forcats::fct_relevel("WT"))

# All data ---------------------------------------------------------------
BetaCells_Islets_wIFN_data = BetaCells_Islets_wIFN %>%
  inner_join(editing_data)

BetaCells_Islets_wIFN_data_long_a2g = BetaCells_Islets_wIFN %>%
  inner_join(editing_data_long_a2g)

mouseData_long = mouse_info %>%
  inner_join(editing_data_long %>%
               # no tandem computed for mouse
               filter(`Editing Index Type` != "Tandem Index"))

mouseData_long_a2g = mouse_info %>%
  inner_join(editing_data_long_a2g %>%
               # no tandem computed for mouse
               filter(`Editing Index Type` != "Tandem Index"))


# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")


# Supp. Fig. 4A ----------------------------------------------------------------

suppFig3a =  BetaCells_Islets_wIFN_data_long_a2g %>%
  ggplot(aes(x = Timepoint, y = EditingIndex, fill = `Editing Index Type`)) +
  geom_boxplot(outliers = F) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  ggh4x::facet_nested(~ Data + Treatment, scales = "free", space = "free") +
  theme_custom() +
  theme(panel.spacing.x = unit(c(0.25, 1, 0.25), "lines")) +
  expand_limits(y = 0) +
  scale_fill_manual(values = c(index_colors, "white")) +
  ylab("Editing index") +
  ggtitle("Effect of IFN-treatment on editing index") +
  theme(legend.title = element_blank())


# Supp. Fig. 4B ----------------------------------------------------------------
ratio_table_suppFig3b = BetaCells_Islets_wIFN_data_long_a2g %>%
  group_by(Timepoint, Condition, Data, `Editing Index Type`) %>%
  summarise(EditingIndex = mean(EditingIndex)) %>%
  pivot_wider(names_from = Condition, values_from = EditingIndex) %>%
  mutate(Ratio = Treatment / Control)

# Compute sd of ratios
stats_suppFig3b <- BetaCells_Islets_wIFN_data_long_a2g %>%
  # group by each condition, genotype and editing index
  group_by(Timepoint, Condition, Data, `Editing Index Type`) %>%
  # mean and sd of each group
  summarise(mean_editing = mean(EditingIndex, na.rm = TRUE),
            sd_editing   = sd(EditingIndex,   na.rm = TRUE),
            n = n()) %>%
  # compute normalized sd and normalized sd squared
  mutate(norm_sd_editing = sd_editing / mean_editing,
         norm_sd_editing_squared = norm_sd_editing^2,
         # divide by n-1
         final = norm_sd_editing_squared / n) %>%
  # group by same, NOT by treatment
  group_by(Timepoint, `Editing Index Type`) %>%
  # compute the square root of the sum of normalized sd squared
  summarise(sqrt_sum_of_norm_sd_editing_squared = sqrt(sum(final))) %>%
  # add overall mean per group
  inner_join(ratio_table_suppFig3b %>%
               ungroup() %>%
               select(-Control, -Treatment)) %>%
  mutate(final_sd = sqrt_sum_of_norm_sd_editing_squared * Ratio)



# ratio between conditions
suppFig3b = ratio_table_suppFig3b %>%
  inner_join(stats_suppFig3b) %>%
  ggplot(aes(x = Timepoint, y = Ratio, fill = `Editing Index Type`)) +
  geom_col(position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = Ratio-final_sd, ymax = Ratio+final_sd), position = position_dodge(width = 0.9), width = 0.2) +
  facet_grid(~ Data, scales = "free", space = "free") +
  theme_custom() +
  theme(panel.spacing.x = unit(c(1), "lines")) +
  expand_limits(y = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values =  c(index_colors, "white")) +
  ggtitle("Effect of cytokine and IFN treatments on cytoplasmic editing",
          subtitle = "Human") +
  ylab("IFN-treated to control index ratio") +
  theme(legend.title = element_blank())


# Supp. Fig. 4C ----------------------------------------------------------------
suppFig3c = mouseData_long %>%
  ggplot(aes(x = Mismatch, y = EditingIndex, fill = Mismatch)) +
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  # ggh4x::facet_nested(`Editing Index Type`~ Data + Tissue, scales = "free", space = "free") +
  facet_grid(Data + Tissue~`Editing Index Type`, space = "free") +
  theme_custom(legend_position = "none") +
  scale_y_log10() +
  # theme(panel.spacing.x = unit(c(0.25, 1), "lines")) +
  expand_limits(y = 0) +
  scale_fill_manual(values = chosen_colors_extended) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),) +
  ylab("Editing index") +
  guides(fill = guide_legend(nrow = 1)) +
  ggtitle("Editing signal", 
          subtitle = "Mouse") +
  theme(axis.title.x = element_blank())

# Join --------------------------------------------------------------------
library(cowplot)
suppFig3 <- plot_grid(plot_grid(suppFig3a, suppFig3b,
                                labels="AUTO", ncol = 2, nrow = 1, align = 'hv', axis = "tblr", label_size=18),
                      plot_grid(suppFig3c,
                                labels=c("C"), ncol = 1, nrow = 1, align = 'hv', axis = "tblr", label_size=18), 
                      rel_heights = c(1.2, 2),
                      labels="", ncol = 1, nrow = 2, align = 'hv', axis = "tblr", label_size=18)
save_plot(file.path(out_plots,"SuppFig3.1.pdf"), suppFig3, ncol = 1, nrow = 2, base_height = 8, base_width = 14)
save_plot(file.path(out_plots,"SuppFig3.1.png"), suppFig3, ncol = 1, nrow = 2, base_height = 8, base_width = 14)




