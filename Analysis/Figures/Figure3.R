
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
global_editing_index = map_dfr(c("ExtraDatasets/Summary/AluEditingIndex.csv", 
                                 "MouseDatasets/Summary/AluEditingIndex.csv", 
                                 "FractionationDatasets/Summary/AluEditingIndex.csv"),
                               ~fread(.x, stringsAsFactors = F)) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Global")) 

cytoplasmic_editing_index = map_dfr(c("ExtraDatasets/Summary/UTR3EditingIndex.csv", 
                                      "MouseDatasets/Summary/UTR3EditingIndex.csv", 
                                      "FractionationDatasets/Summary/UTR3EditingIndex.csv"),
                                    ~fread(.x, stringsAsFactors = F)) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Cytoplasmic"))

# tandem_cytoindex = fread("ExtraDatasets/Summary/TandemUTR3EditingIndex.csv", stringsAsFactors = F) %>%
#   rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Tandem"))

tandem_no_inverted_cytoindex = map_dfr(c("ExtraDatasets/Summary/TandemNoInverted3UTREditingIndex.csv", 
                                         "MouseDatasets/Summary/TandemNoInverted3UTREditingIndex.csv", 
                                         "FractionationDatasets/Summary/TandemNoInverted3UTREditingIndex.csv"),
                                       ~fread(.x, stringsAsFactors = F)) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_TandemNoInverted")) 

editing_data = global_editing_index %>% 
  inner_join(cytoplasmic_editing_index) %>%
  # left_join(tandem_cytoindex) %>%
  inner_join(tandem_no_inverted_cytoindex) 

editing_data_long = editing_data %>%
  pivot_longer(cols = contains("EditingIndex"),names_to = c("Mismatch", "EI_type"), names_sep = "_", values_to = "EditingIndex") %>%
  mutate(Mismatch = str_remove(Mismatch, "EditingIndex"),
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "Global", 
                                                    "CEI" = "Cytoplasmic",
                                                    # "Tandem Index (with inverted)" = "Tandem",
                                                    "Tandem Index" = "TandemNoInverted") %>% 
           forcats::fct_relevel("CEI"))

editing_data_long_a2g = editing_data_long %>% 
  filter(Mismatch == "A2G")


# Classification ----------------------------------------------------------


classification = fread("ExtraDatasets/Metadata/PRJNA386593_ADARp150KO_wIFNTreatment.csv", stringsAsFactors = F) %>%
  rename(Sample = Run,
         Genotype = genotype, 
         Treatment = treatment) %>%
  select(Sample, Genotype, Treatment) %>%
  mutate(Name = "ADARp150 KO with IFN Treatment",
         Data = "PRJNA386593") %>%
  bind_rows(fread("ExtraDatasets/Metadata/PRJNA1090949_ADARp150Reconstitution.csv", stringsAsFactors = F) %>%
              inner_join(fread("ExtraDatasets/Metadata/PRJNA1090949_ADARp150Reconstitution.GEO.csv", stringsAsFactors = F)) %>%
              rename(Sample = Run,
                     Genotype = Title,
                     Treatment = treatment) %>%
              select(Sample, Genotype, Treatment) %>%
              filter(Treatment == "None") %>%
              mutate(Name = "ADARp150 Reconstitution",
                     Data = "PRJNA1090949",
                     Genotype = str_remove(Genotype, "[\\+\\-]FCS"))) %>%
  mutate(Treatment = case_match(Treatment, 
                                "None" ~ "Untreated", 
                                "Interferon beta treated" ~ "IFN-beta treated", 
                                "mock treated" ~ "Mock treated") %>%
           factor(levels = c("Untreated", "Mock treated", "IFN-beta treated")),
         Genotype = Genotype %>% #str_replace_all("\\+", "\n+") %>% 
           str_replace_all("wildtype", "Wild type") %>%
           str_trim() %>%
           factor(levels = c("Wild type", "ADAR1p150 KO", "ADAR1 KO",
                             "Non-targeting shRNA", "ADAR1 shRNA + ADARp150", "ADAR1 shRNA + vector only")),
         Data = factor(Data, levels = c("PRJNA386593", "PRJNA1090949")))

BetaCells_Islets_wIFN = fread("ExtraDatasets/Metadata/BetaCells_Islets_wIFN.csv", stringsAsFactors = F) %>%
  filter(SourceName == "EndoC-BH1 cells") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("2h", "8h", "24h", "48h")),
         Treatment = case_when(Dataset == "Cytokines" & Condition  == "Treatment" ~ "Cytokine",
                               Dataset == "IFN-Alpha" & Condition  == "Treatment" ~ "IFN-Alpha Treated",
                               .default = "Control"),
         Data = case_when(Dataset == "Cytokines" ~ "PRJNA564614",
                          Dataset == "IFN-Alpha" ~ "PRJNA550411"))

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
                     Tissue = "Editing Deficient MEFs")) %>%
  mutate(Genotype = Genotype %>%
           case_match("WT"~"Wild type",
                      "KO"~"Adarp110 KO", 
                      "dKO"~"Adar110 KO&\nAdar2 KO", 
                      .default = Genotype) %>%
           factor(levels = c("Wild type", "Adarp110 KO", "Adar110 KO&\nAdar2 KO", 
                             "RFP", "ADAR1p110", "ADAR1p150")))

fractionation_info = fread("FractionationDatasets/Metadata/PRJNA434426_100_SE_Stranded_hg38_hiPCS.csv", stringsAsFactors = F) %>%
  rename(Sample = Run) %>%
  select(Sample, Fraction) %>%
  mutate(Fraction = str_to_title(Fraction) %>%
           factor(levels = c("Cytoplasm", "Nucleus", "Chromatin")),
         Data = "PRJNA434426")

# All data ---------------------------------------------------------------
allData = editing_data %>%
  inner_join(classification)

allData_long_a2g = editing_data_long_a2g %>%
  inner_join(classification)

BetaCells_Islets_wIFN_data = BetaCells_Islets_wIFN %>%
  inner_join(editing_data)

BetaCells_Islets_wIFN_data_long_a2g = BetaCells_Islets_wIFN %>%
  inner_join(editing_data_long_a2g)

mouseData_long_a2g = mouse_info %>%
  inner_join(editing_data_long_a2g)

fractionationData_long_a2g = fractionation_info %>%
  inner_join(editing_data_long_a2g)

# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")


# Fig. 4A ----------------------------------------------------------------

# biological conditions 
fig3a = allData_long_a2g %>%
  mutate(Treatment = forcats::fct_collapse(Treatment, "Untreated \\ mock treated" = c("Untreated", "Mock treated"))) %>%
  ggplot(aes(x = Genotype, y = EditingIndex, fill = Treatment)) +
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  ggh4x::facet_nested(~ Data + `Editing Index Type`, scales = "free", space = "free") +
                      # strip  = ggh4x::strip_nested(
                      #   background_x = ggh4x::elem_list_rect(fill = c(rep("grey", 2), 
                      #                                                 rep("grey85", 6)),
                      #                                        size   = c(rep(1, 2), 
                      #                                                   rep(0.5, 6))))) +
  theme_custom() +
  expand_limits(y = 0) +
  scale_fill_manual(values = c(case_control_colors)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing.x = unit(c(0.25, 0.25, 1, 0.25, 0.25), "lines")) +
  ylab("Editing index") +
  ggtitle("Effect of ADARp150 on editing index", 
          subtitle = "Human")    +
  theme(axis.title.x = element_blank())

# Fig. 4B ----------------------------------------------------------------

# compute std for ratios
curr_data_fig3b = allData_long_a2g %>%
  filter(Data == "PRJNA386593") %>%
  mutate(Treatment = if_else(grepl(x = Treatment, "None|Mock"), "Control", "Treatment"),
         `Editing Index Type` = str_replace(`Editing Index Type`, 
                                            pattern = " ", replacement = "\n"))

ratio_table_fig3b = curr_data_fig3b %>%
  group_by(Genotype, Treatment, Data, Name, `Editing Index Type`) %>%
  summarise(EditingIndex = mean(EditingIndex)) %>%
  pivot_wider(names_from = Treatment, values_from = EditingIndex) %>%
  filter(!is.na(Treatment)) %>%
  mutate(Ratio = Treatment / Control)

# Compute sd of ratios
stats_fig3b <- curr_data_fig3b %>%
  # group by each condition, genotype and editing index
  group_by(Genotype, Treatment, Data, Name, `Editing Index Type`) %>%
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
  group_by(Genotype, Data, Name, `Editing Index Type`) %>%
  # compute the square root of the sum of normalized sd squared
  summarise(sqrt_sum_of_norm_sd_editing_squared = sqrt(sum(final))) %>%
  # add overall mean per group
  inner_join(ratio_table_fig3b %>%
               ungroup() %>%
               select(Genotype, Data, Name, `Editing Index Type`, Ratio)) %>%
  mutate(final_sd = sqrt_sum_of_norm_sd_editing_squared * Ratio)


# ratio between conditions
fig3b = ratio_table_fig3b %>%
  inner_join(stats_fig3b) %>%
  ggplot(aes(x = Genotype, y = Ratio, fill = `Editing Index Type`)) +
  geom_col(position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = Ratio-final_sd, ymax = Ratio+final_sd), position = position_dodge(width = 0.9), width = 0.2) +
  facet_grid(.~ Data + Name, scales = "free", space = "free") +
  theme_custom() +
  expand_limits(y = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = c(index_colors, "white", "#6E8B3D")) +
  guides(fill = guide_legend(title = "Editing\nIndex\nType")) +
  ggtitle("Effect of ADARp150 on IFN-induced editing",
          subtitle = "Human") +
  ylab("IFN-treated to mock-treated index ratio") +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank())

# Fig. 4C ----------------------------------------------------------------
fig3c =  BetaCells_Islets_wIFN_data_long_a2g %>%
  ggplot(aes(x = Timepoint, y = EditingIndex, fill = `Editing Index Type`)) +
  geom_boxplot(outliers = F) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  ggh4x::facet_nested(~ Data + Treatment, scales = "free", space = "free") +
                      # strip  = ggh4x::strip_nested(
                      #   background_x = ggh4x::elem_list_rect(fill = c(rep("grey", 2), 
                      #                                                 rep("grey85", 6)),
                      #                                        size   = c(rep(1, 2), 
                      #                                                   rep(0.5, 6))))) +
  theme_custom() +
  theme(panel.spacing.x = unit(c(0.25, 1, 0.25), "lines")) +
  expand_limits(y = 0) +
  scale_fill_manual(values =  c(index_colors, "white", "#6E8B3D")) +
  ylab("Editing index") +
  ggtitle("Effect of IFN-treatment on editing index",
          subtitle = "Human") +
  theme(legend.title = element_blank())


# Fig. 4D ----------------------------------------------------------------
ratio_table_fig3d = BetaCells_Islets_wIFN_data_long_a2g %>%
  group_by(Timepoint, Condition, Data, `Editing Index Type`) %>%
  summarise(EditingIndex = mean(EditingIndex)) %>%
  pivot_wider(names_from = Condition, values_from = EditingIndex) %>%
  mutate(Ratio = Treatment / Control)

# Compute sd of ratios
stats_fig3d <- BetaCells_Islets_wIFN_data_long_a2g %>%
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
  inner_join(ratio_table_fig3d %>%
               ungroup() %>%
               select(-Control, -Treatment)) %>%
  mutate(final_sd = sqrt_sum_of_norm_sd_editing_squared * Ratio)



# ratio between conditions
fig3d = ratio_table_fig3d %>%
  inner_join(stats_fig3d) %>%
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


# Fig. 4E -----------------------------------------------------------------
fig3e = mouseData_long_a2g %>%
  ggplot(aes(x = Genotype, y = EditingIndex, fill = `Editing Index Type`)) +
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  ggh4x::facet_nested(~ Data + Tissue, scales = "free", space = "free") +
  theme_custom() +
  theme(panel.spacing.x = unit(c(0.25, 1), "lines")) +
  expand_limits(y = 0) +
  scale_fill_manual(values = c(index_colors, "white")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),) +
  ylab("Editing index") +
  ggtitle("Effect of Adarp150 on editing index", 
          subtitle = "Mouse")  +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1, 0))


# Fig. 4F -----------------------------------------------------------------
fig3f = fractionationData_long_a2g %>%
  ggplot(aes(x = Fraction, y = EditingIndex, fill = `Editing Index Type`)) +
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  facet_grid(~ Data, scales = "free", space = "free") +
  theme_custom(legend_position = "none") +
  expand_limits(y = 0) +
  scale_fill_manual(values = c(index_colors, "white")) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),) +
  ggtitle("Effect of cell fraction on editing index", 
          subtitle = "Human") +
  ylab("Editing index") +
  theme(legend.title = element_blank())

# Join --------------------------------------------------------------------
library(cowplot)
fig3 <- plot_grid(plot_grid(fig3a, fig3b, rel_widths = c(2, 1),
                            labels=c("A", "B"), ncol = 2, nrow = 1, align = 'hv', axis = "lb", label_size=18),
                  plot_grid(fig3c, fig3d,
                            labels=c("C", "D"), ncol = 2, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  plot_grid(fig3e, fig3f, rel_widths = c(2, 1),
                            labels=c("E", "F"), ncol = 2, nrow = 1, align = 'hv', label_size=18, axis = "lb"),
                  rel_heights = c(1.2, 1, 1),
                  labels=c("", "", ""), ncol = 1, nrow = 3, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"Figure3.pdf"), fig3, ncol = 1, nrow = 3, base_height = 7, base_width = 14)
save_plot(file.path(out_plots,"Figure3.png"), fig3, ncol = 1, nrow = 3, base_height = 7, base_width = 14)




