## Autoimmune Datasets Analysis Script
## Based on Sambuseq Figure 3 Template

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

### create directories ------------------------------------------------------
dir.create(out_plots, recursive = T)

# Load intermediate (Editing Indices) -------------------------------------
global_editing_index = fread("AEI.csv", stringsAsFactors = F) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Global")) 

cytoplasmic_editing_index = fread("CEI.csv", stringsAsFactors = F) %>%
  rename_with(.cols = contains("EditingIndex"), ~paste0(.x, "_Cytoplasmic"))

editing_data = global_editing_index %>% 
  inner_join(cytoplasmic_editing_index, by = "Sample")

editing_data_long = editing_data %>%
  pivot_longer(cols = contains("EditingIndex"), names_to = c("Mismatch", "EI_type"), names_sep = "_", values_to = "EditingIndex") %>%
  mutate(Mismatch = str_remove(Mismatch, "EditingIndex"),
         "Editing Index Type" = EI_type,
         `Editing Index Type` = forcats::fct_recode(`Editing Index Type`, "AEI" = "Global", 
                                                    "CEI" = "Cytoplasmic")%>% 
           forcats::fct_relevel("CEI"))

editing_data_long_a2g = editing_data_long %>% 
  filter(Mismatch == "A2G")


# Classification (Autoimmune Metadata) ------------------------------------

classification = fread("combined_autoimmune_metadata_info.csv", stringsAsFactors = F) %>%
  mutate(Condition = case_match(Condition, 
                                "Control" ~ "Healthy",
                                "Healthy Healthy" ~ "Healthy",
                                .default = Condition) %>%
           factor(levels = c("Healthy", "Multiple sclerosis",
                             "SICCA", "Sjogren's syndrome", "Rheumatoid Arthritis", "Lupus",
                             "Non-lesional Atopic Dermatitis", "Lesional Atopic Dermatitis", "Chronic lesion Atopic Dermatitis",
                             "Non-lesional Psoriasis","Lesional Psoriasis",
                             "Crohn's disease", "Ulcerative colitis")))

# All data joined with classifications ------------------------------------
allData = editing_data %>%
  inner_join(classification)

allData_long_a2g = editing_data_long_a2g %>%
  inner_join(classification)

# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")

# Supp. Fig. 5.4A ----------------------------------------------------------------
control_color = case_control_colors[1]
case_color = case_control_colors[2]
intermediate_color = chosen_colors[3]
other_colors = chosen_colors_extended[3:5]

# biological conditions 
suppFig5.4a = allData_long_a2g %>%
  filter(`Editing Index Type` == "CEI")%>%
  ggplot(aes(x = Condition, y = EditingIndex, fill = Condition)) +
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  ggh4x::facet_nested(~ Data, scales = "free", space = "free") +
  theme_custom(legend_position = "none") +
  expand_limits(y = 0) +
  scale_fill_manual(values = c(control_color, case_color,
                               intermediate_color, case_color, case_color, case_color,
                               intermediate_color, case_color, other_colors[2], intermediate_color, case_color,
                               case_color, case_color)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Editing index") +
  ggtitle("CEI in autoimmune diseases", 
          subtitle = "Human")    +
  theme(axis.title.x = element_blank())

suppFig5.4b = allData_long_a2g %>%
  filter(`Editing Index Type` == "CEI")%>%
  ggplot(aes(x = Condition, y = EditingIndex, fill = Condition)) +
  geom_boxplot(outliers = F)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), shape = 21) +
  ggh4x::facet_nested(~ Data, scales = "free", space = "free") +
  theme_custom(legend_position = "none") +
  expand_limits(y = 0) +
  scale_fill_manual(values = c(control_color, case_color,
                               intermediate_color, case_color, case_color, case_color,
                               intermediate_color, case_color, other_colors[2], intermediate_color, case_color,
                               case_color, case_color)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Editing index") +
  ggtitle("AEI in autoimmune diseases", 
          subtitle = "Human")    +
  theme(axis.title.x = element_blank())



# Join --------------------------------------------------------------------
library(cowplot)
suppFig3 <- plot_grid(suppFig5.4a,
                      suppFig5.4b,
                      labels="AUTO", ncol = 1, nrow = 2, align = 'hv', axis = "tblr", label_size=18)
save_plot(file.path(out_plots,"SuppFig5.4.pdf"), suppFig3, ncol = 1, nrow = 2, base_height = 8, base_width = 12)
save_plot(file.path(out_plots,"SuppFig5.4.png"), suppFig3, ncol = 1, nrow = 2, base_height = 8, base_width = 12)

