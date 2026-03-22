##sambusaq chronic vs acute####
##sambuseq_all_control_vs_disaese
####sambuseq control vs diseases

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


# Load intermediate -------------------------------------------------------
cytoplasmic_editing_index = fread("PRJNA386593_ADARp150KO_wIFNTreatment_aluelement_length_treatment_genotype_tissue.csv", stringsAsFactors = F) %>%
  filter(Genotype == "wildtype") %>%
  mutate(Treatment = case_match(Treatment, 
                                "Interferon beta treated" ~ "IFN-beta treated", 
                                "mock treated" ~ "Mock treated") %>%
           factor(levels = c("Mock treated", "IFN-beta treated")))

# Plot --------------------------------------------------------------------
source("custom_theme_and_colors.R")

# Supp. Fig. 3.2A ------------------------------------------------------------

suppFig3.2a = cytoplasmic_editing_index %>%
  mutate(total = sum(IndexedMismatchesOfA2G),
         contrib= IndexedMismatchesOfA2G / total,
         AluElement = tidytext::reorder_within(AluElement, contrib, Treatment)) %>%
  ggplot(aes(x = forcats::fct_reorder(AluElement, contrib), y = contrib, color = Treatment)) +
  geom_point() +
  facet_grid(.~Treatment, space = "free", scale = "free") +
  theme_custom(legend_position = "none") +
  ggtitle("Contribution of elements to CEI",
          subtitle = "Pooled data") +
  labs(fill = "Number of Samples",
       x = expression("CEI " * italic("Alu") * " element"),
       y = "A2G Mismatches") +
  scale_color_manual(values = case_control_colors) +
  tidytext::scale_x_reordered() +
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_log10(labels = scales::label_percent(drop0trailing = TRUE))


# Supp. Fig. B ------------------------------------------------------------
# coverage
covered_elements = cytoplasmic_editing_index %>%
  select(AluElement, Treatment, A2GEditingIndex, TotalCoverageAtAPositions) %>%
  pivot_wider(names_from = Treatment, values_from = c(A2GEditingIndex, TotalCoverageAtAPositions)) %>%
  mutate(across(starts_with("TotalCoverageAtAPositions"), ~replace_na(.x, 0)),
         A2GEditingIndex_IFN_to_mock_ratio = `A2GEditingIndex_IFN-beta treated`/(`A2GEditingIndex_Mock treated`+0.1),
         EnoughCoverage = `TotalCoverageAtAPositions_IFN-beta treated` >= 50 & `TotalCoverageAtAPositions_Mock treated` >= 50) %>% 
  filter(EnoughCoverage) %>%
  pull(AluElement)

suppFig3.2b = cytoplasmic_editing_index %>%
  filter(AluElement%in%covered_elements) %>%
  select(AluElement, Treatment, A2GEditingIndex) %>%
  pivot_wider(names_from = Treatment, values_from = A2GEditingIndex) %>%
  mutate(A2GEditingIndex_IFN_to_mock_ratio = (`IFN-beta treated`+0.1)/(`Mock treated`+0.1),
         Direction = if_else( `IFN-beta treated` > `Mock treated`, "Increased in IFN treatment",
                             if_else(`Mock treated` == 0 & `IFN-beta treated` == 0, "Unchanged", "Decreased in IFN treatment"))) %>%
  ggplot(aes(x = forcats::fct_reorder(AluElement, A2GEditingIndex_IFN_to_mock_ratio), y = A2GEditingIndex_IFN_to_mock_ratio, color = Direction)) +
  geom_point() +
  theme_custom() +
  ggtitle("Changes in editing of CEI elements during IFN treatment",
          subtitle = "For elements with pooled coverage >= 50 in both groups") +
  labs(x = expression("CEI " * italic("Alu") * " element"),
       y = "IFN A2G Editing / Mock A2G Editing") +
  scale_color_manual(values = chosen_colors) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
        # axis.text.x = element_blank(),
        axis.text.x = element_blank() ,
        axis.ticks.x = element_blank()) +
  theme(legend.title = element_blank())


# Supp. Fig. C ------------------------------------------------------------
psuedo_index = cytoplasmic_editing_index %>%
  filter(AluElement%in%covered_elements) %>%
  mutate(TotalCoverageAtAPositions = IndexedCanonicalOfA2G+IndexedMismatchesOfA2G) %>%
  select(AluElement, Treatment, A2GEditingIndex, TotalCoverageAtAPositions) %>%
  pivot_wider(names_from = Treatment, values_from = c(A2GEditingIndex, TotalCoverageAtAPositions)) %>%
  # compute mismatches derived from other group coverage 
  mutate(psuedo_mm1 = `A2GEditingIndex_IFN-beta treated`*`TotalCoverageAtAPositions_Mock treated`/100,
         psuedo_mm2 = `A2GEditingIndex_Mock treated`*`TotalCoverageAtAPositions_IFN-beta treated`/100) %>%
  summarise(across(starts_with("TotalCoverageAtAPositions") | starts_with("psuedo_mm"), sum)) %>%
  mutate(psuedo_index1 = 100 * psuedo_mm1 / (`TotalCoverageAtAPositions_Mock treated` + psuedo_mm1),
         psuedo_index2 = 100 * psuedo_mm2 / (`TotalCoverageAtAPositions_IFN-beta treated` + psuedo_mm2)) %>%
  select(starts_with("psuedo_index")) %>%
  pivot_longer(cols = starts_with("psuedo_index"), names_to = "Treatment", values_to = "index") 

pooled_index = cytoplasmic_editing_index %>%
  filter(AluElement%in%covered_elements) %>%
  group_by(Treatment) %>%
  summarise(across(c(IndexedMismatchesOfA2G, IndexedCanonicalOfA2G), sum)) %>%
  mutate(index = 100 * IndexedMismatchesOfA2G / (IndexedCanonicalOfA2G + IndexedMismatchesOfA2G)) %>%
  select(Treatment, index)

suppFig3.2c = psuedo_index%>%
  bind_rows(pooled_index) %>%
  mutate(Type = if_else(grepl(x = Treatment, "psuedo_index"), "Pseudo-index", "CEI"),
         Treatment = case_match(Treatment,
                                "psuedo_index1" ~ "IFN editing and\nmock expression",
                                "psuedo_index2" ~ "Mock editing and\nIFN expression",
                                "Mock treated" ~ "CEI, Mock treated",
                                "IFN-beta treated" ~ "CEI, IFN-beta treated") %>%
           factor(levels = c("CEI, Mock treated",
                             "Mock editing and\nIFN expression",
                             "IFN editing and\nmock expression",
                             "CEI, IFN-beta treated"))) %>%
  ggplot(aes(x = Treatment, y = index, fill = Treatment)) +
  geom_col(color = "black") +
  theme_custom(legend_position = "none") +
  ggtitle("CEI and psuedo-editing index of pooled data per element",
          subtitle = "For elements with pooled coverage >= 50 in both groups") +
  labs(y = "Index") +
  scale_fill_manual(values = chosen_colors) +
  theme(axis.title.x = element_blank())



# Join --------------------------------------------------------------------
library(cowplot)
suppFig3 <- plot_grid(plot_grid(suppFig3.2a, # rel_heights = c(0.8, 1), 
                                labels=c("A"), ncol = 1, nrow = 1, align = 'hv', label_size=18), 
                      plot_grid(suppFig3.2b, suppFig3.2c, 
                                # rel_widths = c(1.2, 1),
                                labels=c("B", "C"), ncol = 2, nrow = 1, align = 'hv', label_size=18), 
                  # rel_heights = c(0.8, 1), 
                  labels=c("", ""), ncol = 1, nrow = 2, align = 'hv', label_size=18)
save_plot(file.path(out_plots,"SuppFig3.2.pdf"), suppFig3, ncol = 2, nrow = 3, base_height = 4, base_width = 8)
save_plot(file.path(out_plots,"SuppFig3.2.revision.png"), suppFig3, ncol = 2, nrow = 3, base_height = 4, base_width = 8)



# *****************************************************************************************
