rm(list = ls(all = TRUE))
# SETUP ===========================================================================================
# *************************************************************************************************

library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(GEOquery)



# Setup -------------------------------------------------------------------
### Output directory --------------------------------------------------------
out_boxplots = "boxplots"
out_scatterplots = "scatterplots"
out_correlations = "correlations"
out_barplots = "barplots"
out_stats="stats"
out_tables = "tables"


### Consts ------------------------------------------------------------------

#### designs -----
other_design_annotation = "OtherDesign"
control_annotation = "Control"

case_control_projects = c("PRJNA1012708", "PRJNA227074", "PRJNA266572", "PRJNA285798", "PRJNA285953", "PRJNA320361",
                          "PRJNA352062", "PRJNA369684", "PRJNA422124", "PRJNA447934", "PRJNA450428", "PRJNA467446",
                          "PRJNA470512", "PRJNA485140", "PRJNA539863", "PRJNA544126", "PRJNA551288", "PRJNA552599",
                          "PRJNA562638", "PRJNA597049", "PRJNA600939", "PRJNA638819", "PRJNA647880", "PRJNA656180",
                          "PRJNA723500", "PRJNA768419", "PRJNA789591", "PRJNA797845", "PRJNA882083", "PRJNA918345")

longitudinal_projects = c("PRJNA1012708", "PRJNA266572", "PRJNA285953", "PRJNA352062", "PRJNA369684", "PRJNA400331",
                          "PRJNA422124", "PRJNA467446", "PRJNA470512", "PRJNA482564", "PRJNA552599", "PRJNA560793",
                          "PRJNA647880", "PRJNA692462", "PRJNA723500", "PRJNA797845", "PRJNA862866", "PRJNA882083")

other_design_projects = c("PRJNA400331", "PRJNA482564", "PRJNA494963", "PRJNA560793", "PRJNA638819", "PRJNA683803",
                          "PRJNA692462", "PRJNA727526", "PRJNA798677", "PRJNA862866", "PRJNA911908", "PRJNA944738")

comorbidity_projects = c("PRJNA470512", "PRJNA683803", "PRJNA727526", "PRJNA798677", "PRJNA918345")
treatment_projects = c("PRJNA560793", "PRJNA683803")
vaccination_projects = c("PRJNA482564", "PRJNA692462")

#### technical reps ----
# BioProject, id_col
technical_repeats_projects = list(c(BioProject="PRJNA227074", Column="GEO_Accession..exp."),
                                  c(BioProject="PRJNA352062", Column="Sample_code"),
                                  c(BioProject="PRJNA551288", Column="Library_Name"),
                                  c(BioProject="PRJNA560793", Column="GEO_Accession..exp."),
                                  c(BioProject="PRJNA638819", Column="source_name"),
                                  c(BioProject="PRJNA656180", Column="GEO_Accession..exp."),
                                  c(BioProject="PRJNA485140", Column="human_donor"),
                                  c(BioProject="PRJNA600939", Column="GEO_Accession..exp."),
                                  c(BioProject="PRJNA862866", Column="GEO_Accession..exp."),
                                  c(BioProject="PRJNA735648", Column="GEO_Accession..exp."),
                                  c(BioProject="PRJNA390289", Column="TechRepID"),
                                  c(BioProject="PRJNA776746", Column="Donor"),
                                  c(BioProject="PRJNA649006", Column="donor_id"),
                                  c(BioProject="PRJNA882083", Column="Sample.Name"),
                                  c(BioProject="PRJNA918345", Column="SubjectCode"),
                                  c(BioProject="PRJNA285798", Column="name")
                                  ) %>%
  do.call(rbind, .) %>% as.data.frame()

#### paired samples ----
# BioProject, id_col
paired_samples_projects = list(c(BioProject="PRJNA227074", Column="ID"),
                               c(BioProject="PRJNA320361", Column="Individual"),
                               c(BioProject="PRJNA450428", Column="SubjectID"),
                               c(BioProject="PRJNA551288", Column="Library_Name"), # same as technical reps.
                               c(BioProject="PRJNA560793", Column="Subject"),
                               c(BioProject="PRJNA638819", Column="subject_id"),
                               c(BioProject="PRJNA723500", Column="volunteer_code"),
                               c(BioProject="PRJNA447934", Column="donor_id"),
                               c(BioProject="PRJNA485140", Column="human_donor"),
                               c(BioProject="PRJNA544126", Column="Donor"),
                               c(BioProject="PRJNA400331", Column="SubjectID"),
                               c(BioProject="PRJNA390289", Column="PatientID"),
                               c(BioProject="PRJNA695511", Column="Subject_ID"),
                               c(BioProject="PRJNA649006", Column="donor_id"),
                               c(BioProject="PRJNA736947", Column="sampleID"),
                               c(BioProject="PRJNA776746", Column="Donor"), # same as technical reps.
                               c(BioProject="PRJNA918345", Column="SubjectCode")  # same as technical reps.? not clear
                               ) %>% 
  do.call(rbind, .) %>% as.data.frame()

#### multiple tissues datasets ----
# see if need to split by that - YES
multiple_tissue_projects = list(c(BioProject="PRJNA450428", Column="source_name"),
                                c(BioProject="PRJNA552599", Column="source_name"),
                                c(BioProject="PRJNA562638", Column="Host"),
                                c(BioProject="PRJNA447934", Column="cell_type"),
                                c(BioProject="PRJNA949617", Column="cell_type"))%>%
  do.call(rbind, .) %>% as.data.frame()
# consider PRJNA597049?

#### different read length datasets ----
different_readLength_projects = list(c(BioProject="PRJNA400331", Column="Read.Leangth"),
                                    c(BioProject="PRJNA485140", Column="Read.Leangth"),
                                    c(BioProject="PRJNA683803", Column="Read.Leangth"),
                                    c(BioProject="PRJNA768419", Column="Read.Leangth"),
                                    c(BioProject="PRJNA944738", Column="Read.Leangth"))%>%
  do.call(rbind, .) %>% as.data.frame()


### create directories ------------------------------------------------------
#create directories, if do not exist
dir.create(out_boxplots, recursive = T)
dir.create(out_scatterplots, recursive = T)
dir.create(out_correlations, recursive = T)
dir.create(out_barplots, recursive = T)
dir.create(out_stats, recursive = T)
dir.create(out_tables, recursive = T)

# Annotations ----------------------------------------------
### get GEO annotations for data 1 ------
# PRJNA691474
gse <- getGEO(filename="Metadata/GSE164643_family.soft")
pertusis_data<-(GSMList(gse))
result_list <- lapply(pertusis_data, function(obj) c(obj@header$title, obj@header$geo_accession))
result_df <- do.call(rbind, result_list)
result_df<-as.data.frame(result_df)
colnames(result_df)<-c("states","geo_accession")
result_df = result_df %>%
  mutate(states= str_remove(states, pattern = ".*: "))

### get GEO annotations for data 2 ------
# PRJNA739691
gse <- getGEO(filename="Metadata/GSE178583_family.soft")
rabies_data<-(GSMList(gse))
result_list <- lapply(rabies_data, function(obj) c(obj@header$title, obj@header$geo_accession, obj@header$treatment_protocol_ch1))
result_df2 <- do.call(rbind, result_list)
result_df2<-as.data.frame(result_df2)
colnames(result_df2)<-c("title","geo_accession", "treatment_protocol_ch1")
result_df2$title <- sub(".*:\\s*", "", result_df2$title)
result_df2 = result_df2 %>%
  separate(title, into = c("Infection", "Time2", "V1", "V2", "Rep"), sep =" ") %>%
  select("Infection", "Time2","geo_accession")

### get GEO annotations for data 3 ------
# PRJNA285798 (repeats)
gse <- getGEO(filename="Metadata/GSE69529_family.soft")
diarrhae_data<-(GSMList(gse))
result_list <- lapply(diarrhae_data, function(obj) c(obj@header$characteristics_ch1, obj@header$title, obj@header$geo_accession))
result_df3 <- do.call(rbind, result_list)
result_df3<-as.data.frame(result_df3)
result_df3 <- rbind(result_df3 %>%
                       filter(grepl(x = V1, pattern = "age")) %>%
                       rename("name" = "V8", "geo_accession" = "V9") %>%
                       select("name", "geo_accession"),
                     result_df3 %>%
                       filter(grepl(x = V1, pattern = "group")) %>%
                       rename("name" = "V6", "geo_accession" = "V7") %>%
                       select("name", "geo_accession")) %>%
  mutate(name = str_remove(name, ", repeat"))

### get GEO annotations for data 4 ------
# PRJNA470512
gse <- getGEO(filename="Metadata/GSE114192_family.soft")
TB_DM<-(GSMList(gse))
result_list <- lapply(TB_DM, function(obj) c(obj@header$characteristics_ch1[2]))
result_df4 <- do.call(rbind, result_list)
result_df4<-as.data.frame(result_df4)
result_df4$geo_accession<-names((GSMList(gse)))
colnames(result_df4)<-c("states_tb_db","geo_accession")
result_df4$states_tb_db <- sub(".*:\\s*", "", result_df4$states_tb_db)

### get GEO annotations for data 5 ------
#	PRJNA369684
gse <- getGEO(filename="Metadata/GSE94438_family.soft")
tb2_data<-(GSMList(gse))
result_list <- lapply(tb2_data, function(obj) c(obj@header$characteristics_ch1, obj@header$geo_accession))
result_df5 <- do.call(rbind, result_list)
result_df5<-as.data.frame(result_df5)
colnames(result_df5)<-c("tissue", "code", "subjectid", "site", "age", "gender", "group", "time_from_exposure_months", "time_to_tb_months", "geo_accession")
result_df5 = result_df5 %>%
  mutate(across(everything(), ~ str_remove(pattern = ".*: ", .x))) %>%
  select("group", "time_from_exposure_months",  "geo_accession")

### get GEO annotations for data 6 ------
# PRJNA695511
gse <- getGEO(filename="Metadata/GSE165708_family.soft")
tb_hiv_data<-(GSMList(gse))
result_list <- lapply(tb_hiv_data, function(obj) c(obj@header$description, obj@header$geo_accession))
result_df6 <- do.call(rbind, result_list)
result_df6<-as.data.frame(result_df6)
colnames(result_df6)<-c("States", "v1", "v2", "geo_accession")
result_df6 = result_df6 %>%
  separate(States, into = c("State", "tissue", "test"), sep = ", ") %>%
  mutate(Subject_ID = State,
         State = str_remove(State, pattern = " sample.*"),
         test = str_remove(test, pattern = "Mycobacterium tuberculosis ")) %>%
  select("State", "Subject_ID", "test",  "geo_accession")

### get GEO annotations for data 7 ------
# PRJNA560793
gse <- getGEO(filename="Metadata/GSE135965_family.soft")
leishmania_data<-(GSMList(gse))
result_list <- lapply(leishmania_data, function(obj) c(obj@header$title, obj@header$geo_accession))
result_df7 <- do.call(rbind, result_list)
result_df7<-as.data.frame(result_df7)
colnames(result_df7)<-c("States", "geo_accession")
result_df7 = result_df7 %>%
  separate(States, into = c("Subject", "state"), sep = "CD4") %>%
  select("Subject", "state",  "geo_accession")

### get GEO annotations for data 8 ------
#	PRJNA390289
gse <- getGEO(filename="Metadata/GSE99992_family.soft")
Chikungunya_data<-(GSMList(gse))
result_list <- lapply(Chikungunya_data, function(obj) c(obj@header$title, obj@header$geo_accession))
result_df8 <- do.call(rbind, result_list)
result_df8<-as.data.frame(result_df8)
colnames(result_df8)<-c("States", "geo_accession")
result_df8 = result_df8 %>%
  separate(States, into = c("PatientID", "sample_state", "TechRepID"), sep = ", ") %>%
  mutate(TechRepID = paste0(PatientID, "_", sample_state) %>% str_replace_all(pattern = " ", replacement = "_")) %>%
  select("PatientID", "sample_state", "TechRepID",  "geo_accession")

### get GEO annotations for data 9 ------
#	PRJNA736947
gse <- getGEO(filename="Metadata/GSE177040_family.soft")
monocytes_stimulation_data<-(GSMList(gse))
result_list <- lapply(monocytes_stimulation_data, function(obj) c(obj@header$title, obj@header$geo_accession))
result_df9 <- do.call(rbind, result_list)
result_df9<-as.data.frame(result_df9)
colnames(result_df9)<-c("d", "geo_accession")
result_df9 = result_df9 %>%
  separate(col = d, into = c("sampleID", "stimulation_type", "time"), remove = T)

### get GEO annotations for data 10 ------
# PRJNA918345
gse <- getGEO(filename="Metadata/GSE222129_family.soft")
Rhinovirus_data<-(GSMList(gse))
result_list <- lapply(Rhinovirus_data, function(obj) c(obj@header$title, obj@header$geo_accession))
result_df10 <- do.call(rbind, result_list)
result_df10<-as.data.frame(result_df10)
colnames(result_df10)<-c("title", "geo_accession")
result_df10 = result_df10 %>%
  separate(title, into = c("Strain", "state", "v1", "lane"), sep = "_") %>%
  mutate(SubjectCode = paste0(state, "_", v1)) %>%
  select("SubjectCode", "Strain",  "geo_accession")

### get GEO annotations for data 11 ------
# PRJNA467446
gse <- getGEO(filename="Metadata/GSE114180_family.soft")
candida_Data<-(GSMList(gse))
result_list <- lapply(candida_Data, function(obj) c(obj@header$title, obj@header$characteristics_ch1, obj@header$geo_accession))
result_df11 <- do.call(rbind, result_list[])
result_df11<-as.data.frame(result_df11)
result_df11 = rbind(result_df11 %>%
                      filter(grepl(x = V3, pattern = "time")) %>%
                      select(V2,V4) %>%
                      rename("strain" = V2, "geo_accession" = V4),
                    result_df11 %>%
                      filter(!grepl(x = V3, pattern = "time")) %>%
                      select(V3) %>%
                      rename("geo_accession" = V3) %>%
                      mutate(strain = "Uninfected")) %>%
  mutate(strain = str_remove(strain, pattern = ": .*"))

# add to full annotations
sample_info<-read.csv("combined_data_with_all_sra_samples_and_info.fixedAnnotations.csv") 

# add annotation for new samples
sample_info<-sample_info%>%
  merge(result_df, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df2, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T)  %>%
  merge(result_df3, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T)  %>%
  merge(result_df4, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df5, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df6, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df7, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df8, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df9, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T) %>%
  merge(result_df10, by.x="Library.Name", by.y="geo_accession", all.x = T) %>%
  merge(result_df11, by.x="GEO_Accession..exp.", by.y="geo_accession", all.x = T)


# fix annotations
sample_info<-sample_info%>%
  mutate(OrganismType=forcats::fct_collapse(OrganismType, Parasites= c("parasite", "parasites")))%>%
  mutate(OrganismType=forcats::fct_collapse(OrganismType, Bacterial= c("bacterial", "bacterial ", "bacteria")))%>%
  mutate(OrganismType=forcats::fct_collapse(OrganismType, All= c("ALL","all","all (not fungal)")))%>%
  mutate(OrganismType=forcats::fct_collapse(OrganismType, Virus= c("Virus","virus")))%>%
  mutate(OrganismType=forcats::fct_collapse(OrganismType, `Bacteria+Virus`= c("bacteria and virus","Virus and bacteria")))%>%
  mutate(OrganismType=forcats::fct_collapse(OrganismType, Unkown= c("Unkown","Unkown. its unkown what causes teh disease, its 100% unkown and thus considered a very good disease to research due to theneed of information.")))


# Fix Classifications -----------------------------------------------------

### Keep original classification - QC ------
manually_annotated_disease_status<-sample_info%>%
  mutate(Disease_Status_original_annotation = case_when(
    BioProject == "PRJNA691474" ~ states, #
    BioProject == "PRJNA470512" ~ states_tb_db, #
    BioProject == "PRJNA700069" ~ culture,
    BioProject == "PRJNA918345" ~ treatment,
    BioProject == "PRJNA649006" ~ paste(outcome, Timepoint, sep = " | "),
    BioProject == "PRJNA656180" ~ pathogen,
    BioProject == "PRJNA352062" ~ paste(disease_state, Time, sep = " | "),
    BioProject == "PRJNA369684" ~ paste(Group, time_from_exposure_months, sep = " | "),
    BioProject == "PRJNA400331" ~ paste(Group, Timepoint, qft, sep = " | "),
    BioProject == "PRJNA422124" ~ Group,
    BioProject == "PRJNA447934" ~ paste(infection,cell_type, sep = " | "),
    BioProject == "PRJNA695511" ~ paste(State, test, sep = " | "),
    BioProject == "PRJNA750782" ~ paste(disease_state, tb_status, Timepoint, sep = " | "),
    BioProject == "PRJNA797845" ~ distinct_stages_of_tb,
    BioProject == "PRJNA692462" ~ paste(carriage, Timepoint, sep = " | "),
    BioProject == "PRJNA776746" ~ paste(condition, infection, Time_point, sep = " | "),
    BioProject == "PRJNA723500" ~ source_name,
    BioProject == "PRJNA227074" ~ pf_infection_status,
    BioProject == "PRJNA551288" ~ infected_with.healthy_control,
    BioProject == "PRJNA660611" ~ status,
    BioProject == "PRJNA1012708" ~ disease_state,
    BioProject == "PRJNA638819" ~ Time_point,
    BioProject == "PRJNA768419" ~ paste(disease_state, collection_site, sep = "|"),
    BioProject == "PRJNA552599" ~ paste(patient_group, cell_type, sepsis_stage, sep = " | "),
    BioProject == "PRJNA647880" ~ status,
    BioProject == "PRJNA285798" ~ paste(organisms, Group, sep = " | "),
    BioProject == "PRJNA600939" ~ status,
    BioProject == "PRJNA320361" ~ infection,
    BioProject == "PRJNA735648" ~ paste(PHENOTYPE, Time_point, daysreltofirsttimepoin, daysreltofirsttimepoint, sep = " | "),
    BioProject == "PRJNA467446" ~ paste(Time_point, strain, sep = " | "),
    BioProject == "PRJNA485140" ~ paste(viral_infection,time_post_infection_viral_infection, viral_strain, sep = " | "),
    BioProject == "PRJNA544126" ~ treatment, #Human iris pigment epithelial cell (HIPE)
    BioProject == "PRJNA562638" ~ paste(infection, Host, sep = " | "),
    BioProject == "PRJNA266572" ~ paste(disease_state, Time, sep = " | "),
    BioProject == "PRJNA736947" ~ paste(stimulation, Time_point, sep = " | "),
    BioProject == "PRJNA597049" ~ Schisto_status,
    BioProject == "PRJNA285953" ~ paste(Group, history_of_brucellosis, sep = " | "),
    BioProject == "PRJNA539863" ~ source_name,
    BioProject == "PRJNA686397" ~ paste(disease_state, Time_point, sep = " | "),
    BioProject == "PRJNA911908" ~ diagnosis,
    BioProject == "PRJNA789591" ~ Group,
    BioProject == "PRJNA390289" ~ Timepoint,
    BioProject == "PRJNA450428" ~ paste(Visit, source_name, sep = " | "),
    BioProject == "PRJNA482564" ~ paste(Timepoint, vaccine, carriage_status, "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA494963" ~ clinical_form,
    BioProject == "PRJNA560793" ~ status,
    BioProject == "PRJNA662344" ~ paste(Group, days_from_att < 0, sep = " | "),
    BioProject == "PRJNA683803" ~ paste(treatment, outcome, grouping, "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA727526" ~ paste(ETHNICITY, gender, hbv_status, treatment, viral_load_group, "CHECK FOR MORE", sep = " | "), 
    BioProject == "PRJNA739691" ~ paste(Infection, Time2, sep = " | "),
    BioProject == "PRJNA862866" ~ paste(disease_state, study_participant, Time, sep = " | "),
    BioProject == "PRJNA798677" ~ paste(Timepoint, disease_category, field_site, "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA882083" ~ paste(Time, hcv.control, sep = " | "),
    BioProject == "PRJNA944738" ~ paste(sex, disease_state, time_point_in_exercise, "CHECK FOR MORE", sep = " | "), 
    BioProject == "PRJNA949617" ~ paste(treatment, cell_type, sep = " | "),
    TRUE ~ "NOT RUN"  # Keep the original value if no conditions are met
  ),
  ### add colnames -----
  Disease_Status_original_annotation_columns = case_when(
    BioProject == "PRJNA691474" ~ "states", #
    BioProject == "PRJNA470512" ~ "states_tb_db", #
    BioProject == "PRJNA700069" ~ "culture",
    BioProject == "PRJNA918345" ~ "treatment",
    BioProject == "PRJNA649006" ~ paste("outcome", "Timepoint", sep = " | "),
    BioProject == "PRJNA656180" ~ "pathogen",
    BioProject == "PRJNA352062" ~ paste("disease_state", "Time", sep = " | "),
    BioProject == "PRJNA369684" ~ paste("Group", "time_from_exposure_months", sep = " | "),
    BioProject == "PRJNA400331" ~ paste("Group", "Timepoint", "qft", sep = " | "),
    BioProject == "PRJNA422124" ~ "Group",
    BioProject == "PRJNA447934" ~ paste("infection","cell_type", sep = " | "),
    BioProject == "PRJNA695511" ~ paste("State", "test", sep = " | "),
    BioProject == "PRJNA750782" ~ paste("disease_state", "tb_status", "Timepoint", sep = " | "),
    BioProject == "PRJNA797845" ~ "distinct_stages_of_tb",
    BioProject == "PRJNA692462" ~ paste("carriage", "Timepoint", sep = " | "),
    BioProject == "PRJNA776746" ~ paste("condition", "infection", "Time_point", sep = " | "),
    BioProject == "PRJNA723500" ~ "source_name",
    BioProject == "PRJNA227074" ~ "pf_infection_status",
    BioProject == "PRJNA551288" ~ "infected_with.healthy_control",
    BioProject == "PRJNA660611" ~ "status",
    BioProject == "PRJNA1012708" ~ "disease_state",
    BioProject == "PRJNA638819" ~ "Time_point",
    BioProject == "PRJNA768419" ~ paste("disease_state", "collection_site", sep = "|"),
    BioProject == "PRJNA552599" ~ paste("patient_group", "cell_type", "sepsis_stage", sep = " | "),
    BioProject == "PRJNA647880" ~ "status",
    BioProject == "PRJNA285798" ~ paste("organisms", "Group", sep = " | "),
    BioProject == "PRJNA600939" ~ "status",
    BioProject == "PRJNA320361" ~ "infection",
    BioProject == "PRJNA735648" ~ paste("PHENOTYPE", "Time_point", "daysreltofirsttimepoin", "daysreltofirsttimepoint", sep = " | "),
    BioProject == "PRJNA467446" ~ paste("Time_point", "strain", sep = " | "),
    BioProject == "PRJNA485140" ~ paste("viral_infection","time_post_infection_viral_infection", "viral_strain", sep = " | "),
    BioProject == "PRJNA544126" ~ "treatment", #Human iris pigment epithelial cell (HIPE)
    BioProject == "PRJNA562638" ~ paste("infection", "Host", sep = " | "),
    BioProject == "PRJNA266572" ~ paste("disease_state", "Time", sep = " | "),
    BioProject == "PRJNA736947" ~ paste("stimulation", "Time_point", sep = " | "),
    BioProject == "PRJNA597049" ~ "Schisto_status",
    BioProject == "PRJNA285953" ~ paste("Group", "history_of_brucellosis", sep = " | "),
    BioProject == "PRJNA539863" ~ "source_name",
    BioProject == "PRJNA686397" ~ paste("disease_state", "Time_point", sep = " | "),
    BioProject == "PRJNA911908" ~ "diagnosis",
    BioProject == "PRJNA789591" ~ "Group",
    BioProject == "PRJNA390289" ~ "Timepoint",
    BioProject == "PRJNA450428" ~ paste("Visit", "source_name", sep = " | "),
    BioProject == "PRJNA482564" ~ paste("Timepoint", "vaccine", "carriage_status", "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA494963" ~ "clinical_form",
    BioProject == "PRJNA508744" ~ paste("Donor", "Incubation_time", "treatment_condition", "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA560793" ~ "status",
    BioProject == "PRJNA662344" ~ paste("Group", "days_from_att < 0", sep = " | "),
    BioProject == "PRJNA683803" ~ paste("treatment", "outcome", "grouping", "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA727526" ~ paste("ETHNICITY", "gender", "hbv_status", "treatment", "viral_load_group", "CHECK FOR MORE", sep = " | "), 
    BioProject == "PRJNA739691" ~ paste("Infection", "Time2", sep = " | "),
    BioProject == "PRJNA862866" ~ paste("disease_state", "study_participant", "Time", sep = " | "),
    BioProject == "PRJNA798677" ~ paste("Timepoint", "disease_category", "field_site", "CHECK FOR MORE", sep = " | "),
    BioProject == "PRJNA882083" ~ paste("Time", "hcv.control", sep = " | "),
    BioProject == "PRJNA944738" ~ paste("sex", "disease_state", "time_point_in_exercise", "CHECK FOR MORE", sep = " | "), 
    BioProject == "PRJNA949617" ~ paste("treatment", "cell_type", sep = " | "),
    TRUE ~ "NOT RUN"  # Keep the original value if no conditions are met
  )) %>%
  
  
  ### add classifications ----
  mutate(Disease_Status_manual_annotation = case_when(
    BioProject == "PRJNA691474" & !grepl("unstimulated", states, ignore.case = TRUE)  ~ "Pertussis", #
    BioProject == "PRJNA691474" & grepl("unstimulated", states, ignore.case = TRUE)  ~ control_annotation, #
    BioProject == "PRJNA470512" & states_tb_db == "TB_DM" ~ "TB", #
    BioProject == "PRJNA470512" & states_tb_db == "TB_only" ~ "TB",
    BioProject == "PRJNA470512" & states_tb_db == "TB_IH" ~ "TB", #
    BioProject == "PRJNA470512" & states_tb_db == "DM_only" ~ control_annotation,
    BioProject == "PRJNA470512" & states_tb_db == "IH" ~ control_annotation, #
    BioProject == "PRJNA470512" & states_tb_db == "Healthy_Control" ~ control_annotation,
    BioProject == "PRJNA700069" & culture == "rhinovirus stimulated" ~ "RV", #
    BioProject == "PRJNA700069" & culture == "Unstimulated" ~ control_annotation,
    BioProject == "PRJNA918345" & treatment  == "Rhinovirus B48 infection" ~ "RV B48", #
    BioProject == "PRJNA918345" & treatment  == "Rhinovirus C15 infection" ~ "RV C15",
    BioProject == "PRJNA918345" & treatment  == "none" ~ control_annotation,
    BioProject == "PRJNA649006" & outcome  == "Infected" & Timepoint == "D0" ~ control_annotation, #
    BioProject == "PRJNA649006" & outcome  == "Infected" & Timepoint == "D3" ~ "RSV", #
    BioProject == "PRJNA656180" & pathogen  == "negative"  ~ control_annotation,
    BioProject == "PRJNA656180" & pathogen  != "negative"  ~ pathogen,
    
    
    BioProject == "PRJNA352062" & disease_state  == "TB Subjects" & Time   == "DX"  ~ "TB", #
    BioProject == "PRJNA352062" & disease_state  == "Healthy Controls" & Time   == "DX"  ~ control_annotation,
    BioProject == "PRJNA422124" & Group  == "Active_TB"  ~ "TB",
    BioProject == "PRJNA422124" & Group  == "Control"  ~ control_annotation,
    BioProject == "PRJNA447934" & infection    == "Mtb"  ~ "TB", #
    BioProject == "PRJNA447934" & infection    == "uninf"  ~ control_annotation,
    BioProject == "PRJNA750782" &  tb_status=="TB" &  Timepoint=="0" ~ "TB", #
    BioProject == "PRJNA750782" &  tb_status=="non TB" &  Timepoint=="0" ~ control_annotation,
    BioProject == "PRJNA797845" & distinct_stages_of_tb     == "Active TB_baseline"  ~ "TB", #
    BioProject == "PRJNA797845" & distinct_stages_of_tb     == "HH controls_baseline"  ~ control_annotation,
    BioProject == "PRJNA692462" & carriage  == "POS" & Timepoint    == "baseline" ~ "S. pneumoniae", #
    BioProject == "PRJNA692462" & carriage  == "NEG" & Timepoint    == "baseline" ~ control_annotation,
    BioProject == "PRJNA776746" & condition == "NHB" & infection == "influenza A virus" & Time_point == "18 hr"  ~ "Influenza", #
    BioProject == "PRJNA776746" & condition == "NHB" & infection == "control" & Time_point == "18 hr" ~ control_annotation,
    BioProject == "PRJNA723500" & source_name == "whole blood before first inf (baseline)"  ~ control_annotation,
    BioProject == "PRJNA723500" & source_name == "whole blood after first inf (diagnosis)"  ~ "Malaria",
    BioProject == "PRJNA227074" & pf_infection_status == "U"  ~ control_annotation,
    BioProject == "PRJNA227074" & pf_infection_status == "I"  ~ "Malaria",
    BioProject == "PRJNA551288" & infected_with.healthy_control == "Enterovirus"  ~ "Enterovirus",
    BioProject == "PRJNA551288" & infected_with.healthy_control == control_annotation  ~ control_annotation,
    BioProject == "PRJNA660611" & status == "DNA virus infected patient (Adenovirus, Cytomegalovirus, Ebstein-Barr virus, Herpes Simplex virus)"  ~ "DNA virus",
    BioProject == "PRJNA660611" & status == "other respiratory RNA virus infected patient (Parainfluenza virus and Respiratory Syncytial virus)"  ~ "Parainfluenza,RSV",
    BioProject == "PRJNA660611" & status == "Entero/Rhinovirus infected patient"  ~ "Entero/Rhinovirus",
    BioProject == "PRJNA660611" & status == "Metapneumovirus infected patient"  ~ "Metapneumovirus",
    BioProject == "PRJNA660611" & status == "Influenza virus infected patient"  ~ "Influenza",
    BioProject == "PRJNA660611" & status == "Dengue virus infected patient"  ~ "Dengue",
    BioProject == "PRJNA660611" & status == "healthy_ctrl"  ~ control_annotation,
    BioProject == "PRJNA1012708" & disease_state  == "Healthy control"  ~ control_annotation,
    BioProject == "PRJNA1012708" & disease_state  == "Herpes Zoster"  ~ "Herpes Zoster",
    BioProject == "PRJNA638819" & Time_point   == "day 0"  ~ control_annotation,
    BioProject == "PRJNA638819" & Time_point   == "day 8"  ~ "Dengue",
    BioProject == "PRJNA768419" & disease_state == "healthy"  ~ control_annotation,
    BioProject == "PRJNA768419" & disease_state == "sepsis"  ~ "Sepsis",
    BioProject == "PRJNA552599" & patient_group == "Healthy"  ~ control_annotation,
    BioProject == "PRJNA552599" & patient_group == "Sepsis" & sepsis_stage == "Early" ~ "Sepsis",
    BioProject == "PRJNA647880" & status == "Hlty"  ~ control_annotation,
    BioProject == "PRJNA647880" & status == "Seps_P"  ~ "Sepsis",
    BioProject == "PRJNA285798" & organisms == "Negative" & Group == "Control" ~ control_annotation,
    BioProject == "PRJNA285798" & organisms == "Adenovirus" & Group == "Diarrhea" ~ "Adenovirus",
    BioProject == "PRJNA285798" & organisms == "DAEC" & Group == "Diarrhea" ~ "Diffusely-adherent E.coli",
    BioProject == "PRJNA285798" & organisms == "EAEC" & Group == "Diarrhea" ~ "Enteroaggregative E.coli",
    BioProject == "PRJNA285798" & organisms == "EPEC" & Group == "Diarrhea" ~ "Enteropathogenic E.coli",
    BioProject == "PRJNA285798" & organisms == "Norovirus" & Group == "Diarrhea" ~ "Norovirus",
    BioProject == "PRJNA285798" & organisms == "Rotavirus" & Group == "Diarrhea" ~ "Rotavirus",
    BioProject == "PRJNA285798" & organisms == "Salmonella" & Group == "Diarrhea" ~ "Salmonella",
    BioProject == "PRJNA285798" & organisms == "Shigella" & Group == "Diarrhea" ~ "Shigella",
    BioProject == "PRJNA600939" & status  == "Healthy"  ~ control_annotation,
    BioProject == "PRJNA600939" & status  == "Ebola virus disease survivors"  ~ "Ebola",
    BioProject == "PRJNA320361" & infection == "Non-infected"  ~ control_annotation,
    BioProject == "PRJNA320361" & infection  == "Salmonella"  ~ "Salmonella",
    BioProject == "PRJNA320361" & infection  == "Listeria"  ~ "Listeria",
    BioProject == "PRJNA735648" & PHENOTYPE == "Healthy"  ~ control_annotation,
    BioProject == "PRJNA735648" & PHENOTYPE == "Candidemia"  & (daysreltofirsttimepoin == 0 | daysreltofirsttimepoint == 0 | (is.na(daysreltofirsttimepoin) & is.na(daysreltofirsttimepoint))) ~ "Candida",
    BioProject == "PRJNA735648" & PHENOTYPE == "Bacterial"  ~ "Bacterial",
    BioProject == "PRJNA735648" & PHENOTYPE == "Viral"  ~ "Viral",
    BioProject == "PRJNA735648" & PHENOTYPE == "Mixed Candida/bacterial" & (daysreltofirsttimepoin == 0 | daysreltofirsttimepoint == 0 | (is.na(daysreltofirsttimepoin) & is.na(daysreltofirsttimepoint)))  ~ "Mixed Candida/Bacterial",
    BioProject == "PRJNA467446" & strain == "Uninfected" ~ control_annotation,
    BioProject == "PRJNA467446" & strain != "Uninfected" & Time_point == "0 minutes"  ~ "Candida",
    BioProject == "PRJNA485140" & viral_infection  == "mock infected" & time_post_infection_viral_infection == "24h" ~ control_annotation,
    BioProject == "PRJNA485140" & viral_infection  == "ZIKV-infected" & time_post_infection_viral_infection == "24h" & viral_strain != "ZIKV-SD001" ~ "Zika",
    BioProject == "PRJNA544126" & treatment == "Mock infection"  ~ control_annotation,
    BioProject == "PRJNA544126" & treatment == "Zika virus infection"  ~ "Zika", #Human iris pigment epithelial cell (HIPE)
    BioProject == "PRJNA562638" & infection == "Mock"  ~ control_annotation,
    BioProject == "PRJNA562638" & infection == "Vibrio vulnificus (MO6-24) wild type" ~ "V. vulnificus",
    BioProject == "PRJNA266572" & disease_state  == "Lyme disease" & Time == "acute Lyme pre-treatment (V1)"  ~ "Lyme",
    BioProject == "PRJNA266572" & disease_state  == "Healthy controls"  ~ control_annotation,
    BioProject == "PRJNA736947" & stimulation == "Control"  & Time_point == "6h" ~ control_annotation,
    BioProject == "PRJNA736947" & stimulation == "A. fumigatus" & Time_point == "6h"  ~ "A. fumigatus (Fungal)",
    BioProject == "PRJNA736947" & stimulation == "N. meningitidis" & Time_point == "6h" ~ "N. meningitidis (Gram -)",
    BioProject == "PRJNA736947" & stimulation == "S. aureus"  & Time_point == "6h" ~ "S. aureus (Gram +)",
    BioProject == "PRJNA597049" & Schisto_status == "none"  ~ control_annotation,
    BioProject == "PRJNA597049" & Schisto_status != "none"  ~ Schisto_status,
    BioProject == "PRJNA285953" & Group == "Leishmaniasis"   ~ "Leishmaniasis",
    BioProject == "PRJNA285953" & Group == "Control"   ~ control_annotation,
    BioProject == "PRJNA285953" & Group == "Brucellosis" & history_of_brucellosis == "Acute brucellosis"   ~ "Brucellosis",
    BioProject == "PRJNA539863" & source_name == "Healthy control_PBMCs"  ~ control_annotation,
    BioProject == "PRJNA539863" & source_name == "Q Fever Fatigue Syndrom (QFS)_PBMCs"  ~ "Q Fever",
    BioProject == "PRJNA686397" & Time_point == "HHC"  ~ control_annotation,
    BioProject == "PRJNA686397" & Time_point  == "Second"  ~ "Leprosy",
    BioProject == "PRJNA789591" & Group == "non-chagasic control"  ~ control_annotation,
    BioProject == "PRJNA789591" & Group == "chronic chagasic cardiomyopathy"  ~ "chronic chagasic cardiomyopathy",
    BioProject == "PRJNA662344" & Group == "Control" ~ control_annotation,
    BioProject == "PRJNA662344" & Group == "PTB" & days_from_att <0 ~ "TB",
    BioProject == "PRJNA450428" & Visit == "AV" ~ "AV Bronchitis",
    BioProject == "PRJNA450428" & Visit == "CV" ~ control_annotation,
    BioProject == "PRJNA739691" & Infection == "Tha-infected" & Time2 == "24hr" ~ "Rabies",
    BioProject == "PRJNA739691" & Infection == "non-infected" & Time2 == "24hr" ~ control_annotation,
    BioProject == "PRJNA882083" & Time == "acute" & hcv.control == control_annotation ~ control_annotation,
    BioProject == "PRJNA882083" & Time == "acute" & hcv.control == "Hepatitis C" ~ "Hepatitis C",
    BioProject == "PRJNA949617" & cell_type != "CD14+ primary monocytes" & treatment == "ARPE19-VZV ORF23 cell lysate batch 3" ~ "Chickenpox",
    BioProject == "PRJNA949617" & cell_type != "CD14+ primary monocytes" & treatment == "unstimulated" ~ control_annotation,
    BioProject == "PRJNA695511" & State == "Healthy control" & test == "non-challenged library" ~ control_annotation,
    BioProject == "PRJNA695511" & State == "Healthy control" & test == "challenged library" ~ "TB",
    BioProject == "PRJNA390289" & Timepoint == "1-2d post symptom onset" ~ "Chikungunya",
    BioProject == "PRJNA390289" & Timepoint == "15-17d post symptom onset" ~ control_annotation,
    TRUE ~ other_design_annotation ))%>%
  ### add special traits ----
  mutate(paired_sample_id = case_when(
    BioProject == "PRJNA227074" ~ ID %>% as.character(),
    BioProject == "PRJNA320361" ~ Individual %>% as.character(),
    BioProject == "PRJNA450428" ~ SubjectID %>% as.character(),
    BioProject == "PRJNA551288" ~ Library_Name %>% str_extract(pattern = "[A-Za-z0-9]+"),
    BioProject == "PRJNA560793" ~ Subject %>% as.character(),
    BioProject == "PRJNA638819" ~ subject_id %>% as.character(),
    BioProject == "PRJNA723500" ~ volunteer_code %>% as.character(),
    BioProject == "PRJNA447934" ~ donor_id %>% as.character(),
    BioProject == "PRJNA485140" ~ human_donor %>% as.character(),
    BioProject == "PRJNA544126" ~ Donor %>% as.character(),
    BioProject == "PRJNA400331" ~ SubjectID %>% as.character(),
    BioProject == "PRJNA390289" ~ PatientID %>% as.character(),
    BioProject == "PRJNA695511" ~ Subject_ID %>% as.character(),
    BioProject == "PRJNA649006" ~ donor_id %>% as.character(),
    BioProject == "PRJNA736947" ~ sampleID %>% as.character(),
    BioProject == "PRJNA776746" ~ Donor %>% as.character(),
    BioProject == "PRJNA918345" ~ SubjectCode %>% as.character(),
    .default = NA_character_  
  ),
  technical_repeats_id = case_when(
    BioProject %in% c("PRJNA227074","PRJNA560793","PRJNA656180",
                      "PRJNA600939","PRJNA862866", "PRJNA735648") ~ `GEO_Accession..exp.` %>% as.character(),
    BioProject == "PRJNA352062" ~ Sample_code %>% as.character(),
    BioProject == "PRJNA551288" ~ Library_Name %>% as.character(),
    BioProject == "PRJNA638819" ~ source_name %>% as.character(),
    BioProject == "PRJNA390289" ~ TechRepID %>% as.character(),
    BioProject == "PRJNA776746" ~ Donor %>% as.character(),
    BioProject == "PRJNA485140" ~ human_donor %>% as.character(),
    BioProject == "PRJNA649006" ~ donor_id %>% as.character(),
    BioProject == "PRJNA918345" ~ SubjectCode %>% as.character(),
    BioProject == "PRJNA882083" ~ Sample.Name %>% as.character(),
    BioProject == "PRJNA285798" ~ name %>% as.character(),
    .default = NA_character_  
  ),
  multiple_tissue_id = case_when(
    BioProject == "PRJNA450428" | BioProject == "PRJNA552599" ~ source_name,
    BioProject == "PRJNA562638" ~ Host,
    BioProject == "PRJNA447934" | BioProject == "PRJNA949617" ~ cell_type,
    .default = NA_character_  
  )) %>%
  
  select(Disease_Status_manual_annotation, Disease_Status_original_annotation, Disease_Status_original_annotation_columns,BioProject,Run,Tissue,OrganismType,DiseaseType,Disease, GEO_Accession..exp., paired_sample_id, technical_repeats_id, multiple_tissue_id)



# table(manually_annotated_disease_status$Disease_Status_manual_annotation)

manually_annotated_disease_status<-manually_annotated_disease_status%>%
  mutate(Hypothesis = case_when(
    BioProject == "PRJNA691474" ~ "Pertussis", #
    BioProject == "PRJNA470512" ~ "TB_DM_IH", #
    BioProject == "PRJNA700069" ~ "RV1", #
    BioProject == "PRJNA918345" ~ "RV2", #
    BioProject == "PRJNA649006" ~ "RSV", #
    BioProject == "PRJNA656180" ~ "RSV_RV", #
    BioProject == "PRJNA352062" ~ "TB1",
    BioProject == "PRJNA369684" ~ "TB2",
    BioProject == "PRJNA400331" ~ "TB3", #
    BioProject == "PRJNA422124" ~ "TB_Stages", #
    BioProject == "PRJNA447934" ~ "TB4", #
    BioProject == "PRJNA695511" ~ "TB5",
    BioProject == "PRJNA750782" ~ "TB6", #
    BioProject == "PRJNA797845" ~ "TB7", #
    BioProject == "PRJNA692462" ~ "S. pneumoniae", #
    BioProject == "PRJNA776746" ~ "Influenza A_PolyIC",
    ##not URI#
    BioProject == "PRJNA723500" ~ "Malaria",
    BioProject == "PRJNA227074" ~ "Malaria2",
    BioProject == "PRJNA551288" ~ "EV",
    BioProject == "PRJNA660611" ~ "Mixed0",
    BioProject == "PRJNA1012708" ~ "HZ",
    BioProject == "PRJNA638819" ~ "Dengue",
    BioProject == "PRJNA768419" ~ "Sepsis",
    BioProject == "PRJNA647880" ~ "Sepsis2",
    BioProject == "PRJNA552599" ~ "Sepsis3",
    BioProject == "PRJNA285798" ~ "Diarhea",
    BioProject == "PRJNA320361" ~ "Diarhea2",
    BioProject == "PRJNA600939" ~ "Ebola",
    BioProject == "PRJNA735648" ~ "Mixed",
    BioProject == "PRJNA467446" ~ "Candida",
    BioProject == "PRJNA485140" ~ "Zika",
    BioProject == "PRJNA544126" ~ "Zika2",
    BioProject == "PRJNA562638" ~ "V.F",
    BioProject == "PRJNA266572" ~ "Lyme",
    BioProject == "PRJNA736947" ~ "Mixed2",
    BioProject == "PRJNA597049" ~ "Shcisto",
    BioProject == "PRJNA285953" ~ "Leishmaniasis/Brucellosis",
    BioProject == "PRJNA539863" ~ "Q",
    BioProject == "PRJNA686397" ~ "Leprosy",
    BioProject == "PRJNA911908" ~ "Myostisis",
    BioProject == "PRJNA789591" ~ "Chagas",
    BioProject == "PRJNA450428" ~ "Viral Bronchiolitis",
    BioProject == "PRJNA662344" ~ "TB_Stages2",
    BioProject == "PRJNA739691" ~ "Rabies",
    BioProject == "PRJNA882083" ~ "Hepatitis",
    BioProject == "PRJNA390289" ~ "Chikungunya",
    BioProject == "PRJNA494963" ~ "Leprosy2",
    BioProject == "PRJNA683803" ~ "Immunodeficiency",
    BioProject == "PRJNA727526" ~ "Hepatitis2",
    BioProject == "PRJNA798677" ~ "TB_DM2",
    BioProject == "PRJNA944738" ~ "ME/CFS",
    BioProject == "PRJNA482564" ~ "Pneumonia",
    BioProject == "PRJNA560793" ~ "Visceral leishmaniasis",
    BioProject == "PRJNA862866" ~ "Malaria_Cycle",
    BioProject == "PRJNA949617" ~ "Chickenpox",
    TRUE ~ other_design_annotation 
  ))


# add dataset features ----------------------------------------------------


manually_annotated_disease_status<-manually_annotated_disease_status%>%
  mutate(Tissue=ifelse(Tissue=="Human airway epithelial cells","HAEC",Tissue))%>%
  mutate(Tissue=ifelse(Tissue=="Left ventricular wall of heart","LVW",Tissue))%>%
  mutate(Hypothesis=paste(Hypothesis, Tissue,sep=" ")) %>%
  filter(Disease_Status_original_annotation != "NOT RUN")

# add whether there are technical replicates, paired samples, study type
manually_annotated_disease_status <- manually_annotated_disease_status %>%
  # annotate according to constant BioProject annotations
  mutate(PairedData = if_else(BioProject %in% paired_samples_projects$BioProject,
                                      "Yes", 
                                      "No"),
         TechicalReplicates = if_else(BioProject %in% technical_repeats_projects$BioProject,
                                      "Yes", 
                                      "No"),
         MultipleTissuesData = if_else(BioProject %in% multiple_tissue_projects$BioProject,
                                      "Yes", 
                                      "No"),
         DifferentReadLengthData = if_else(BioProject %in% different_readLength_projects$BioProject,
                                           "Yes", 
                                           "No"),
         # annotate according to if other design or not (where other design means groups who were not annotated as case-control)
         StudyDesign = if_else(Disease_Status_manual_annotation == other_design_annotation,
                               other_design_annotation,
                               "Case-Control"),
         Longitudinal_StudyDesign = if_else(BioProject %in% longitudinal_projects,
                                           "Yes",
                                           "No"),
         Comorbidity_StudyDesign = if_else(BioProject %in% comorbidity_projects,
                                           "Yes",
                                           "No"),
         Treatment_StudyDesign = if_else(BioProject %in% treatment_projects,
                                         "Yes",
                                         "No"),
         Vaccination_StudyDesign = if_else(BioProject %in% vaccination_projects,
                                         "Yes",
                                         "No"))

# # Comparisons -------------------------------------------------------------
comparisons = map_dfr(.x = manually_annotated_disease_status %>%
          filter(Disease_Status_manual_annotation != other_design_annotation) %>%
          pull(BioProject) %>%
          unique(), .f = ~ (manually_annotated_disease_status %>% 
                              filter(BioProject==.x) %>%
                              pull(Disease_Status_manual_annotation) %>% 
                              unique() %>%
                              combn(2, simplify = F) %>% 
                              do.call(rbind, .) %>% 
                              as.data.frame() %>%
                              filter(V1 == control_annotation | V2 == control_annotation,
                                     V1 != other_design_annotation, V2 != other_design_annotation) %>%
                              mutate(BioProject = .x))) %>%
  mutate(Case = if_else(V1 == control_annotation, V2, V1),
         Control = control_annotation) %>%
  select(-V1, -V2)


# Write -------------------------------------------------------------------


fwrite(manually_annotated_disease_status,
       file = file.path(out_tables, "Classifications.Metadata.samples.new.v2.nonewlines.csv"),
       quote = T, row.names = F, scipen = 999)

fwrite(manually_annotated_disease_status %>%
         select(-Run, -GEO_Accession..exp.) %>%
         group_by(across(everything())) %>%
         tally %>%
         # distinct %>%
         arrange(BioProject),
       file = file.path(out_tables, "Classifications.Metadata.new.v2.nonewlines.csv"),
       quote = T, row.names = F, scipen = 999)

fwrite(comparisons %>%
         inner_join(manually_annotated_disease_status %>%
                      distinct(BioProject, PairedData, TechicalReplicates, MultipleTissuesData)),
       file = file.path(out_tables, "Classifications.Metadata.comparisons.new.v2.nonewlines.csv"),
       quote = T, row.names = F, scipen = 999)



#************************************************************************************************************