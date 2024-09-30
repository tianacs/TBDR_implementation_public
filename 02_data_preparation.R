# Data preparation and cleaning for Nanopore tNGS manuscript
# 17.07.2024

# Install pacman if needed
if (!require("pacman"))install.packages ("pacman")
# load packages
pacman::p_load(pacman, tidyverse, here, openxlsx, janitor, caret)

# Epi2Me output for each sequencing run contains a JSON file, containing all the run
# results and analysis interpretation, and a summary CSV file with some of the predicted
# resistances. I have extracted relevant information and merged the output of the different
# runs with separate short python scripts.

# Here we clean the sequencing data, by removing controls, selecting the appropriate results,
# if sequencing was repeated for failed samples, cleaning up names, specifying additional information

# Prepare reference test results (local DST) ----
# Online Trackers (updated regularly)
nicd_tracker_df <- readxl::read_excel(
  here ("data/NICD_Sample_Tracker_20240904.xlsx"), na = "N/A") %>% 
  janitor::clean_names() 

# Check for duplicates
nicd_tracker_df %>% 
  group_by(lab_id) %>% 
  mutate (count = n()) %>% 
  filter (count > 1) #%>% view()

# CLEAN and MODIFY data for consistency 
nicd_tracker_df <- nicd_tracker_df %>% 
  # Detect (and remove) duplicates
  group_by(lab_id) %>% 
  distinct(lab_id, .keep_all = TRUE)

# RedCap data (export regulary)
redcap_df <- read_csv (here("data//1734TuberculosisPati-TBtNGSsequencing_DATA_2024-09-25_0908.csv")) %>% 
  janitor::clean_names() 
# Get the labelled extract
redcap_labs_df <- read_csv (here("data//1734TuberculosisPati-TBtNGSsequencing_DATA_LABELS_2024-09-25_0909.csv")) %>% 
  janitor::clean_names()


## Xpert MTB/RIF Ultra ----
cidrz_xpert <- redcap_labs_df %>% 
  filter(event_name == "Baseline (Arm 1: RSA, Zambia)",
         data_access_group == "Global_Zambia") %>% 
  rename(sequencing_id = record_id,
         #lab_xpert_sample_type = xpert_tb_test_1_type_of_sample,
         lab_xpert_mtb_result = xpert_tb_test_1_mtb_result,
         lab_xpert_mtb_category = xpert_tb_test_1_result_category,
         lab_xpert_mtb_ct = xpert_tb_test_1_result_ct_cycle_threshold_value,
         lab_xpert_rif = xpert_tb_test_1_rif_resistance) %>% 
  select(sequencing_id, starts_with("lab"))  %>% 
  mutate(across(c(lab_xpert_mtb_result, lab_xpert_mtb_category, lab_xpert_rif), ~ tolower(.)))

# One sample not on Redcap but recorded in Tracker: 3511-225
corrected_values <- list (
  sequencing_id = "3511-225", 
  lab_xpert_mtb_result = "detected (mtb+)",
  lab_xpert_mtb_category = "low",
  lab_xpert_mtb_ct = 20,
  lab_xpert_rif = "not detected"
)
cidrz_xpert[cidrz_xpert$sequencing_id == "3511-225", ] <- corrected_values

nicd_xpert <- nicd_tracker_df %>% 
  mutate(sample_volume = case_when (tngs == FALSE ~ "insufficient")) %>% 
  select (lab_id, study_id, 
          starts_with( "lab_xpert_ultra"),
          -lab_xpert_ultra_sample_type,
          sample_volume) %>% 
  rename (sequencing_id = lab_id,
          #lab_xpert_sample_type = lab_xpert_ultra_sample_type,
          lab_xpert_mtb_result = lab_xpert_ultra_mtb_result,
          lab_xpert_mtb_category = lab_xpert_ultra_mtb_result_category,
          lab_xpert_mtb_ct = lab_xpert_ultra_mtb_ct_value,
          lab_xpert_rif = lab_xpert_ultra_rif_resistance,) %>% 
  mutate(across(c(lab_xpert_mtb_result, lab_xpert_mtb_category, lab_xpert_rif), ~ tolower(.)))

xpert_df <- bind_rows(cidrz_xpert, nicd_xpert) %>% 
  mutate(across(c(lab_xpert_mtb_result, lab_xpert_mtb_category, lab_xpert_rif), ~ as.factor(.)),
         # Add information to sequencing results on eligibility
         is_eligible = case_when (
           lab_xpert_mtb_result == "not detected (mtb-)" ~ FALSE,
           lab_xpert_mtb_category %in% c("trace", "very low") ~ FALSE,
           lab_xpert_mtb_category %in% c("low", "medium", "high") ~ TRUE
         ))

# Eligible samples have a Xpert Ultra bacterial load category of "high", "medium" or "low"
xpert_mini <- xpert_df %>% 
  select (sequencing_id,
          lab_xpert_mtb_result,
          lab_xpert_mtb_category, 
          sample_volume,
          is_eligible)
  
# CHECK THIS WITH CIDRZ - 3511-308 and 3511-367 are "MTB not detected" but have a GXPU category of "Medium" or "Low"

## Xpert MTB/XDR information ----
nicd_xdr_df <- nicd_tracker_df %>% 
  select (lab_id, study_id, red_cap_id,
          starts_with( "xpert_xdr")) %>% 
  rename (sequencing_id = lab_id)

## Phenotypic DST results ----
# 
recode_function_tbcult <- function(x) {
  case_when(
    x == 1 ~ "Positive MTB",
    x == 2 ~ "Positive NTM",
    x == 3 ~ "Contaminated",
    x == 4 ~ "Pending",
    x == 0 ~ "Negative (sterile)",
    x == 99	~ "Unknown",
    TRUE ~ NA 
  ) 
}

recode_function_dst <- function(x) {
  case_when(
    x == 1 ~ "Sensitive", 
    x == 2 ~ "Resistant",
    x == 0 ~ "Failed", 
    x == 99 ~ "Unknown",
    TRUE ~ NA
  )
}

redcap_dst_df <- redcap_df %>% 
  filter(redcap_event_name == "baseline_arm_1",
         redcap_data_access_group == "global_zambia") %>% 
  select( record_id,
          mb_tbcult_t1_result,
          mb_drug_inh,
          mb_drug_rif,
          mb_drug_pza,
          mb_drug_emb,
          mb_drug_sm) %>%
  # applying the recode_function defined above across a set of columns
  mutate(mb_tbcult_t1_result = recode_function_tbcult(mb_tbcult_t1_result),
         across (starts_with ("mb_drug"), ~ recode_function_dst(.))) %>% 
  rename (sequencing_id = record_id, 
          mgit_result = mb_tbcult_t1_result)

# Preparing the DST results from the NICD tracker
nicd_dst_df <- nicd_tracker_df %>% 
  select(lab_id,
         starts_with ("culture_")) %>% 
  rename(sequencing_id = lab_id,
         mgit_result = culture_result,
         mb_drug_inh = culture_inh_resistance,
         mb_drug_rif = culture_rif_resistance,
         mb_drug_emb = culture_emb,
         mb_drug_pza = culture_pza_resistance) %>% 
  # Correct spelling mistake
  mutate (mgit_result = case_when (mgit_result == "Positve MTB" ~ "Positive MTB",
                                   TRUE ~ mgit_result))

dst_df <- bind_rows(redcap_dst_df, nicd_dst_df)

## Save DST datasets ---- 

#write_csv(xpert_df, here("01_data/clean_xpert.csv"))
#write_csv(nicd_xdr_df, here("01_data/clean_xdr.csv"))
#write_csv(dst_df, here("01_data/clean_dst.csv"))


# Prepare index test (TBDR) resuls ----
## Load data from new output ----
df_csv  <- read_csv("data/wf_tb_amr_v2.0.0-alpha4_csv.csv")
df_json <- read_csv("data/wf_tb_amr_v2.0.0-alpha4_json.csv")
df_spol <- read_csv("data/wf_tb_amr_v2.0.0-alpha4_spoligo.csv")

df_full <-
  full_join (
    df_csv,
    df_json %>% select (-sample_pass, -sample_coverage),
    by = c("sample" = "sample_id", "experiment")
  ) %>%
  # rename drug columns
  rename_with(.cols = ends_with (".x"),
              .fn = ~ str_replace(., ".x", "_variant")) %>%
  rename_with(.cols = ends_with (".y"),
              .fn = ~ str_replace(., ".y", "_result"))

df_full <-
  full_join(df_full,
            df_spol,
             by = c("sample" = "sample_id", "experiment"))

## Add experiment information and correct sample_type info ----
seq_df <- df_full %>%
  rename(sample_id = sample) %>% 
  mutate(site = case_when(
    str_detect(sample_id, "3511") ~ "CIDRZ",
    str_detect(sample_id, "YA") ~ "NICD"
  ),) %>%
  # Remove any controls and other samples
  filter(!is.na (site),
         !experiment %in% c("IeDEA_JNB_20231101_EXP011_25GT") # run was repeated with data available
         ) %>%
  # Correct sample names
  mutate(sample_id = case_when(
    str_detect(sample_id, "-sediment",) ~ gsub("-sediment", "_x_sediment", sample_id),
    str_detect(sample_id, "-sputum",) ~ gsub("-sputum", "_x_sputum", sample_id),
    TRUE ~ sample_id
  )) %>%
  # Cleaning sample name by separating sequencing id, study id and sample type from name
  separate(
    sample_id,
    into = c("sequencing_id", "study_id", "sample_type"),
    sep = "_",
    remove = FALSE
  ) %>%
  # Correcting sample names:
  mutate(
    sequencing_id = case_when(
      sequencing_id == "3511-083" ~ "3511-83",
      sequencing_id == "3511-003" ~ "3511-3",
      TRUE ~ sequencing_id
    ),
    study_id = case_when(study_id == "x" ~ NA,
                          TRUE ~ study_id),
    sample_type = case_when(
      sample_type == "sediments" ~ "sediment",
      str_detect(sample_type, "sediment") ~ "sediment",
      str_detect(sample_type, "sputum") ~ "sputum",
      TRUE ~ sample_type
    ),
    study_id = as.numeric(study_id)
  ) %>%
  mutate(
    sample_type = case_when(
      experiment %in% c(
        "IeDEA_JNB_20230828_EXP007_25GT",
        "IeDEA_JNB_20230828_EXP008_25GT",
        "IeDEA_LSK_20230823_25GT_SEDIMENTS_EXP013",
        "IeDEA_LSK_20231025_25GT_SEDIMENTS_EXP016",
        "IeDEA_JNB_20240703_EXP018_25GT"
      ) ~ "sediment",
      # EXP015_sediments from CIDRZ, all are sediments except for 3511-310
      experiment == "IeDEA_LSK_20230907_25GT_SEDIMENTS_EXP015" &
        sequencing_id == "3511-310" & is.na(sample_type) ~ "sputum",
      experiment == "IeDEA_LSK_20230907_25GT_SEDIMENTS_EXP015" &
        is.na(sample_type) ~ "sediment",
      # MANUAL correction: Exp20 sample ID file was updated
      experiment == "IeDEA_LSK_20231220_25GT_EXP020" &
        sequencing_id %in% c("3511-254", "3511-255", "3511-258", "3511-261") ~ "sediment",
      is.na(sample_type) ~ ("sputum"),
      TRUE ~ sample_type
    )
  ) %>%
  # Remove unnecessary columns
  select(-sample_id)

## Adding additional information ----
# Count number of failed targets w/o qc sample (competitive internal control) or hsp65
targets_missed <- seq_df %>%
  select (sequencing_id,
          sample_type,
          experiment,
          qc_status,
          failed_targets) %>%
  separate_rows (failed_targets, sep = ";") %>%
  filter (
    !str_detect (failed_targets, "qc"),
    !str_detect (failed_targets, "hsp65"),
    failed_targets != "-"
  ) %>%
  group_by(sequencing_id, sample_type, experiment, qc_status) %>%
  summarise (dr_targets.missed = n())

# Add number of missed targets (targets_missed) to data
seq_df <- seq_df %>%
  left_join (targets_missed,
             by = c("sequencing_id", "sample_type", "experiment", "qc_status")) %>%
  mutate (
    # Remove QC target from failed targets
    #failed_targets2 = str_remove_all(failed_targets, ";*qc:\\d+\\.\\d+;*"),
    # Count number of failed targets
    targets_missed = case_when (
      qc_status == "fail" ~ NA,
      #if sample is "fail", targets_missed is NA
      is.na(failed_targets) ~ 0,
      failed_targets == "" ~ 0,
      TRUE ~ str_count(failed_targets, ";") +
        1
    ),
    # Specify if sample is complete or not
    status_detailed = case_when (
      qc_status == "fail" ~ "fail",
      dr_targets.missed != 0 ~ "incomplete",
      TRUE ~ "complete"
    ),
    dr_targets.missed = case_when (
      status_detailed == "complete" ~ 0,
      qc_status == "fail" ~ NA,
      TRUE ~ dr_targets.missed
    )
  )

## Group and only select most successful sample per sample ----
seq_unique <- seq_df %>%
  group_by(sequencing_id, sample_type) %>% 
  
  # Rank samples based on my criteria
        # 1) Prioritise samples with a pc_control outcome of TRUE
        # 2) Select the sample with the fewest targets missed
        # 3) If they are equal - select the "latest" sequencing run
  arrange(desc(pc_control), 
          targets_missed, 
          desc(experiment)) %>% 
  slice(1) %>% 
  ungroup () %>% 
  ## Add eligibility (GXPU) information)
  left_join (xpert_mini, by = "sequencing_id")
  

## Save cleaned dataset of sequencing data of eligible samples  ----
tngs_df <- seq_unique %>%
  filter(is_eligible | experiment == "IeDEA_JNB_20240703_EXP018_25GT",
          pc_control,
          ntc_control) %>% 
  mutate(supplementary = experiment == "IeDEA_JNB_20240703_EXP018_25GT") %>% 
  ## Add if paired (sediment and sputum available)
  group_by(sequencing_id) %>%
  mutate (is_paired = n_distinct (sample_type) > 1) %>%
  ungroup ()
  
# TO FOLLOW UP: Two experiments from South Africa passed at the sites, but the positive control
# failed with the new pipeline - why?

tngs_df %>% summarise(n_distinct(sequencing_id))

# Cleanup and save files
write_csv (tngs_df, file = "data/01_clean_tngs.csv")

# DST results in long format ----

## TBDR tNGS data ----
tngs_long <- tngs_df %>%
  select(
    sequencing_id,
    site,
    sample_type,
    qc_status,
    sample_coverage,
    ends_with (c("_result", "_variant")),
    -lab_xpert_mtb_result,
  ) %>% 
  # Create a row for each drug
  pivot_longer(
    cols = ends_with(c("_result", "_variant")),
    names_to = c("drug", "type"),
    names_sep = "_"
  ) %>%
  pivot_wider(names_from = type,
               values_from = value) %>%
  mutate(result = case_when (qc_status == "fail" ~ "fail",
                              sample_coverage == FALSE ~ "fail",
                              result == "no_resistance_detected" ~ "sensitive",
                              TRUE ~ result),
          variant = case_when(result == "fail" ~ "fail", 
                              result == "undetermined" ~ "undetermined",
                              TRUE ~ variant),
          drug = case_when (drug == "ETO" ~ "ETH",
                            TRUE ~ drug)
  )

write_csv(tngs_long, here("data/03_clean_tngs_long.csv"))

## Xpert MTB/RIF Ultra results ----
xpert_rif <- xpert_df %>%
  select(sequencing_id,
         lab_xpert_rif) %>%
  mutate(
    lab_xpert_rif = case_when (lab_xpert_rif == "detected" ~ "resistant",
                               lab_xpert_rif == "not detected" ~ "sensitive",
                               TRUE ~ lab_xpert_rif),
    drug = "RIF")



## BD BacTec MGIT phenotypic DST ----
dst_long <- dst_df %>%
  mutate_at(vars(starts_with ("mb_drug")), tolower) %>%
  # Rename columns
  rename(
    INH = mb_drug_inh,
    RIF = mb_drug_rif,
    EMB = mb_drug_emb,
    PZA = mb_drug_pza,
    STM = mb_drug_sm
  ) %>%
  # Create individual rows per drug
  pivot_longer(
    cols = c(INH, RIF, EMB, PZA, STM),
    names_to = "drug",
    values_to = "dst_result"
  )

## Add EXP018 pDST information
pDST_df.raw <- readxl::read_xlsx("data/IeDEA_JNB_20240703_EXP018_25GT_pDST_WGS.xlsx") 

exp018_pDST_long <- pDST_df.raw %>% 
  rename(sequencing_id = alias,
         mgit_result = `Culture result`) %>% 
  select(-barcode) %>% 
  pivot_longer(cols = -(c(sequencing_id, mgit_result)),names_to = "drug", values_to = "dst_result") %>%
  mutate(drug = str_sub(drug, 2),
         # Rename drugs
         drug = case_when(drug == "LEV" ~ "LFX",
                          TRUE ~ drug),
         wgs_variant = case_when(dst_result %in% c("R", "S") ~ "-",
                                 TRUE ~ dst_result),
         dst_result = case_when(dst_result == "R" ~ "resistant",
                                dst_result == "S" ~ "sensitive",
                                TRUE ~ "resistant"),
         mgit_result = case_when(mgit_result == "Pos" ~ "Positive MTB",
                                 TRUE ~ mgit_result)) 

dst_long.combined <- 
  bind_rows(dst_long, exp018_pDST_long)

rm(pDST_df.raw)

## Xpert MTB/XDR ----
xdr_long <- nicd_xdr_df %>%
  #mutate_at (vars(starts_with ("mb_drug")), tolower) %>%
  # Rename columns
  rename(
    xdr_culture = xpert_xdr_mtb_result_category,
    INH = xpert_xdr_inh_resistance,
    # For crossreferencing with tNGS data - renaming it to this group
    LFX = xpert_xdr_flq_resistance,
    AMK = xpert_xdr_amk_resistance,
    KAN = xpert_xdr_kan_resitance,
    CAP = xpert_xdr_cap_resistance,
    ETH = xpert_xdr_eth_resistance
  ) %>%
  # Added this for combined result, since XDR reports FQs 
  mutate (MFX = LFX) %>% 
  pivot_longer(
    cols = c("INH", "LFX", "MFX", "AMK", "KAN", "CAP", "ETH"),
    names_to = "drug",
    values_to = "xdr_result"
  ) %>%
  mutate(
    xdr_result = case_when(
      xdr_result == "Not detected" ~ "sensitive",
      xdr_result == "Detected" ~ "resistant",
      TRUE ~ xdr_result
    )
  ) %>% 
  ungroup()



# Combine index and reference test results for DTA analysis ----
## For index test (TBDR) need one row per sample and drug (results from sediment and sputum in wide format) ----
## Results either from sputum, from sediment
tngs_dta <- tngs_long %>%
  # Compare sediment and sputum
  pivot_wider(
    names_from = sample_type,
    values_from = c(qc_status, sample_coverage, result, variant),
    names_sep = "."
  ) 

## Combine and save DTA data ----
dta_df <- tngs_dta %>%
  left_join(xpert_rif, by = c("sequencing_id", "drug")) %>%
  left_join(dst_long.combined, by = c("sequencing_id", "drug")) %>%
  left_join(xdr_long %>% select (-red_cap_id, -study_id),
             by = c("sequencing_id", "drug"))

write_csv(dta_df, here("data/04_clean_dta_data.csv"))

# Demographic and clinical data (Table xxx) ----
## Sequencing and record ID matching ----
id.key  <- read_csv (here("data/id.key.nicd_v2.csv"))

# Add ID info to sequencing information
eligible.seq <- tngs_df %>%
  select (sequencing_id, site, sample_type, status_detailed) %>%
  # Pivot for one line per sample
  pivot_wider(names_from = sample_type,
              values_from = status_detailed)

eligible.seq <- eligible.seq %>%
  # Match RedCap IDs with Lab IDs (sequencing ID)
  left_join (id.key, by = c("sequencing_id" = "lab_id")) %>%
  # For Zambia (3511-xxx), sequencing ID is the redcap ID
  # For SA (YAxxxxxxx), some clinical data manually collected and encoded with sequencing ID
  mutate (record_id = case_when (is.na (record_id) ~ sequencing_id,
                                 TRUE ~ record_id))

## Load all available clinical information ----
redcap_demog_df <-
  redcap_labs_df %>% 
  filter (event_name == "Baseline (Arm 1: RSA, Zambia)") %>%
  select (
    record_id,
    study_identifier_if_applicable,
    patient_recruited_as_part_of,
    age_calculated,
    sex_at_birth,
    hiv_status,
    history_of_active_tb,
    data_access_group
  ) %>%
  mutate (age = trunc (age_calculated))

additional.clinical.info <-
  # Additional information collected from file review
  read_csv (here ("data/additional_clinical_data.csv")) %>%
  select (-REDCAPID,
          -StudyID) %>%
  rename (
    record_id = lab_id,
    age = Age,
    sex_at_birth = Gender,
    hiv_status = `HIV Status`,
    history_of_active_tb = `Previous TB`,
    file = File
  )

clinical.data <-
  bind_rows(redcap_demog_df, additional.clinical.info) %>%
  mutate (
    study_identifier_if_applicable =
      case_when (
        !is.na (study_id_HE2RO) ~ study_id_HE2RO,
        TRUE ~ study_identifier_if_applicable
      ),
    hiv_status =
      case_when (
        hiv_status == "Negative" ~ "Negative HIV-Test",
        hiv_status == "Positive" ~ "Positive HIV-Test",
        is.na (hiv_status) ~ "Unknown",
        TRUE ~ hiv_status
      )
  ) %>%
  select (-study_id_HE2RO)

## Combining and saving demographic and clinical data ----
demographic_data <- eligible.seq %>%
  left_join(clinical.data, by = "record_id")

## Manual additions (follow up by Denise 13.09.2024)
demographic_data <- demographic_data %>% 
  mutate(sex_at_birth = case_when(sequencing_id == "YA00494121" ~ "Male",
                                  sequencing_id == "YA00515952" ~ "Male",
                                  sequencing_id == "YA00515953" ~ "Male",
                                  sequencing_id == "YA00517146" ~ "Female",
                                  TRUE ~ sex_at_birth))


write_csv(demographic_data, here("data/02_clean_demog.csv"))


# Clean up ----
rm(df_csv, df_json, df_species_spol, 
   targets_missed,
   nicd_tracker_df,
   cidrz_xpert, nicd_xpert, nicd_dst_df,
   redcap_df, redcap_demog_df, redcap_dst_df,
   id.key,
   eligible.seq,
   clinical.data,
   additional.clinical.info,
   recode_function_dst, recode_function_tbcult,redcap_labs_df)

print ("All lines executed")



