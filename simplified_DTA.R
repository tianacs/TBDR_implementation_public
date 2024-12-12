# Simplified DTA analysis including stratification

# Combining the previous information from simplified_DTA.R and temp_stratified_analysis.R

# Function to calculate TP, FP, FN, FN counts
calculate_counts <- function(data, index_test, ref_test, group = NULL){
  count_data <- data %>% 
    # use "group by" conditionally
    {if(!is.null(group)) group_by(., !!sym(group)) else .} %>% 
    summarise(TP =sum(!!sym(index_test) =="resistant" & !!sym(ref_test) == "resistant", na.rm = T),
              FP =sum(!!sym(index_test) =="resistant" & !!sym(ref_test) == "sensitive", na.rm = T),
              FN =sum(!!sym(index_test) =="sensitive" & !!sym(ref_test) == "resistant", na.rm = T),
              TN =sum(!!sym(index_test) =="sensitive" & !!sym(ref_test) == "sensitive", na.rm = T),
              .groups = "drop" # Ensures the output is not grouped
    )
  # count_data <- count_data %>% 
  #   mutate(test = row_name)
  return(count_data)
}

# Define order of TB-drugs
drug_order <- c("RIF", "INH", "EMB", "PZA",
                "STM", "AMK", "CAP", "KAN", 
                "ETH", "LFX", "MXF", "BDQ", 
                "CFZ", "LZD", "DLM", "PMD")

# Define the column headers
headers <- c("Reference Test (drug)", "TP", "FP", "FN", "TN", 
             "Sensitivity (95% CI)", "Specificity (95% CI)", 
             "Sensitivity (95% CI)", "Specificity (95% CI)")


# To be adjusted in 03_analysis.Rmd or in 02_data_preparation.R
test_df <- dta_df %>% 
  mutate(result.sediment = case_when(result.sediment == "resistant" ~ "resistant",
                                     result.sediment == "sensitive" ~ "sensitive",
                                     TRUE ~ NA),
         result.sputum = case_when(result.sputum == "resistant" ~ "resistant",
                                   result.sputum == "sensitive" ~ "sensitive",
                                   TRUE ~ NA),
         result.xpert = case_when(lab_xpert_rif == "resistant" ~ "resistant",
                                  lab_xpert_rif == "sensitive" ~ "sensitive",
                                  TRUE ~ NA),
         result.dst = case_when(dst_result == "resistant" ~ "resistant",
                                dst_result == "sensitive" ~ "sensitive",
                                TRUE ~ NA),
         result.xdr = case_when(xdr_result == "resistant" ~ "resistant",
                                xdr_result == "sensitive" ~ "sensitive",
                                TRUE ~ NA),
         result.combined = case_when(
           # 1. if any test is "resistant" return "resistant"
           lab_xpert_rif == "resistant" | dst_result == "resistant" | xdr_result == "resistant" ~ "resistant",
           # 2. if any result is "sensitive" return "sensitive
           lab_xpert_rif == "sensitive" | dst_result == "sensitive" | xdr_result == "sensitive" ~ "sensitive"
         )) %>% 
  select(sequencing_id,
         drug,
         starts_with("result."))


# Create Figure 4: DTA between TBDR (sedi) compared to separate tests
combined.overall <- calculate_counts(test_df, index_test = "result.sediment", ref_test = "result.combined") %>% 
  mutate(test = "Total")
combined.drugs <- calculate_counts(test_df, index_test = "result.sediment", ref_test = "result.combined", group = "drug") %>% 
  mutate(drug = factor(drug, levels = drug_order),
         test = drug) %>% 
  arrange(drug)

forest_data <- bind_rows(
  combined.overall, 
  combined.drugs
)

Forest (forest_data, 
        study = forest_data$test,
        se.axis = c(0,1), sp.axis = c(0,1),
        col.headers = headers)

# Create a supplementary figure to stratify by GXPU
b <- test_df %>%
  left_join(xpert_df %>% select(sequencing_id,
                                lab_xpert_mtb_category),
            by = "sequencing_id") %>% 
  mutate(lab_xpert_mtb_category = case_when(lab_xpert_mtb_category == "high" ~ "high",
                                            lab_xpert_mtb_category == "medium" ~ "medium",
                                            TRUE ~ "low or lower"),
         lab_xpert_mtb_category = factor(lab_xpert_mtb_category, levels = c("high", 
                                                                            "medium", 
                                                                            "low or lower")))

b.counts <- b %>%
  filter(drug == "RIF") %>% 
  calculate_counts("result.sediment", "result.combined", "lab_xpert_mtb_category") %>% 
  mutate(test = as.character(lab_xpert_mtb_category)) 

Forest (b.counts, 
        study = b.counts$test,
        se.axis = c(0,1), sp.axis = c(0,1),
        col.headers = headers)


b.counts <- b %>%
  filter(lab_xpert_mtb_category == "low or lower") %>% 
  calculate_counts("result.sediment", "result.combined", group = "drug") %>% 
  mutate(drug = factor(drug, levels = drug_order),
         test = as.character(drug)) %>% 
  arrange(drug)

Forest (b.counts, 
        study = b.counts$test,
        se.axis = c(0,1), sp.axis = c(0,1),
        col.headers = headers)





