# Simplified DTA analysis

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


# Create counts
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
  return(count_data)
}

# Simplified code for DTA tables ----
drug_order <- c("RIF", "INH", "EMB", "PZA",
                "STM", "AMK", "CAP", "KAN", 
                "ETH", "LFX", "MXF", "BDQ", 
                "CFZ", "LZD", "DLM", "PMD")


# Create forest plot
forest_data <- bind_rows(
  calculate_counts(test_df, "result.sediment", "result.combined") %>% 
    mutate(test = "Total"),
  calculate_counts(test_df, "result.sediment", "result.combined", "drug") %>% 
    mutate(drug = factor(drug, levels = drug_order),
           test = drug) %>% 
    arrange(drug)) 

# Define the column headers:
headers <- c("Reference Test (drug)", "TP", "FP", "FN", "TN", 
             "Sensitivity (95% CI)", "Specificity (95% CI)", 
             "Sensitivity (95% CI)", "Specificity (95% CI)")

Forest (forest_data, 
        study = forest_data$test,
        se.axis = c(0,1), sp.axis = c(0,1),
        col.headers = headers)


# For GXPU 
b <- test_df %>%
  left_join(xpert_df %>% select(sequencing_id,
                                lab_xpert_mtb_category),
            by = "sequencing_id") %>%
  calculate_counts("result.sediment", "result.combined", "lab_xpert_mtb_category") %>% 
  mutate(test = as.character(lab_xpert_mtb_category))


# Define the column headers:
headers <- c("Reference Test (drug)", "TP", "FP", "FN", "TN", 
             "Sensitivity (95% CI)", "Specificity (95% CI)", 
             "Sensitivity (95% CI)", "Specificity (95% CI)")

Forest (b, 
        study = b$test,
        se.axis = c(0,1), sp.axis = c(0,1),
        col.headers = headers)

library(epiR)

# Function to calculate sensitivity, specificity, and 95% CI
calculate_metrics <- function(data, index_col, reference_col) {
  # Create contingency table
  tbl <- table(data[[index_col]], data[[reference_col]])
  
  print(tbl)
  
  # Check if the table has at least one value in all required cells
  if (all(dim(tbl) == c(2, 2)) && all(rowSums(tbl) > 0) && all(colSums(tbl) > 0)) {
    # Calculate sensitivity, specificity, and CI
    result <- epiR::epi.tests(tbl)
    
    # Extract the relevant statistics from the detail data frame
    detail <- result$detail
    
    sensitivity <- detail %>% filter(statistic == "se")
    specificity <- detail %>% filter(statistic == "sp")
    
    # Create a data frame with the extracted values
    metrics <- data.frame(
      Sensitivity = sensitivity$est,
      Sensitivity_Lower_CI = sensitivity$lower,
      Sensitivity_Upper_CI = sensitivity$upper,
      Specificity = specificity$est,
      Specificity_Lower_CI = specificity$lower,
      Specificity_Upper_CI = specificity$upper
    )
    
    # Return with all numeric values rounded to 2 decimal points
    return(metrics %>% mutate(across(everything(), ~ round(.x, 2))))
  }
  else {
    # Return NA for invalid tables
    return(data.frame(
      Sensitivity = NA,
      Sensitivity_Lower_CI = NA,
      Sensitivity_Upper_CI = NA,
      Specificity = NA,
      Specificity_Lower_CI = NA,
      Specificity_Upper_CI = NA
    ))
  }}


test_df %>% 
  calculate_metrics(index_col = "result.sediment", reference_col = "result.xpert")

calculate_metrics(test_df, index_col = "result.sediment", reference_col = "result.dst")

calculate_metrics(test_df, index_col = "result.sediment", reference_col = "result.combined")




# Grouping doesn't work yet ----
# Group by category and loop through groups
results <- test_df %>%
  left_join(xpert_df %>% select(sequencing_id,
                                lab_xpert_mtb_category),
            by = "sequencing_id") %>%
  group_by(lab_xpert_mtb_category) %>% view()
  group_split() %>% 
  lapply(function(group_data){
    calculate_metrics(group_data, index_col = "results.sediment", reference_col = "results.dst")
  }) %>%
  bind_rows(.id = "Category")

# Print results
print(results)

# Group by category and loop through groups
results <- test_df %>%
  group_by(drug) %>%
  group_split() %>%
  lapply(function(group_data){
    calculate_metrics(group_data, index_col = "results.sediment", reference_col = "results.combined")
  }) %>%
  bind_rows(.id = "Category")


test_df %>% 
  filter(str_detect(sequencing_id, "3511")) %>% 
  calculate_metrics(index_col = "results.sediment", reference_col = "results.xpert")
