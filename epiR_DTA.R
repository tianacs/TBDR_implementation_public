# Simplified DTA analysis


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
