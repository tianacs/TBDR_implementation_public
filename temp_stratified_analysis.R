combined_ref.cidrz <- combined_ref %>% 
  filter(str_detect(sequencing_id, "3511"))

combined_ref.nicd <- combined_ref %>% 
  filter(str_detect(sequencing_id, "YA"))

test <- combined_ref


contingency_ref <- test %>% 
  select (drug, result, combined_result) %>% 
  filter (!(result %in% c("undetermined", "fail")),
          !is.na (combined_result), 
          !is.na (result)) %>% 
  mutate (result = as.factor (result), 
          dst_result = as.factor (combined_result))

drug_data_ref <- contingency_ref %>% 
  group_by(drug) %>% 
  summarise(TP=sum(result =="resistant" & dst_result =="resistant"),
            FP=sum(result =="resistant" & dst_result =="sensitive"),
            FN=sum(result =="sensitive" & dst_result =="resistant"),
            TN=sum(result =="sensitive" & dst_result =="sensitive")) %>% 
  #filter(!(TP == 0 & FN == 0),
  # Remove STM since not pDST tested for consistently)
  #drug != "STM") %>% 
  arrange(factor(drug, levels = c("RIF", "INH", "EMB", "PZA"))) 

overall_dst_data_ref <- contingency_ref %>% 
  summarise(TP= sum(result == "resistant" & dst_result == "resistant"),
            FP= sum(result == "resistant" & dst_result == "sensitive"),
            FN= sum(result == "sensitive" & dst_result == "resistant"),
            TN= sum(result == "sensitive" & dst_result == "sensitive")) %>% 
  mutate(drug = "Overall") 


combined_dta_table_ref <- bind_rows(overall_dst_data_ref,
                                    drug_data_ref) %>% 
  mutate(n = TP+TN+FP+FN) %>% 
  mutate(drug = recode(drug, 
                       RIF = "rifampicin (RIF)",
                       INH = "isoniazid (INH)",
                       EMB = "ethambutol (EMB)",
                       PZA = "pyrazinamide (PZA)",
                       STM = "streptomycin (STM)",
                       AMK = "amikacin (AMK)",
                       CAP = "capreomycin (CAP)",
                       KAN = "kanamycin (KAN)",
                       ETH = "ethionamide (ETH)",
                       LFX = "levofloxacin (LFX)",
                       MXF = "moxifoxacin (MFX)",
                       BDQ = "bedaquiline (BDQ)",
                       CFZ = "clofazimine (CFZ)",
                       LZD = "linezolid (LZD)",
                       DLM = "delamanid (DLM)"
  ))

# Draw forest plot
Forest(combined_dta_table_ref, 
       study = combined_dta_table_ref$drug,
       se.axis = c(0,1), sp.axis = c(0,1),
       col.headers = headers)

