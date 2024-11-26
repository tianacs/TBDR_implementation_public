dta_edit <- dta_df %>% 
  rename (result = result.sediment) %>% 
  mutate(drug = factor (drug, levels = drug_order))

simplified_dta <- dta_edit %>% 
  select(sequencing_id, 
         site, 
         drug, 
         result, 
         lab_xpert_rif,
         dst_result, 
         xdr_result) %>% 
  mutate(result.clean = case_when(result == "resistant" ~ "resistant",
                                  result == "sensitive" ~ "sensitive"))


tab <- table(simplified_dta$result.clean, simplified_dta$lab_xpert_rif)

library(epiR)

result <- epi.tests(tab)
result
result$detail[3,2]

simp.cidrz <- simplified_dta %>% 
  filter(site == "CIDRZ")

tab.cidrz <- table(simp.cidrz$result.clean, simp.cidrz$lab_xpert_rif)
epi.tests(tab.cidrz)

simp.nicd <- simplified_dta %>% 
  filter(site == "NICD")

tab.nicd <- table(simp.nicd$result.clean, simp.nicd$lab_xpert_rif)
epi.tests(tab.nicd)
