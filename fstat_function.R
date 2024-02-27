#Thank you to Marina Vabistsevits for producing the make_mvmr_input function and Emma Hazelwood for testing and refining the function.

library(dplyr)
library(vroom)
library(MVMR)
library(TwoSampleMR)
library(readr)
library(tibble)
library(tidyr)


# Read in make MVMR function ----------------------------------------------

make_mvmr_input <- function(exposure_dat, outcome.id.mrbase="ebi-a-GCST006464", outcome.data=NULL){
  # provide exposure_dat created in the same way as for TwoSampleMR 
  # also specify the outcome argument [only ONE!] (MR-base ID or full gwas data in .outcome format)
  
  # extract SNPs for both exposures from outcome dataset
  # (for the selected option mr.base or local outcome data)
  if (!is.null(outcome.id.mrbase)) {
    # if mrbase.id is provided
    outcome_dat <- extract_outcome_data(snps = unique(exposure_dat$SNP),
                                        outcomes = outcome.id.mrbase)
  } else if (!is.null(outcome.data)){
    # if outcome df is provided
    outcome_dat <- outcome.data %>% filter(SNP %in% exposure_dat$SNP)
  }
  
  # harmonise datasets
  exposure_dat <- exposure_dat %>% mutate(id.exposure = exposure)
  outcome_harmonised <- mv_harmonise_data(exposure_dat, outcome_dat)
  
  exposures_order <- colnames(outcome_harmonised$exposure_beta)
  
  # Create variables for the analysis 
  
  ### works for many exposures
  no_exp = dim(outcome_harmonised$exposure_beta)[2] # count exposures
  
  # add beta/se names
  colnames(outcome_harmonised$exposure_beta) <- paste0("betaX", 1:no_exp)
  colnames(outcome_harmonised$exposure_se) <- paste0("seX", 1:no_exp)
  
  XGs <-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'), 
                  as.data.frame(outcome_harmonised$exposure_se)   %>%rownames_to_column('SNP'), 
                  by = "SNP")
  
  YG <- data.frame(beta.outcome = outcome_harmonised$outcome_beta,
                   se.outcome = outcome_harmonised$outcome_se) %>% 
    mutate(SNP = XGs$SNP)
  
  
  return(list(YG = YG,
              XGs = XGs,
              exposures = exposures_order))
}

# Load example datasets ---------------------------------------------------
#Give OpenGWAS IDs
exposure_ids<-c("ebi-a-GCST90092808", "ebi-a-GCST90092809", "ebi-a-GCST90092822", "ebi-a-GCST90092838", "ebi-a-GCST90092883", "ebi-a-GCST90092992", "ebi-a-GCST90093000")
outcome_id<-"ebi-a-GCST006464"

#Extract exposures for MVMR
exposure_dat <- mv_extract_exposures(exposure_ids)

#Use above function to format
mvmr_input<-make_mvmr_input(exposure_dat=exposure_dat,outcome.id.mrbase=outcome_id)


# Build a function to drop exposures until all F stats > 10 ---------------
remove_F_stat<-function(exposure_ids,outcome_id,all_pairwise_correlations,sres){
  
  vector <- seq_along(exposure_ids)
  df <- data.frame(exposure_ids, vector)
  
  df1<-df
  
  while (min(sres)<10){
    print("Minimum F stat less than 10. Removing an exposure and trying again")
    df1 <- df1[-c(which.min(sres)),]
    exposure_dat1 <- mv_extract_exposures(df1$exposure_ids)
    mvmr_input1 <- make_mvmr_input(exposure_dat1, outcome.id.mrbase=outcome_id)
    mvmr_out1 <- format_mvmr(BXGs = mvmr_input1$XGs %>% select(contains("beta")),  
                             BYG = mvmr_input1$YG$beta.outcome,                        
                             seBXGs = mvmr_input1$XGs %>% select(contains("se")),      
                             seBYG = mvmr_input1$YG$se.outcome,                        
                             RSID = mvmr_input1$XGs$SNP)
    cormat_1 <- all_pairwise_correlations[which(colnames(all_pairwise_correlations) %in% df1$exposure_ids), which(colnames(all_pairwise_correlations) %in% df1$exposure_ids)]
    se_matrix1 <- mvmr_out1 %>% as_tibble() %>% select(contains("sebetaX")) %>% as.data.frame()
    phenocov_mvmr1 <- phenocov_mvmr(Pcov = cormat_1, seBXGs = se_matrix1)
    sres <- strength_mvmr(r_input=mvmr_out1, gencov=phenocov_mvmr1)
  }
  
  return(list(mvmr_input = mvmr_input1,
              mvmr_out = mvmr_out1,
              cormat=cormat_1,
              se_matrix=se_matrix1,
              phenocov_mvmr=phenocov_mvmr1,
              sres=sres))
         
}


# Try out function --------------------------------------------------------

#Prep MVMR data
mvmr_out <- format_mvmr(BXGs = mvmr_input$XGs %>% select(contains("beta")),  # exposure betas
                        BYG = mvmr_input$YG$beta.outcome,                     # outcome beta
                        seBXGs = mvmr_input$XGs %>% select(contains("se")),  # exposure SEs
                        seBYG = mvmr_input$YG$se.outcome,                     # outcome SEs
                        RSID = mvmr_input$XGs$SNP)                            # SNPs


#Read in genetic correlations
pairwise <- read_csv("initial_correlations.csv")
all_pairwise_correlations <- as.matrix(pairwise)
se_matrix <- mvmr_out %>% as_tibble() %>% select(contains("sebetaX")) %>% as.data.frame()
colnames(all_pairwise_correlations)<-exposure_ids

#Calculate initial conditional F stats
phenocov_mvmr <- phenocov_mvmr(Pcov = all_pairwise_correlations, seBXGs = se_matrix)
sres <- strength_mvmr(r_input=mvmr_out, gencov=phenocov_mvmr)

mvmr_dat<-remove_F_stat(exposure_ids=exposure_ids,outcome_id=outcome_id,all_pairwise_correlations=all_pairwise_correlations,sres=sres)
