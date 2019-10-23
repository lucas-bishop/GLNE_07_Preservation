# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Sourcing if needed ------------------------------------------------------

source("code/R/tidyOtuData.R")



# Functions ---------------------------------------------------------------

# Function for calculating the relative abundance of Otu reads
calc_otu_rel_abund <- function(otu_shared_tax_df){
  
  # Uses tidied share file as input
  rel_abund <- otu_shared_tax_df %>% 
    group_by(sample) %>% # Groups by sample col
    mutate(total_count = sum(count)/6, # Calculating total reads per sample and appending to df (need to divide by 6 because 6 tax levels per sample)
           rel_abund = count/total_count*100) %>%  # Calculating rel abund for each Otu and appending to df
    ungroup()
  
  return(rel_abund)
  
}



# Function for aggregating relative abundances for Otus with same tax classification
aggregate_otu_rel_abund <- function(otu_rel_abund_df, tax_level_chr){
    
  # Creating df of aggregated relative abundances
  agg_otu_rel_abund <- otu_rel_abund_df %>% 
    filter(level == tax_level_chr) %>% 
    group_by(sample, classification) %>% # Aggregates classifications based on cage and day
    summarize(agg_rel_abund = sum(rel_abund)) %>% # Calculating aggregated rel abund
    ungroup()
    
  return(agg_otu_rel_abund)
  
}



# Function for finding top tax classifications for a given tax level
get_top_otu_tax <- function(agg_otu_rel_abund_df){
  
  # Finding top Otu tax classifications
  top_otu_tax <- agg_otu_rel_abund_df %>% 
    group_by(classification) %>%
    summarize(median_rel_abund = median(agg_rel_abund)) %>% # Calculating median aggregated relative abundance across all samples
    top_n(median_rel_abund, n = 9) %>% # Choosing the top Otus based on median rel abund across samples
    # filter(median_rel_abund > 0) %>% # Only choosing Otus with non-zero median relative abundances
    pull(classification) # Creating vector of Otu tax classifications
  
  return(top_otu_tax)
  
}



# Function to reassign classifications to "Other" if not within the top tax for a given tax level
get_top_otu_data <- function(otu_shared_tax_df, tax_level_chr){
  
  # Calculating relative abundances for each Otu
  otu_rel_abund <- calc_otu_rel_abund(otu_shared_tax_df)
  # Aggregating relative abundances for a given tax level and Otus with same tax classification
  agg_otu_rel_abund <- aggregate_otu_rel_abund(otu_rel_abund, tax_level_chr)
  # Pulling top tax classifications based on median aggregated relative abundances
  top_otu_tax <- get_top_otu_tax(agg_otu_rel_abund)
  
  # Reassigning any tax classification as "Other" if not in the top_otu_tax list
  # If the df consists of samples (no controls)...
  if (any(str_detect(colnames(otu_rel_abund), "^cage$"))){
    
    top_otu_data <- agg_otu_rel_abund %>% 
      mutate(classification = case_when(!(classification %in% top_otu_tax) ~ "Other",
                                        TRUE ~ classification)) %>% 
      group_by(cage, day, classification) %>% # Need to collapse all of the newly labelled "Other" tax into single entries
      summarize(agg_rel_abund = sum(agg_rel_abund)) %>% # Adds up all entries for "Other" classification
      ungroup()
  
  # If the df consists of controls (no samples)...
  } else {
    
    top_otu_data <- agg_otu_rel_abund %>% 
      mutate(classification = case_when(!(classification %in% top_otu_tax) ~ "Other",
                                        TRUE ~ classification)) %>% 
      group_by(sample, classification) %>% # Need to collapse all of the newly labelled "Other" tax into single entries
      summarize(agg_rel_abund = sum(agg_rel_abund)) %>% # Adds up all entries for "Other" classification
      ungroup()
    
  }
  
  return(top_otu_data)
  
}



# Examples ----------------------------------------------------------------

# # Calculating Otu relative abundances
# otu_rel_abund <- calc_otu_rel_abund(otu_shared_tax)

# # Aggregating relative abundances for Otus with same tax
# agg_rel_abund <- aggregate_otu_rel_abund(otu_rel_abund, "genus")

# # Creating vector of top Otu taxa for a given level
# top_otu_tax <- get_top_otu_tax(agg_rel_abund)

# # Creating top tax df with all tax not in top list labelled as "Other" for downstream analysis/plotting
# get_top_otu_data(otu_shared_tax, "genus")

