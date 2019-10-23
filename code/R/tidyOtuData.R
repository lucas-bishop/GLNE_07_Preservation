# Loading dependencies ----------------------------------------------------

library(tidyverse)



# Sourcing if needed ------------------------------------------------------



# Functions ---------------------------------------------------------------

# Function for cleaning up raw shared files from mothur
get_tidy_shared <- function(mothur_shared_df){
  
  # Tidying mothur control shared file
  shared <- mothur_shared_df %>% 
    select(sample = Group, matches("Otu\\d+")) %>% # Renaming 'Group' col and selecting Otu cols
    gather(key = otu, value = count, contains("Otu")) # Tidying df

  return(shared)
  
}



# Function for cleaning up raw taxonomy files from mothur
get_tidy_tax <- function(mothur_tax_df){
  
  # Tidying mothur taxonomy file and separating classifications
  tax <- mothur_tax_df %>% 
    rename_all(tolower) %>% # Formatting colnames
    mutate(taxonomy = str_replace_all(taxonomy, "(\\(\\d+\\)|\\;$)", "")) %>% # Removing confidence scores and trailing ';'
    separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = ";") %>% # Breaking tax into levels
    gather(key = level, value = classification, -otu, -size) # Tidying up
  
  return(tax)
  
}



# Function for cleaning up raw alpha diversity files from mothur (sobs and shannon)
get_tidy_alpha <- function(mothur_alpha_df){
  
  alpha <- mothur_alpha_df %>% 
    select(group, method, sobs, shannon)
  
  return(alpha)
  
}



# Function for calculating relative abundances then combining tidied shared and taxonomy files into one long df
get_otu_shared_tax <- function(mothur_shared_df, mothur_tax_df){

  # Tidying shared file and calculating relative abundances
  shared <- get_tidy_shared(mothur_shared_df)

  # Tidying tax file
  tax <- get_tidy_tax(mothur_tax_df)

  # Combining the tidied dfs
  shared_tax <- inner_join(shared, tax)

  return(shared_tax)

}



# Examples ----------------------------------------------------------------

# # Getting tidied shared file
# otu_shared <- get_tidy_shared(mothur_shared)

# # Getting tidied taxonomy file
# otu_tax <- get_tidy_tax(mothur_tax)

# # Getting tidied alpha diversity file
# otu_alpha <- get_tidy_alpha(mothur_alpha)

# # Combining the shared and taxonomy files for downstream analysis
# otu_shared_tax <- get_otu_shared_tax(mothur_shared, mothur_tax)

