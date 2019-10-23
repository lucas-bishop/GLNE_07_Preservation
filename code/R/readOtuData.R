# Loading dependencies ----------------------------------------------------

# Loading dependencies
library(tidyverse)



# Command line interface --------------------------------------------------

# Allowing script to run from the command line
args <- commandArgs(trailingOnly = TRUE)

# Assigning command line inputs
mothur_sub_shared_filename <- args[1] # Subsampled shared file from mothur
mothur_tax_filename <- args[2] # Sample taxonomy file from mothur


read_tsv("data/mothur/process/mock.final.shared")
# Functions ---------------------------------------------------------------

if (!is.na(mothur_sub_shared_filename)){
  mothur_sub_shared <- read_tsv(mothur_sub_shared_filename)
}

mothur_shared_filename <- "data/mothur/process/final.0.03.subsample.shared"


# Function for reading in subsampled shared file from mothur
read_mothur_shared <- function(mothur_sub_shared_filename = mothur_shared_filename){
  
  mothur_sub_shared <- read_tsv(mothur_sub_shared_filename)
  
  return(mothur_sub_shared)
  
}

read_mothur_shared()

# Function for reading in taxonomy file from mothur
read_mothur_tax <- function(mothur_tax_filename){
  
  if(!is.na(mothur_tax_filename)){
    
    mothur_tax <- read_tsv(mothur_tax_filename)
    
  } else {
    
    mothur_tax <- read_tsv("data/mothur/process/final.taxonomy")
    
  }
  
  return(mothur_tax)
  
}






read_mothur_shared(mothur_sub_shared_filename)
read_mothur_shared("data/mothur/process/final.0.03.subsample.shared")

read_mothur_tax(mothur_tax_filename)



# Reading in subsampled shared file
mothur_sub_shared <- read_tsv("data/mothur/process/final.0.03.subsample.shared")

# Reading in the taxonomy file
mothur_tax <- read_tsv("data/mothur/process/final.taxonomy")



