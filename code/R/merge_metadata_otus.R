######################################################################
# Author: Begum Topcuoglu
# Date: 2019-11-06
# Title: Merge meta-data and OTU information
######################################################################

######################################################################
# Description: 

# This script will read in data from excel files generated for 
# incubation experiments. Each excel sheet is for 30 samples 
# incubated in 3 conditions (A:freeze upon arrival,
# B:incubate at RT 24h, C:incubate at RT on preservation 
# buffer 24h). There are 90 rows in each excel sheet to 
# represent 30 samples in triplicates.

#     - Sample ID
#     - Tube number
#     - Diagnosis information used for randomization


# The script does:
#     1: Take each excel sheet and change the tube number to include
#       the plate number information as well.

#     2. Put all the plates in one dataframe.

#     3. Read in the meta-data file that has more clinical information
#     for each patient. (The sample_IDs are the same as the excel sheets

#     4. Merge detailed info about SRN status for each patient. 



######################################################################
# Dependencies and Outputs: 

# Be in the project directory.

# The outputs are:
#   (1) 
#   (2) 
######################################################################


################### IMPORT LIBRARIES  ################################
library(tidyverse)
library(readxl)
######################################################################



######################################################################
#     1: Take each excel sheet and change the tube number to include
#       the plate number information as well.
#     2. Put all the plates in one dataframe.
######################################################################
files <- list.files(path= 'data/process', pattern='.*_setup.xlsx', full.names = TRUE)

i <- 0
plates_total = data.frame()
for(file_name in files){
  i <- i + 1
  plate <- read_excel(file_name) %>% 
    select(sample_ID, dx, tube_no)
  plate$tube_no <- paste0("P", i, "S", plate$tube_no)
  plates_total <- rbind(plates_total , plate)
}

files <- list.files(path= 'data/process', pattern='.*_setup_final.xlsx', full.names = TRUE)

i <- 6
rest_plates_total = data.frame()
plate_total = rbind(plates_total, rest_plates_total)
for(file_name in files){
  i <- i + 1
  plate <- read_excel(file_name) %>% 
    select(sample_ID, dx, tube_no)
  plate$tube_no <- paste0("P", i, "S", plate$tube_no)
  plates_total <- rbind(plates_total , plate)
}


######################################################################
# 3. Read in the meta-data file that has more clinical information
#   for each patient. (The sample_IDs are the same as the excel sheets)
######################################################################
## Group High Risk Normal with Normal samples
meta <- read.delim('data/process/clinical.txt', 
                   header=T, 
                   sep='\t', 
                   na.strings=c("","NA")) %>%
  select(sample_ID, 
         crfDx, 
         cdePolyps, 
         crfClscpyRsltsPlpsAdNum, 
         crfClscpyRsltsPlpsAdSzMin, 
         crfClscpyRsltsPlpsAdSzMax) %>% 
  distinct() %>% 
  filter(crfDx != "Pending") 

######################################################################
# 4. Merge detailed info about SRN status for each patient. 
######################################################################
# If patients have more than 3 adenomas, adenomas larger than 1cm then
# call those as advanced adenomas. 

# High-risk normal patients are defined as normal.

full <- inner_join(meta, plates_total, by="sample_ID") %>% 
  mutate(diagnosis = case_when(crfClscpyRsltsPlpsAdNum > 3 & crfDx == "Adenoma" ~ "AdvAdenoma",
                        crfClscpyRsltsPlpsAdSzMin > 1 & crfDx == "Adenoma" ~ "AdvAdenoma",
                        crfClscpyRsltsPlpsAdSzMin > 1 & crfDx == "Adenoma" ~ "AdvAdenoma",
                        crfDx == "High Risk Normal" ~ "Normal",
                        crfDx == "Cancer" ~ "Carcinoma",
                        crfDx == "Normal" ~ "Normal", 
                        TRUE ~ "Adenoma")) %>% 
  select(diagnosis, sample_ID, tube_no) %>% 
  mutate(tube_no = str_replace_all(tube_no, "P1S", ""))

######################################################################
# 5. Bring in data from mothur analyses to merge with metadata
######################################################################

# First define a function to extract the A, B, C conditions from the 
# end of the tube number. Make that a seperate column to define which 
# condition that tube was incubated as. 

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

axes <- read.table('../GLNE_07_Preservation/data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.axes',
                        header=T,
                        row.names=1)

conditions <-  axes %>% 
  mutate(condition = substrRight(rownames(.), 1)) %>% 
  select(condition)

design <- read.table('../GLNE_07_Preservation/data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.axes', header=T) %>% 
  bind_cols(conditions) %>% 
  select(group, condition) %>% 
  mutate(condition = case_when(condition=="A" ~ "Frozen",
                               condition=="B" ~ "RoomTemp24h",
                               condition=="C" ~ "RoomTempBuffer24h")) %>% 
  write_tsv("../GLNE_07_Preservation/data/mothur/process/condition.design")

# Now we have a column for conditions.

distances <- axes %>% 
  mutate(tube_no = rownames(.)) %>% 
  bind_cols(conditions)

# Merge the mothur file with the metadata
merged <- inner_join(distances, full, by="tube_no")



