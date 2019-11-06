library(tidyverse)
library(readxl)


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

## Group High Risk Normal with Normal samples
meta <- read.delim('data/process/clinical.txt', header=T, sep='\t', na.strings=c("","NA")) %>%
  select(sample_ID, crfDx, cdePolyps, crfClscpyRsltsPlpsAdNum, crfClscpyRsltsPlpsAdSzMin, crfClscpyRsltsPlpsAdSzMax) %>% 
  distinct() %>% 
  filter(crfDx != "Pending") 


full <- inner_join(meta, plates_total, by="sample_ID") %>% 
  mutate(diagnosis = case_when(crfClscpyRsltsPlpsAdNum > 3 & crfDx == "Adenoma" ~ "AdvAdenoma",
                        crfClscpyRsltsPlpsAdSzMin > 1 & crfDx == "Adenoma" ~ "AdvAdenoma",
                        crfClscpyRsltsPlpsAdSzMin > 1 & crfDx == "Adenoma" ~ "AdvAdenoma",
                        crfDx == "High Risk Normal" ~ "Normal",
                        crfDx == "Cancer" ~ "Carcinoma",
                        crfDx == "Normal" ~ "Normal", 
                        TRUE ~ "Adenoma")) %>% 
  select(diagnosis, sample_ID, tube_no)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

conditions <- read.table('../GLNE_07_Preservation/data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.axes',
                        header=T,
                        row.names=1) %>%
  mutate(condition = substrRight(rownames(.), 1)) %>% 
  select(condition)

distances <- read.table('../GLNE_07_Preservation/data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.axes',
                        header=T,
                        row.names=1) %>%
  mutate(tube_no = rownames(.)) %>% 
  bind_cols(conditions)


full <- inner_join(full, distances, by="tube_no") %>% 
  head(6)

  
ggplot(full) +
  aes(x=axis1, y=axis2, col=condition) +
  geom_point()
