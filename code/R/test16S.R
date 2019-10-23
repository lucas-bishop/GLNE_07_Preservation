# Test script for verifying script functions/pipelines

library(tidyverse)



# Reading in raw data files
mothur_shared <- read_tsv("data/process/final.0.03.subsample.shared")
mothur_tax <- read_tsv("data/process/final.taxonomy")


source("code/R/tidyOtuData.R")
# code/R/tidyOtuData.R
otu_shared <- get_tidy_shared(mothur_shared)
otu_tax <- get_tidy_tax(mothur_tax)
otu_shared_tax <- get_otu_shared_tax(mothur_shared, mothur_tax)


source("code/R/calcTopOtu.R")
# code/R/calcTopOtu.R
otu_rel_abund <- calc_otu_rel_abund(otu_shared_tax)
agg_rel_abund <- aggregate_otu_rel_abund(otu_rel_abund, "genus")
# get_top_otu_tax(agg_rel_abund)
get_top_otu_data(otu_shared_tax, "genus")

summary(otu_shared_tax)

# mock_rel_abund <- calc_otu_rel_abund(mock_shared_tax)
# agg_mock_rel_abund <- aggregate_otu_rel_abund(mock_rel_abund, "class")
# get_top_otu_tax(agg_mock_rel_abund)
# get_top_otu_data(mock_shared_tax, "genus")



source("code/R/plotOtuTax.R")
# code/R/plotOtuTax.R
# plot_top_otu_bar(otu_shared_tax, "genus")
plot_top_otu_area(otu_shared_tax, "genus")

# plot_top_otu_bar(mock_shared_tax, "genus")
# plot_top_otu_area(mock_shared_tax, "genus")



# source("code/R/plotOtuDiv.R")
# # code/R/plotOtuDiv.R
# plot_otu_alpha_line(otu_alpha, "shannon")
# plot_otu_alpha_sina(otu_alpha, "sobs", "median")


# OTU plotting ------------------------------------------------------------

library(cowplot)
library(scales)

test <- otu_shared_tax %>%
  get_top_otu_data("genus")

# Function for calculating the derivative of change in rel abund by comparing to the day before
calc_otu_dxdt <- function(top_otu_data_df) {

  # Calculating dxdt across time series
  dxdt_df <- top_otu_data_df %>% 
    group_by(cage, classification) %>% # Groups for calculating derivatives
    arrange(day) %>% # Arranging the rows by time just in case
    mutate(dxdt = case_when(day == min(day) ~ 0, # Sets dxdt = 0 for the first time point as baseline
                            TRUE ~ (agg_rel_abund - lag(agg_rel_abund, 1))/(day - lag(day, 1)))) %>% # Calculates dxdt ((change rel abund)/(change time)) from day to day
    ungroup()

  # Returning the calculated values attached to df
  return(dxdt_df)

}

# Testing function
calc_otu_dxdt(test) %>% 
  arrange(cage, classification)



# Function for creating line plots of 16S relative abundances from top OTUs
plot_otu_line <- function(otu_shared_tax_df, tax_level_chr) {
  
  # Generating the top Otu df; anything not in top list for given tax level is labelled as "Other"
  top_plot_data <- get_top_otu_data(otu_shared_tax_df, tax_level_chr) %>% 
    mutate(classification = str_replace_all(classification, "_", " "), # Formatting legend classification labels immediately before plotting
           classification = str_replace(classification, "(unclassified)", "(\\1)")) # More label formatting
  
  # Plotting a stacked bar chart of the Otu taxa data
  line_plot <- top_plot_data %>% 
    ggplot(aes(x = day, y = agg_rel_abund, color = cage)) + # Setting up the plotting conditions
    annotate("segment", x = c(3,8), xend = c(3, 8),
             y = c(0,0), yend = c(100,100),
             lwd = 0.1, alpha = 0.75) +
    annotate("text", x = 3, y = 106, label = "Start", size = 3) +
    annotate("text", x = 8, y = 106, label = "Stop", size = 3) +
    geom_line(lwd=0.65) + # Making line plot
    scale_x_continuous(breaks = seq(1, max(top_plot_data$day), by = 2), expand = c(0,0)) + # Setting breaks to every other day
    scale_y_continuous(expand = c(0,0)) +
    coord_cartesian(xlim = c(1, max(top_plot_data$day)), ylim = c(0,100), clip = "off") + # Makes sure y axis goes from 0 to 100 (some of the samples may equal 100 with floating decimal so would be excluded if limits set)
    facet_wrap(~ classification, nrow = 2, scales = "free") + # Plotting each classification individually
    labs(#title = "16S Abundance",
         x = "Day",
         y = "Relative abundance (%)",
         color = "Group") +
    theme_classic() +
    theme(strip.text.x = element_text(face = "bold", margin = margin(b = 20)), # Increasing spacing between strip labels and chart
          strip.background.x = element_blank(), # Removing strip label box
          panel.spacing = unit(2, "lines"))
  
  return(line_plot)
  
}



# Function for creating line plots of log10 16S relative abundances from top OTUs
plot_otu_log10_line <- function(otu_shared_tax_df, tax_level_chr) {
  
  # Generating the top Otu df; anything not in top list for given tax level is labelled as "Other"
  top_plot_data <- get_top_otu_data(otu_shared_tax_df, tax_level_chr) %>% 
    mutate(log10_agg_rel_abund = log10(agg_rel_abund), # log10 transforming the relative abundances
           log10_agg_rel_abund = case_when(!is.finite(log10_agg_rel_abund) ~ NA_real_, # Replacing any non-finite values with NA to create gaps in geom_line
                                           TRUE ~ log10_agg_rel_abund), # Keep any finite values as is
           classification = str_replace_all(classification, "_", " "), # Formatting legend classification labels immediately before plotting
           classification = str_replace(classification, "(unclassified)", "(\\1)")) # More label formatting

  # Generating nice looking breaks
  change_breaks <- pretty(n = 4, top_plot_data$log10_agg_rel_abund)
  
  # Plotting a stacked bar chart of the Otu taxa data
  line_plot <- top_plot_data %>% 
    ggplot(aes(x = day, y = log10_agg_rel_abund, color = cage, shape = cage)) + # Setting up the plotting conditions
    annotate("segment", x = c(3,8), xend = c(3, 8),
             y = rep(min(change_breaks), 2), yend = rep(max(change_breaks), 2),
             lwd = 0.1, alpha = 0.75) +
    annotate("text", x = 3, y = max(change_breaks)*1.06, label = "Start", size = 3) +
    annotate("text", x = 8, y = max(change_breaks)*1.06, label = "Stop", size = 3) +
    geom_hline(yintercept = 0, linetype = "42") + # Drawing line indicating no change
    geom_point(size = 1, alpha = 0.75) +
    geom_line(lwd = 0.65, alpha = 0.75) + # Making line plot
    scale_x_continuous(breaks = seq(1, max(top_plot_data$day), by = 2), expand = c(0,0)) + # Setting breaks to every other day
    scale_y_continuous(breaks = change_breaks, labels = change_breaks, expand = c(0,0)) +
    coord_cartesian(xlim = c(1, max(top_plot_data$day)), ylim = range(change_breaks), clip = "off") + # Makes sure y axis goes from 0 to 100 (some of the samples may equal 100 with floating decimal so would be excluded if limits set)
    facet_wrap(~ classification, nrow = 2, scales = "free") + # Plotting each classification individually
    labs(#title = "16S Abundance",
      x = "Day",
      y = "Log10 Relative abundance (%)",
      color = "Group") +
    theme_classic() +
    theme(strip.text.x = element_text(face = "bold", margin = margin(b = 20)), # Increasing spacing between strip labels and chart
          strip.background.x = element_blank(), # Removing strip label box
          panel.spacing = unit(2, "lines"))
  
  return(line_plot)
  
}



# Function for plotting line plots of top otu dxdt values across time series
plot_otu_dxdt_line <- function(otu_shared_tax_df, tax_level_chr) {
  
  # Generating the top Otu df; anything not in top list for given tax level is labelled as "Other"
  top_plot_data <- get_top_otu_data(otu_shared_tax_df, tax_level_chr) %>% 
    calc_otu_dxdt() %>% # Calculating derivative of top otus
    mutate(classification = str_replace_all(classification, "_", " "), # Formatting legend classification labels immediately before plotting
           classification = str_replace(classification, "(unclassified)", "(\\1)")) # More label formatting
  
  # Generating nice looking breaks
  change_breaks <- pretty(top_plot_data$dxdt)
  
  # Plotting a stacked bar chart of the Otu taxa data
  dxdt_line_plot <- top_plot_data %>% 
    ggplot(aes(x = day, y = dxdt, color = cage)) + # Setting up the plotting conditions
    annotate("segment", x = c(3,8), xend = c(3, 8),
             y = rep(min(change_breaks), 2), yend = rep(max(change_breaks), 2),
             lwd = 0.1, alpha = 0.75) +
    annotate("text", x = 3, y = max(change_breaks)*1.06, label = "Start", size = 3) +
    annotate("text", x = 8, y = max(change_breaks)*1.06, label = "Stop", size = 3) +
    geom_hline(yintercept = 0, linetype = "42") + # Drawing line indicating no change
    geom_line(lwd=0.65) + # Making line plot
    scale_x_continuous(breaks = seq(1, max(top_plot_data$day), by = 2), expand = c(0,0)) + # Setting x axis breaks to every other day
    scale_y_continuous(breaks = change_breaks, labels = change_breaks, expand = c(0,0)) + # Setting y axis breaks based on 
    coord_cartesian(xlim = c(1, max(top_plot_data$day)), ylim = range(change_breaks), clip = "off") + # Makes sure y axis goes from 0 to 100 (some of the samples may equal 100 with floating decimal so would be excluded if limits set)
    facet_wrap(~ classification, nrow = 2, scales = "free") + # Plotting each classification individually
    labs(#title = "Relative Change in 16S Abundance",
      x = "Day",
      y = "dx/dt",
      color = "Group") +
    theme_classic() +
    theme(strip.text.x = element_text(face = "bold", margin = margin(b = 20)), # Increasing spacing between strip labels and chart
          strip.background.x = element_blank(), # Removing strip label box
          panel.spacing = unit(2, "lines"))
  
  return(dxdt_line_plot)
  
}



plot_otu_line(otu_shared_tax, "genus")
plot_otu_log10_line(otu_shared_tax, "genus")
plot_otu_dxdt_line(otu_shared_tax, "genus")





# otu_line_plot <- plot_otu_line(otu_shared_tax, "genus")
# otu_change_line_plot <- plot_top_otu_change_line(otu_shared_tax, "genus")
# 
# 
# otu_plot2 <- plot_grid(plot_grid(otu_line_plot + theme(legend.position = "none"),
#                     plot.new(),
#                     otu_change_line_plot + theme(legend.position = "none"),
#                     labels = c("A", "", "B"), nrow = 3, rel_heights = c(0.47, 0.06, 0.47)),
#           get_legend(otu_line_plot),
#           nrow = 1,
#           rel_widths = c(0.9, 0.1))
# 
# otu_plot1 <- plot_top_otu_area(otu_shared_tax, "genus")
# 
# ggsave(plot = otu_plot1, filename = "results/figures/otuPlot1.png", width = 16, height = 6, units = "in", dpi = 600)
# 
# ggsave(plot = otu_plot2, filename = "results/figures/otuPlot2.png", width = 16, height = 10, units = "in", dpi = 600)


# test3 %>% 
#   ggplot(aes(x=day, y = agg_rel_abund, color = cage)) +
#   geom_line(lwd=0.65) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   facet_wrap(~ classification, nrow = 2) +
#   labs(title = "Relative 16S abundance",
#        x = "Day",
#        y = "Relative 16S abundance",
#        color = "Group") +
#   theme_classic()


# change_breaks <- pretty(test3$change_rel_abund)
# 
# test3 %>% 
#   ggplot(aes(x=day, y = change_rel_abund, color = cage)) +
#   geom_hline(yintercept = 0, linetype = "42") +
#   geom_line(lwd=0.65) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0), limits = range(change_breaks), breaks = change_breaks, labels = c("ND", change_breaks[-1])) +
#   facet_wrap(~ classification, nrow = 2) +
#   labs(title = "Change in 16S Abundance Relative to Baseline",
#        x = "Day",
#        y = "log2 transformed change in 16S abundance\nrelative to baseline samples",
#        color = "Genus") +
#   theme_classic()

# Alpha diversity ---------------------------------------------------------


mothur_alpha <- read_tsv("data/mothur/process/sample.final.groups.ave-std.summary")

tidy_alpha <- mothur_alpha %>% 
  separate(group, into = c("cage", "day"), sep = "D") %>% 
  mutate(cage = case_when(cage == "82" ~ "Cef1",
                   cage == "83" ~ "Cef2",
                   cage == "84" ~ "Mock"),
         day = as.numeric(day)) %>% 
  filter(method == "ave") %>% 
  select(cage, day, nseqs, coverage, invsimpson, shannon, sobs) %>% 
  gather(key=metric, value=value, nseqs, coverage, invsimpson, shannon, sobs)

# write_tsv(tidy_alpha, "data/mothur/process/sampleAlphaDiversity.tsv")
# 
# alpha_names <- c(
#   'coverage' = "Good's Coverage",
#   'sobs' = "OTU Richness",
#   'shannon' = "Shannon Diversity",
#   'invsimpson' = "Inverse Simpson Diversity"
# )
# 
# alpha_plot1 <- read_tsv("data/mothur/process/sampleAlphaDiversity.tsv") %>%
#   group_by(cage, day, metric) %>% 
#   summarize(mean_value = mean(value),
#             std_dev = sd(value)) %>% 
#   ungroup() %>% 
#   filter(metric != "nseqs") %>%
#   mutate(cage = factor(cage, levels = c("Mock", "Cef1", "Cef2"))) %>% 
#   ggplot(aes(x = day, y = mean_value, group = cage, color = cage)) +
#   geom_line(size=0.5, alpha = 0.75) +
#   # geom_errorbar(aes(ymin = mean_value - std_dev, ymax = mean_value + std_dev), size = 0.25, width = 0.1) +
#   scale_x_continuous(breaks = seq(1, 23, by = 2), expand = c(0,0)) + # Setting breaks to every other day
#   facet_wrap(~ metric, scales="free_y", strip.position = "left", nrow = 1, labeller = labeller(.rows = alpha_names)) +
#   labs(x = "Day",
#        y = NULL,
#        color = "Group") +
#   theme_classic() +
#   theme(panel.spacing = unit(2, "lines"),
#     strip.background = element_blank(),
#     strip.placement="outside")


metrics <- c("sobs", "coverage", "shannon", "invsimpson")

# Function for plotting alpha diversity
plot_alpha_line <- function(tidy_alpha_df, metric_chr_vec) {
  
  # Vector for labelling axes
  alpha_names <- c(
    'coverage' = "Good's Coverage",
    'sobs' = "OTU Richness",
    'shannon' = "Shannon Diversity",
    'invsimpson' = "Inverse Simpson Diversity"
  )
  
  # Iterating through supplied vector of alpha diversity metrics to generate plots and append into single list
  alpha_plots <- list()
  for (i in metric_chr_vec) {
    
    # Filter the data for each metric and calculate average
    alpha_data <- tidy_alpha_df %>%
      filter(metric == i) %>% 
      group_by(cage, day) %>% 
      summarize(mean_value = mean(value),
                std_dev = sd(value)) %>% 
      ungroup() %>% 
      mutate(cage = factor(cage, levels = c("Mock", "Cef1", "Cef2")))
    
    alpha_breaks <- pretty(alpha_data$mean_value)
    
    # Iteratively creating plots based on filtered data
    alpha_plots[[i]] <- alpha_data %>%
      ggplot(aes(x = day, y = mean_value, group = cage, color = cage)) +
      annotate("segment", x = c(3,8), xend = c(3, 8),
               y = rep(min(alpha_breaks), 2), yend = rep(max(alpha_breaks), 2),
               lwd = 0.1, alpha = 0.75) +
      annotate("text", x = 3, y = max(alpha_breaks), label = "Start", size = 3) + # NEED BE CAREFUL WITH RANGE/GOODS COVERAGE VALUES
      annotate("text", x = 8, y = max(alpha_breaks), label = "Stop", size = 3) +
      geom_line(size=0.5, alpha = 0.75) +
      # geom_errorbar(aes(ymin = mean_value - std_dev, ymax = mean_value + std_dev), size = 0.25, width = 0.1) +
      scale_x_continuous(breaks = seq(1, 23, by = 2), expand = c(0,0)) + # Setting breaks to every other day
      scale_y_continuous(breaks = alpha_breaks,
                         expand = c(0,0)) +
      coord_cartesian(xlim = c(1, max(alpha_data$day)), ylim = range(alpha_breaks), clip = "off") + 
      labs(x = NULL,
           y = alpha_names[i],
           color = "Group") +
      theme_classic() +
      theme(plot.margin = unit(c(2,1,1,1), "lines"),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none")
    
  }
  
  # alpha_plots_margins <- list()
  # for (i in metric_chr_vec) {
  #   
  #   if (i == metric_chr_vec[length(metric_chr_vec)]) {
  #     
  #     alpha_plots_margins[[i]] <- alpha_plots[[i]] + theme(plot.margin = unit(c(1,0.5,1,1.825), "lines"))
  #       
  #   } else  if (i == metric_chr_vec[1]) {
  #     
  #     alpha_plots_margins[[i]] <- alpha_plots[[i]] + theme(plot.margin = unit(c(1,0.9,1,0.75), "lines")) 
  #       
  #   } else {
  #     
  #     alpha_plots_margins[[i]] <- alpha_plots[[i]] + theme(plot.margin = unit(c(1,1,1,1.325), "lines"))
  #     
  #   }
  #   
  # }
  
  # Arranging all the plots and attaching shared x axis label and legend
  joined_alpha_plots <- ggdraw(add_sub(plot_grid(plotlist = alpha_plots, nrow = 1),
                                       "Day", vpadding=grid::unit(0,"lines"),
                                       y = 6, x = 0.512, vjust = 4.5, size = 11.5))
  
  # Returning the final plot arrangements
  return(joined_alpha_plots)
  
}


alpha_plot1 <- plot_alpha_line(tidy_alpha, metrics)




# EXTRACT MAIN PLOTTING PORTION AS SEPARATE FUNCTION
# CALCULATE LOG2 BEFORE CALCULATING MEAN/SD


# Function for plotting the relative changes in alpha diversity compared to baseline
plot_alpha_change_line <- function(tidy_alpha_df, metric_chr_vec) {
  
  # Vector for labelling axes
  alpha_change_names <- c(
    'coverage' = "Change in Good's Coverage relative to baseline (log2)",
    'sobs' = "Change in OTU Richness relative to baseline (log2)",
    'shannon' = "Change in Shannon Diversity relative to baseline (log2)",
    'invsimpson' = "Change in Inverse Simpson Diversity relative to baseline (log2)"
  )
  
  # Iterating through supplied vector of alpha diversity metrics to generate plots and append into single list
  alpha_change_plots <- list()
  for (i in metric_chr_vec) {
    
    # Filter the data for each metric and calculate relative change
    alpha_change_data <- tidy_alpha_df %>%
      filter(metric == i) %>% 
      group_by(cage, day) %>% 
      summarize(mean_value = mean(value),
                std_dev = sd(value)) %>% 
      ungroup() %>% 
      group_by(cage) %>% 
      # mutate(change_value = log10(mean_value/mean(mean_value[day %in% c("1", "2", "3")]))) %>% 
      mutate(change_value = log2(mean_value/mean(mean_value[day %in% c("1", "2", "3")]))) %>%
      ungroup() %>% 
      mutate(cage = factor(cage, levels = c("Mock", "Cef1", "Cef2")))
    
    alpha_change_breaks <- pretty(alpha_change_data$change_value)
    
    # Iteratively creating plots based on filtered data
    alpha_change_plots[[i]] <- alpha_change_data %>%
      ggplot(aes(x = day, y = change_value, group = cage, color = cage)) +
      annotate("segment", x = c(3,8), xend = c(3, 8),
               y = rep(min(alpha_change_breaks), 2), yend = rep(max(alpha_change_breaks), 2),
               lwd = 0.1, alpha = 0.75) +
      annotate("text", x = 3, y = max(alpha_change_breaks), label = "Start", size = 3) +
      annotate("text", x = 8, y = max(alpha_change_breaks), label = "Stop", size = 3) +
      geom_hline(yintercept = 0, linetype = "42") +
      geom_line(size=0.5, alpha = 0.75) +
      # geom_errorbar(aes(ymin = mean_value - std_dev, ymax = mean_value + std_dev), size = 0.25, width = 0.1) +
      scale_x_continuous(breaks = seq(1, 23, by = 2), expand = c(0,0)) + # Setting breaks to every other day
      scale_y_continuous(breaks = alpha_change_breaks,
                         expand = c(0,0)) +
      coord_cartesian(xlim = c(1, max(alpha_change_data$day)), ylim = range(alpha_change_breaks), clip = "off") + 
      labs(x = NULL,
           y = alpha_change_names[i],
           color = "Group") +
      theme_classic() +
      theme(plot.margin = unit(c(2,1,1,1), "lines"),
            axis.title.y = element_text(size = 8, margin = margin(r = 8)),
            legend.position = "none")
    
  }
  
  # alpha_change_plots_margins <- list()
  # for (i in metric_chr_vec) {
  #   
  #   if (i == metric_chr_vec[length(metric_chr_vec)]) {
  #     
  #     alpha_change_plots_margins[[i]] <- alpha_change_plots[[i]] + theme(plot.margin = unit(c(1,0.5,1,1.825), "lines"))
  #       
  #   } else  if (i == metric_chr_vec[1]) {
  #     
  #     alpha_change_plots_margins[[i]] <- alpha_change_plots[[i]] + theme(plot.margin = unit(c(1,0.9,1,0.75), "lines")) 
  #       
  #   } else {
  #     
  #     alpha_change_plots_margins[[i]] <- alpha_change_plots[[i]] + theme(plot.margin = unit(c(1,1,1,1.325), "lines"))
  #     
  #   }
  #   
  # }
  
  # Arranging all the plots and attaching shared x axis label and legend
  joined_alpha_plots <- ggdraw(add_sub(plot_grid(plotlist = alpha_change_plots, nrow = 1),
                                                 "Day", vpadding=grid::unit(0,"lines"),
                                                 y = 6, x = 0.512, vjust = 4.5, size = 11.5))
  
  # Returning the final plot arrangements
  return(joined_alpha_plots)
  
}

# for (i in metrics) {
#   if (i == metrics[length(metrics)]) {
#     print("yes")
#   } else {
#     print("no")
#   }
# }


# alpha_plot2 <- plot_alpha_change_line(tidy_alpha, metrics)
# 
# 
# 
# alpha_plots <- plot_grid(plot_grid(alpha_plot1 + theme(legend.position = "none"),
#                     plot.new(),
#                     alpha_plot2 + theme(legend.position = "none"),
#                     labels = c("A", "", "B"), nrow = 3, rel_heights = c(0.47, 0.06, 0.47)),
#           get_legend(otu_line_plot),
#           nrow = 1,
#           rel_widths = c(0.9, 0.1))
# 
# 
# ggsave(plot = alpha_plots, filename = "results/figures/alphaPlots.png", width = 16, height = 8, units = "in", dpi = 600)



# Beta diversity ----------------------------------------------------------

mothur_beta_dist <- "data/mothur/process/sample.final.braycurtis.0.03.lt.ave.dist"

# # Creating function for reading in lower trianglular distance data and tidying it
# get_tidy_dist <- function(dist_filename_chr) {
#   
#   # Finding the row with the max number of entries (columns) and pulling for use setting matrix dimensions
#   n <- max(count.fields(dist_filename_chr))
#   
#   # Reading in the distance matrix
#   # Fills missing data with "NA", skips the first row (only has number of rows), sets 1st col as row names,
#   # and sets col names as 1 through value of n
#   dist_matrix <- as.matrix(read.table(dist_filename_chr, fill = TRUE, skip = 1, row.names = 1, col.names = 1:n))
#   
#   # Setting col names to be same as row names (drops the last rowname because there isn't a column for it since it would be all "NA")
#   colnames(dist_matrix) <- row.names(dist_matrix)[-n]
#   
#   # Tidying up matrix for downstream use
#   dist_df <- as_tibble(dist_matrix, rownames = "row") %>% # Converting matrix to tibble and keeping rownames as 1st col (default drops them)
#     gather(key = column, value = distance, -row) %>% # Gathers all of the columns other than row 
#     select(row, everything()) %>% 
#     filter(!is.na(distance)) # Filters out any rows where the distance has value of "NA" (missing data from original triangle matrix)
#   
#   return(dist_df)
#   
# }
# 
# test<-get_tidy_dist(mothur_beta_dist)

get_tidy_dist <- function(dist_filename_chr) {
  
  # Finding the row with the max number of entries (columns) and pulling for use setting matrix dimensions
  n <- max(count.fields(dist_filename_chr))
  
  # Reading in the distance matrix
  # Fills missing data with "NA", skips the first row (only has number of rows), sets 1st col as row names,
  # and sets col names as 1 through value of n-1. Finally, appends col of NA
  dist_matrix <- cbind(as.matrix(read.table(dist_filename_chr, fill = TRUE, skip = 1, row.names = 1, col.names = 1:n-1)), NA)
  
  # Replacing all NAs with 
  dist_matrix[is.na(dist_matrix)] <- 0
  
  # Setting col names to be same as row names (drops the last rowname because there isn't a column for it since it would be all "NA")
  colnames(dist_matrix) <- row.names(dist_matrix)
  
  # Filling in values in upper triangle using matrix addition (requires m1 rows = m2 cols and m1 cols = m2 rows)
  complete_matrix <- dist_matrix + t(dist_matrix)
  
  # Tidying up matrix for downstream use
  complete_dist_df <- as_tibble(complete_matrix, rownames = "row") %>% # Converting matrix to tibble and keeping rownames as 1st col (default drops them)
    gather(key = column, value = distance, -row) %>% # Gathers all of the columns other than row
    select(row, everything()) %>%
    filter(!is.na(distance)) # Filters out any rows where the distance has value of "NA" (missing data from original triangle matrix)
  
  return(complete_dist_df)
  
}
  
  

# get_tidy_dist(mothur_beta_dist) %>% 
#   write_tsv("data/mothur/process/sampleBetaDiversity.tsv")



mothur_beta_data <- get_tidy_dist(mothur_beta_dist) %>%
  filter(str_detect(row, c("D1$|D2$|D3$"))) %>% # Selecting only rows from time point 1 because want to compare time 2 and 3 to 1
  # filter(!str_detect(column, c("D1$|D2$|D3$"))) %>% # Selecting only rows from time points 2 and 3 because want to compare to time point 1
  separate(row, into = c("cage", "day"), sep = "D") %>% 
  separate(column, into = c("col_cage", "col_day"), sep = "D") %>% 
  filter(cage == col_cage) %>% # Limiting comparisons to within groups
  mutate(cage = case_when(cage == "82" ~ "Cef1",
                          cage == "83" ~ "Cef2",
                          cage == "84" ~ "Mock"),
         cage = factor(cage, levels = c("Mock", "Cef1", "Cef2")),
         col_cage = case_when(col_cage == "82" ~ "Cef1",
                              col_cage == "83" ~ "Cef2",
                              col_cage == "84" ~ "Mock"),
         col_cage = factor(col_cage, levels = c("Mock", "Cef1", "Cef2")),
         day = as.numeric(day),
         col_day = as.numeric(col_day)) %>% 
  group_by(cage, col_day) %>% 
  summarize(mean_dist = mean(distance)) %>% 
  ungroup() 

beta_breaks <- pretty(mothur_beta_data$mean_dist)

# Plotting by day
beta_time_plot <- mothur_beta_data %>% 
  ggplot(aes(x=col_day, y = mean_dist, color=cage)) +
  annotate("segment", x = c(3,8), xend = c(3, 8),
           y = rep(min(beta_breaks), 2), yend = rep(max(beta_breaks), 2),
           lwd = 0.1, alpha = 0.75) +
  annotate("text", x = 3, y = max(beta_breaks)*1.02, label = "Start", size = 3) +
  annotate("text", x = 8, y = max(beta_breaks)*1.02, label = "Stop", size = 3) +  
  geom_line() +
  scale_x_continuous(breaks = seq(1, 23, by = 2), expand = c(0,0)) + # Setting breaks to every other day
  scale_y_continuous(breaks = beta_breaks,
                     expand = c(0,0)) +
  coord_cartesian(xlim = c(1, max(mothur_beta_data$col_day)), ylim = range(beta_breaks), clip = "off") +
  labs(color = "Group",
       x = "Day",
       y = "Bray-Curtis distance relative to baseline") +
  theme_classic() +
  theme(plot.margin = unit(c(2,1,1,1), "lines"))


ggsave(plot = beta_time_plot, filename = "results/figures/betaPlot.png", width = 8, height = 5, units = "in", dpi = 600)


# Ordinations -------------------------------------------------------------

# Reading in axes made using filtered dist file
mothur_pcoa <- read_tsv("data/mothur/process/sample.final.braycurtis.0.03.lt.ave.pcoa.axes")
mothur_pcoa_loadings <- read_tsv("data/mothur/process/sample.final.braycurtis.0.03.lt.ave.pcoa.loadings")


pcoa_data <- mothur_pcoa %>% 
  select(group, axis1, axis2) %>% 
  separate(group, into = c("cage", "day"), sep = "D") %>%  
  mutate(cage = case_when(cage == "82" ~ "Cef1",
                          cage == "83" ~ "Cef2",
                          cage == "84" ~ "Mock"),
         cage = factor(cage, levels = c("Mock", "Cef1", "Cef2")),
         day = as.numeric(day)) %>% 
  arrange(day)


# pcoa_label_positions <- tibble(x=c(0.43, -0.7, -0.62),
#                                y=c(0.05, 0.1, 0.5),
#                                cage = c('Mock', 'Cef1', 'Cef2'),
#                                label= c('Mock', 'Cef1', 'Cef2'))

pcoa_label_positions <- tibble(x=c(0.43, -0.56, -0.37),
                               y=c(0.03, 0.1, 0.5),
                               cage = c('Mock', 'Cef1', 'Cef2'),
                               label= c('Mock', 'Cef1', 'Cef2'))


pcoa_plot <- pcoa_data %>%
  ggplot(aes(x=axis1, y=axis2, color = cage)) +
  geom_path(alpha = 0.5, show.legend = FALSE) +
  geom_point(data = subset(pcoa_data, day == "1" | day == "23"), aes(fill = as.character(day)),
             shape = 21, stroke = 1.5, size = 3, alpha = 0.5) +
  # scale_size_continuous(range = c(0.1,2)) +
  # scale_color_discrete(guide = FALSE) +
  # scale_fill_gradient(low = "#FFFFFF", high = "#000000", limits = c(1,23), breaks = c(1,23)) +
  scale_fill_manual(values = c("1" = "#FFFFFF", "23" = "#000000")) +
  # stat_ellipse(aes(group=cage), show.legend=FALSE) +
  # geom_text(data = pcoa_label_positions, aes(x, y, label = label, color = cage), inherit.aes=FALSE, show.legend=FALSE, size=3) +
  labs(x=paste0("PCo Axis 1 (", format(round(mothur_pcoa_loadings$loading[1], 2), nsmall = 2), "%)"),
       y=paste0("PCo Axis 2 (", format(round(mothur_pcoa_loadings$loading[2], 2), nsmall = 2), "%)"),
       fill = "Day",
       color = "Group") +
  # guides(fill = guide_colorbar(ticks = FALSE)) +
  theme_classic() +
  theme(
    legend.key.height = unit(0.8, "line"),
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(0,0,0,0)
  )

ggsave(plot = pcoa_plot, filename = "results/figures/pcoaPlot.png", width = 8, height = 7, units = "in", dpi = 600)





# Reading in axes made using filtered dist file
mothur_nmds <- read_tsv("data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.axes")
mothur_nmds_stress <- read_tsv("data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.stress")

nmds_data <- mothur_nmds %>% 
  select(group, axis1, axis2) %>% 
  separate(group, into = c("cage", "day"), sep = "D") %>%  
  mutate(cage = case_when(cage == "82" ~ "Cef1",
                          cage == "83" ~ "Cef2",
                          cage == "84" ~ "Mock"),
         cage = factor(cage, levels = c("Mock", "Cef1", "Cef2")),
         day = as.numeric(day)) %>% 
  arrange(day)


nmds_label_positions <- tibble(x=c(-0.37, 0.46, -0.56),
                               y=c(0.46, -0.14, 0.00),
                               cage = c('Mock', 'Cef1', 'Cef2'),
                               label= c('Mock', 'Cef1', 'Cef2'))

# nmds_label_positions <- tibble(x=c(-0.37, 0.94, -1.02),
#                                y=c(0.5, 0.03, 0.1),
#                                cage = c('Mock', 'Cef1', 'Cef2'),
#                                label= c('Mock', 'Cef1', 'Cef2'))


nmds_plot <- nmds_data %>% 
  ggplot(aes(x=axis1, y=axis2, color = cage, fill = day)) +
  geom_path(alpha = 0.5, show.legend = FALSE) +
  geom_point(shape = 21, stroke = 1, size = 3, alpha = 0.75) +
  scale_size_continuous(range = c(0.1,2)) +
  scale_color_discrete(guide = FALSE) +
  scale_fill_gradient(low = "#FFFFFF", high = "#000000", limits = c(1,23), breaks = c(1,23)) +
  annotate("text", x = -0.53, y = -0.7, label = paste0("Stress = ", min(mothur_nmds_stress$Stress))) +
  # stat_ellipse(aes(group=cage), show.legend=FALSE) +
  geom_text(data = nmds_label_positions, aes(x, y, label = label, color = cage), inherit.aes=FALSE, show.legend=FALSE, size=3) +
  labs(x="NMDS Axis 1",
       y="NMDS Axis 2",
       fill = "Day") +
  guides(fill = guide_colorbar(ticks = FALSE)) +
  theme_classic() +
  theme(
    legend.key.height = unit(0.8, "line"),
    legend.key.size = unit(0.25, "cm"),
    legend.margin = margin(0,0,0,0)
  )


plot_grid(plotlist = list(beta_time_plot,
                          plot_grid(plotlist = list(pcoa_plot + theme(legend.position = "none"),
                                                    plot.new(),
                                                    nmds_plot + theme(legend.position = "none"),
                                                    get_legend(pcoa_plot)),
                                    nrow = 1, rel_widths = c(0.44, 0.06, 0.44, 0.06),
                                    labels = c("B", "", "C", ""))),
          nrow = 2, labels = c("A", ""))

# plot_grid(plotlist = list(pcoa_plot + theme(legend.position = "none"),
#                           plot.new(),
#                           nmds_plot + theme(legend.position = "none"),
#                           get_legend(pcoa_plot)),
#           nrow = 1, rel_widths = c(0.44, 0.06, 0.44, 0.06),
#           labels = c("B", "", "C", ""))



ggsave(plot = nmds_plot, filename = "results/figures/nmdsPlot.png", width = 8, height = 7, units = "in", dpi = 600)


all_beta_plots <- plot_grid(plotlist = list(beta_time_plot,
                                            plot_grid(plotlist = list(pcoa_plot + theme(legend.position = "none"),
                                                                      ggplot() + theme_void(),
                                                                      nmds_plot + theme(legend.position = "none"),
                                                                      get_legend(pcoa_plot)),
                                                      nrow = 1, rel_widths = c(0.44, 0.06, 0.44, 0.06),
                                                      labels = c("B", "", "C", ""))),
                            nrow = 2, labels = c("A", ""))

ggsave(plot = all_beta_plots, filename = "results/figures/betaPlots.png", width = 16, height = 12, units = "in", dpi = 600)



# ggplot(aes(x=axis1, y=axis2, color = treatment, shape = time)) +
  # geom_point(alpha = 0.75) +
  # stat_ellipse(aes(group=treatment), show.legend=FALSE) +
  # geom_text(data = nmds_balbc_label_positions, aes(x, y, label=label, color=treatment), inherit.aes=FALSE, show.legend=FALSE, size=3) +
  # annotate("text", x = -0.72, y = -0.8, label = paste0("Stress = ", min(mothur_nmds_stress_balbc$Stress))) +
  # scale_color_manual(guide = 'none', values=brewer.pal(4, "Set1")) +
  # labs(title = "BALB/C",
  #      shape = "Time",
  #      x="NMDS Axis 1",
  #      y="NMDS Axis 2") +
  # theme_classic() +
  # theme(
  #   strip.background.x = element_blank(),
  #   legend.key.height = unit(0.8, "line"),
  #   legend.key.size = unit(0.06, "cm"),
  #   legend.margin = margin(0,0,0,0)
  # )
