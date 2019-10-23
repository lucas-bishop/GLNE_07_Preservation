# Loading dependencies ----------------------------------------------------

library(tidyverse)
library(RColorBrewer)



# Sourcing if needed ------------------------------------------------------

source("code/R/calcTopOtu.R")



# Functions ---------------------------------------------------------------

# Function for creating color palette (expands default color palette if necessary)
get_color_palette <- colorRampPalette(brewer.pal(12, "Paired"))


  
# Function for generating a stacked bar chart of top tax classifications for a given level
plot_top_otu_bar <- function(otu_shared_tax_df, tax_level_chr){
  
  # Generating the top Otu df; anything not in top list for given tax level is labelled as "Other"
  top_plot_data <- get_top_otu_data(otu_shared_tax_df, tax_level_chr) %>% 
    mutate(classification = str_replace_all(classification, "_", " "), # Formatting legend classification labels immediately before plotting
           classification = str_replace(classification, "(unclassified)", "(\\1)")) # More label formatting
  
  # Calculating the number of unique Otu classifications to set number of colors needed for plot
  n_classifications <- top_plot_data %>% 
    distinct(classification) %>% 
    count() %>% 
    pull()
  
  # If the supplied df is from sample data (no controls)...
  if (any(str_detect(colnames(otu_shared_tax_df), "^cage$"))){
    
    # Plotting a stacked bar chart of the Otu taxa data
    stacked_plot <- top_plot_data %>% 
      ggplot(aes(x = day, y = agg_rel_abund, fill = classification)) + # Setting up the plotting conditions
      geom_col() + # Making stacked bar chart
      annotate("segment", x = c(3.5,8.5), xend = c(3.5,8.5),
               y = c(0,0), yend = c(100,100),
               linetype = "42", alpha = 0.75) +
      # annotate("rect", xmin = 3.5, xmax = 8.5, ymin = 0, ymax = 100, fill = "white", alpha = 0.5) +
      scale_fill_manual(values = get_color_palette(n_classifications)) + # Creating color palette
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(breaks = seq(1, max(top_plot_data$day), by = 2)) + # Setting breaks to every other day
      coord_cartesian(ylim = c(0,100)) + # Makes sure y axis goes from 0 to 100 (some of the samples may equal 100 with floating decimal so would be excluded if limits set)
      facet_wrap(~ cage) + # Plotting each cage grouped individually
      labs(x = "Day", y = "Relative abundance (%)",
           title = "16S Abundance", fill = str_to_title(tax_level_chr)) +
      theme_classic() +
      theme(strip.text.x = element_text(margin = margin(b = 10)), # Increasing spacing between strip labels and chart
            strip.background.x = element_blank()) # Removing strip label box
    
  # Or if the supplied df is from controls (no samples)...
  } else {
    
    # Plotting a stacked bar chart of the Otu taxa data
    stacked_plot <- top_plot_data %>% 
      ggplot(aes(x = sample, y = agg_rel_abund, fill = classification)) + # Setting up the plotting conditions
      geom_col() + # Making stacked bar chart
      scale_fill_manual(values = get_color_palette(n_classifications)) + # Creating color palette
      scale_y_continuous(expand = c(0,0)) +
      coord_cartesian(ylim = c(0,100)) + # Makes sure y axis goes from 0 to 100 (some of the samples may equal 100 with floating decimal so would be excluded if limits set)
      labs(x = "Standard", y = "Relative abundance (%)",
           title = "16S Abundance", fill = str_to_title(tax_level_chr)) +
      theme_classic() +
      theme(strip.text.x = element_text(margin = margin(b = 10)), # Increasing spacing between strip labels and chart
            strip.background.x = element_blank()) # Removing strip label box
    
  }
  
  return(stacked_plot)
  
}



# Function for generating an area chart of top tax classifications for a given level (easier to see trends over long time courses)
# NOTE: Only works if `x` is a continuous variable (in this case 'day')
plot_top_otu_area <- function(otu_shared_tax_df, tax_level_chr){
  
  # If the supplied df doesn't have sampling over time...
  if (!any(str_detect(colnames(otu_shared_tax_df), "^day$"))){
    
    stop("This isn't the plot you are looking for. These area plots require multiple data points of the same groups over time so try a different style.")
  
  # But if it does...
  } else {
    
    # Generating the top Otu df; anything not in top list for given tax level is labelled as "Other"
    top_plot_data <- get_top_otu_data(otu_shared_tax_df, tax_level_chr) %>% 
      mutate(classification = str_replace_all(classification, "_", " "), # Formatting legend classification labels immediately before plotting
             classification = str_replace(classification, "(unclassified)", "(\\1)")) # More label formatting
    
    # Calculating the number of unique Otu classifications to set number of colors needed for plot
    n_classifications <- length(unique(top_plot_data$classification))
    
    # Plotting an area plot of the Otu taxa data
    area_plot <- top_plot_data %>% 
      ggplot(aes(x = day, y = agg_rel_abund, fill = classification)) + # Setting up the plotting conditions
      geom_area() + # Making area plot
      annotate("segment", x = c(3,8), xend = c(3,8),
               y = c(0,0), yend = c(100,100),
               lwd = 0.25, alpha = 0.75) +
      annotate("text", x = 3, y = 102, label = "Start", size = 3) +
      annotate("text", x = 8, y = 102, label = "Stop", size = 3) +
      scale_x_continuous(breaks = seq(1, max(top_plot_data$day), by = 2), expand = c(0,0)) + # Setting breaks to every other day
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values = get_color_palette(n_classifications)) + # Creating color palette
      coord_cartesian(xlim = c(1, max(top_plot_data$day)), ylim = c(0,100), clip = "off") + # Makes sure y axis goes from 0 to 100 (some of the samples may equal 100 with floating decimal so would be excluded if limits set)
      facet_wrap(~ cage, nrow = 1, scales = "free") + # Plotting each cage grouped individually
      labs(x = "Day", y = "Relative abundance (%)",
           # title = "16S Abundance",
           fill = str_to_title(tax_level_chr)) +
      theme_classic() +
      theme(strip.text.x = element_text(face = "bold", margin = margin(b = 20)), # Increasing spacing between strip labels and chart
            strip.background.x = element_blank(), # Removing strip label box
            panel.spacing = unit(2, "lines"))
    
    return(area_plot)
    
  }
  
}



# Examples ----------------------------------------------------------------

# # Creating a stacked bar chart of the tax data across samples
# plot_top_otu_bar(otu_shared_tax, "genus")

# # Creating an area chart of the tax data across samples
# plot_top_otu_area(otu_shared_tax, "genus")

