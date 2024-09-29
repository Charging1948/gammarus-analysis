# Required Libraries
library(tidyverse)
library(janitor) # for "clean_names" function
library(ggplot2) # for plotting
library(readxl)
library(ggthemes)
library(dplyr)
library(lubridate)
library(survival)
library(survminer)
library(flexsurv)
library(drc)

# Remove all old stuff from environment
rm(list = ls())

ld_df <- data.frame()
lc_path <- "./results/lc50"

# Set to false to iterate faster
print_plots <- FALSE
save_plots <- FALSE

plot_width <- 13
plot_height <- 8.5

# Specify the subfolder name (e.g., "plots")
subfolder <- "plots"

# Create the subfolder if it doesn't exist
if (!dir.exists(subfolder)) {
  dir.create(subfolder)
}

# Create the subfolder if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Function to read all data inside given excel-sheet out of given file
# also cleans names of headers and converts date column to correct type
process_data <- function(file_path, sheet) {
  
  # Step 1: Read the file and clean the data
  data <- read_excel(file_path, sheet = sheet)
  data <- clean_names(data)
  
  # Step 2: Convert 'date' to Date format
  data$date <- as.Date(data$date)
  
  # Step 3: Return the processed data
  return(data)
}

# Function to process the data
process_temps <- function(data) {
  
  result <- data %>%
    group_by(date, tub) %>%
    summarise(mean_tub = mean(m_temp, na.rm = TRUE)) %>%
    group_by(date) %>%
    summarise(mittelwert_all_tubs = mean(mean_tub, na.rm = TRUE))
  
  return(result)
}

generate_compParm <- function(processed_data, temps, title_file) {
  # Find first and last recorded date of dataset
  first_date <- min(processed_data$date)
  last_date <- max(processed_data$date)
  
  all_sites <- sort(unique(processed_data$site))
  
  all_dates <- seq(from = first_date, to = last_date, by = "day")
  
  # Calculate time until death
  processed_data$time_until_death <- as.numeric(difftime(processed_data$date, first_date, units = "days")) + 1
  
  # All events are deaths (1)
  processed_data$event <- 1
  
  processed_data <- processed_data %>%
    mutate(
      dose = ifelse(treatment == "mt", 1, 0)
    )
  
  # Calculate total number of individuals per site and treatment
  total_individuals <- processed_data %>%
    group_by(site, treatment) %>%
    summarise(total_ind = n(), .groups = 'drop')
  
  all_combinations <- expand.grid(site = unique(processed_data$site),
                                  treatment = unique(processed_data$treatment),
                                  date = all_dates)

  # Calculate the number of individuals per site, treatment, and date
  daily_mortality <- processed_data %>%
    group_by(site, treatment, date) %>%
    summarise(daily_ind = n(), .groups = 'drop') %>%
    left_join(total_individuals, by = c("site", "treatment")) %>%
    left_join(processed_temps, by = "date") %>%
    mutate(daily_mortality_rate = daily_ind / total_ind)
  
  full_daily_mortality <- all_combinations %>%
    left_join(daily_mortality, by = c("site", "treatment", "date")) %>%
    mutate(
      daily_mortality_rate = ifelse(is.na(daily_mortality_rate), 0, daily_mortality_rate),
      daily_ind = ifelse(is.na(daily_ind), 0, daily_ind)
    ) %>%
    group_by(date) %>%
    fill(mittelwert_all_tubs, .direction = "downup") %>%
    ungroup()
  
  full_daily_mortality <- full_daily_mortality %>%
    group_by(site, treatment) %>%
    fill(total_ind, .direction = "downup") %>%
    ungroup()
  
  # Calculate cumulative mortality over time for each site and treatment
  cumul_mortality <- full_daily_mortality %>%
    group_by(site, treatment) %>%
    arrange(date) %>%
    left_join(processed_temps, by = "date", suffix = c("", ".new")) %>%
    mutate(
      mittelwert_all_tubs = mittelwert_all_tubs.new,
      cumulative_mortality = cumsum(daily_mortality_rate),
      date_label = paste(format(date, "%d.%m.%Y"), "\n", round(mittelwert_all_tubs, 2), "°C"),
      dose = ifelse(treatment == "mt", 1, 0)
    ) %>%
    ungroup()
  
  
  splitted_data <- cumul_mortality %>%
    split(cumul_mortality$dose)
  
  for (dataset in splitted_data) {
    current_dose <- dataset$dose[1]
    # Assuming you have multiple sites and want to compare treatment effects across them
    # Fit a dose-response model with site as a factor
    dr_model_site <- drm(cumulative_mortality ~ `mittelwert_all_tubs`, site,fct = LL.2(), data = dataset, type = "binomial")
    
    lc50 <- ED(dr_model_site, 50)
    
    # Split the row names by ":" and extract the second element
    r_names <- sapply(strsplit(rownames(lc50), split = ":"), function(x) x[2])
    
    ld_sub_df <- data.frame(
      t_run = title_file,
      t_dose = current_dose,
      t_site = r_names,
      t_lc = lc50
    ) %>%
      arrange(t_site) %>%
      rename(
        Run = t_run,
        Treatment = t_dose,
        Site = t_site,
        LC50_Temperature = t_lc.Estimate,
        LC50_Std_Error = t_lc.Std..Error
      ) %>%
      mutate(
        LC50_Temperature = round(LC50_Temperature, 2),
        LC50_Std_Error = round(LC50_Std_Error, 2),
        Treatment = ifelse(Treatment == 1, "Yes", "No"),
      )
    
    ld_df <<- rbind(ld_df, ld_sub_df)
    
    compParm(dr_model_site, "e", "/")
  }
}

generate_significance <- function(processed_data, processed_temps, treatment_to_generate = NULL, title_file, site_to_generate = NULL) {
  # Find first and last recorded date of dataset
  first_date <- min(processed_data$date)
  last_date <- max(processed_data$date)
  
  all_sites <- sort(unique(processed_data$site))
  
  all_dates <- seq(from = first_date, to = last_date, by = "day")
  
  # Calculate time until death
  processed_data$time_until_death <- as.numeric(difftime(processed_data$date, first_date, units = "days"))
  
  # All events are deaths (1)
  processed_data$event <- 1
  
  processed_data <- processed_data %>%
    mutate(
      treatment_group = ifelse(treatment == "mt", "thiacloprid", "no thiacloprid")
    )
  
  if (!is.null(site_to_generate)) {
    processed_data <- processed_data %>%
      filter(site == site_to_generate)
  }
  
  if (!is.null(treatment_to_generate)) {
    processed_data <- processed_data %>%
      filter(treatment == treatment_to_generate)
  } else {
    log_rank_test <- survdiff(Surv(time_until_death, event) ~ treatment, data = processed_data)
    
    # View the result of the log-rank test
    print(log_rank_test)
  }
  
  # Visualize Kaplan-Meier survival curves
  if (is.null(treatment_to_generate)) {
    km_fit <- survfit(Surv(time_until_death, event) ~ treatment_group, data = processed_data)
  }
  else {
    km_fit <- survfit(Surv(time_until_death, event) ~ site, data = processed_data)
  }
  g <- ggsurvplot(km_fit,
                  data = processed_data,
                  pval = TRUE,
                  conf.int = TRUE,
                  risk.table = TRUE,
                  ggtheme = theme_minimal() +
                    theme(
                      panel.background = element_rect(fill = "white"),
                      plot.background = element_rect(fill = "white"),
                      panel.grid.major = element_line(color = "lightgray"),
                      panel.grid.minor = element_line(color = "lightgray")
                    )) +
    ggtitle(paste("Kaplan-Meier Survival Plot |", title_file,
                  ifelse(is.null(treatment_to_generate), "", paste("(", treatments_titles[treatment_to_generate], ")")),
                  ifelse(is.null(site_to_generate), "", paste("( Site", site_to_generate, ")"))
                  ))
  
  return(g)
}

generate_plot <- function(processed_data, processed_temps, treatment_to_generate, plot_title) {
  
  # Find first and last recorded date of dataset
  first_date <- min(processed_data$date)
  last_date <- max(processed_data$date)
  
  all_sites <- sort(unique(processed_data$site))
  
  all_dates <- seq(from = first_date, to = last_date, by = "day")
  
  # Filter the data by the selected treatment
  filtered_data <- processed_data %>%
    filter(treatment == treatment_to_generate)
  
  # Calculate total number of individuals per site and treatment
  total_individuals <- processed_data %>%
    filter(treatment == treatment_to_generate) %>%
    group_by(site, treatment) %>%
    summarise(total_ind = n(), .groups = 'drop')
  
  all_combinations <- expand.grid(site = unique(processed_data$site),
                                  treatment = c(treatment_to_generate),
                                  date = all_dates)
  
  # Calculate the number of individuals per site, treatment, and date
  daily_mortality <- filtered_data %>%
    group_by(site, treatment, date) %>%
    summarise(daily_ind = n(), .groups = 'drop') %>%
    left_join(total_individuals, by = c("site", "treatment")) %>%
    left_join(processed_temps, by = "date") %>%
    mutate(daily_mortality_rate = daily_ind / total_ind)
  
  full_daily_mortality <- all_combinations %>%
    left_join(daily_mortality, by = c("site", "treatment", "date")) %>%
    mutate(
      daily_mortality_rate = ifelse(is.na(daily_mortality_rate), 0, daily_mortality_rate),
      daily_ind = ifelse(is.na(daily_ind), 0, daily_ind)
    ) %>%
    group_by(date) %>%
    fill(mittelwert_all_tubs, .direction = "downup") %>%
    ungroup()
  
  full_daily_mortality <- full_daily_mortality %>%
    group_by(site, treatment) %>%
    fill(total_ind, .direction = "downup")
  
  # Calculate cumulative mortality over time for each site and treatment
  cumulative_mortality <- full_daily_mortality %>%
    group_by(site, treatment) %>%
    arrange(date) %>%
    left_join(processed_temps, by = "date", suffix = c("", ".new")) %>%
    mutate(
      mittelwert_all_tubs = mittelwert_all_tubs.new,
      cumulative_mortality = cumsum(daily_mortality_rate),
      date_label = paste(format(date, "%d.%m.%Y"), "\n", round(mittelwert_all_tubs, 2), "°C")
    )
  
  seq_breaks <- seq(from = 0.0, to = 1.0, by = 0.25)
  seq_labels <- paste(seq_breaks * 100, "%")
  
  
  # Generate the plot
  p <- ggplot(cumulative_mortality, aes(x = date, y = cumulative_mortality, color = as.factor(site), group = site, shape = as.factor(site))) +
    geom_line() +
    geom_point(size = 3) +
    scale_shape_manual(values = c(0, 1, 2, 3)) +  # Custom shapes for different sites
    scale_color_manual(values = c("red", "green", "blue", "purple")) +  # Custom colors
    scale_x_date(labels = cumulative_mortality$date_label, breaks = cumulative_mortality$date) +  # Set custom labels and breaks for the x-axis
    scale_y_continuous(labels = seq_labels, breaks = seq_breaks) +  # Set custom labels and breaks for the y-axis
    labs(
      title = plot_title,
      x = "Date / Mean Temperature of all Tubs",
      y = "Cumulative Mortality",
      color = "Site",  # Change legend title for colors
      shape = "Site"   # Change legend title for shapes
    ) + 
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "lightgray"),
      panel.grid.minor = element_line(color = "lightgray")
    )
  
  # Return the plot
  return(p)
}

# List of files and treatments
# INFO: Windows
# files <- c("C:\\Uni_ffm\\Bachelor\\Protokoll\\Rohdaten_1_neu.xlsx",
#            "C:\\Uni_ffm\\Bachelor\\Protokoll\\Rohdaten_2.xlsx",
#            "C:\\Uni_ffm\\Bachelor\\Protokoll\\Rohdaten_3_neu.xlsx")
files <- c("./data/Rohdaten_1_neu.xlsx",
           "./data/Rohdaten_2.xlsx",
           "./data/Rohdaten_3_neu.xlsx")
file_titles <- setNames(c("Winter", "Spring", "Summer"), files)

# All possible treatments
treatments <- c("ot", "mt")

# Create an object where one can access a title of given treatment using `treatments_titles["mt"]`
treatments_titles <- setNames(c("no thiacloprid", "thiacloprid"), treatments)



# Loop over files, process data once, and generate plots for each treatment
for (file in files) {
  # Process the data once per file
  int_data <- process_data(file, ifelse(file == files[3], "Tabelle1 (2)", "Tabelle1"))
  processed_temps <- process_temps(int_data)
  
  int_data <- int_data %>%
    filter(!is.na(individual))
  
  title_file <- file_titles[file]
  
  generate_compParm(int_data, processed_temps, title_file)
  
  for (site in unique(int_data$site)) {
    
    # Generate significance of a specific site, for comparing the treatments
    significance_site <- generate_significance(int_data, processed_temps, title_file = title_file, site_to_generate = site)
    
    if (print_plots) {
      # Display plot
      print(significance_site)
    }
    
    if (save_plots) {
      # Optionally, save plot to file
      output_filename <- paste("significance_", tolower(title_file), "_site-", site, ".png", sep = "")
      ggsave(file.path(subfolder, output_filename), plot = significance_site$plot, height = plot_height, width = plot_width)
    }
  }
  
  # Get the plot that shows the statistical significance for a run, comparing between the two different treatments
  significance <- generate_significance(int_data, processed_temps, title_file = title_file)
  
  if (print_plots) {
    # Display plots
    print(significance)
  }
  
  if (save_plots) {
    # Optionally, save plot to file
    output_filename <- paste("significance_", tolower(title_file), ".png", sep = "")
    ggsave(file.path(subfolder, output_filename), plot = significance$plot, height = plot_height, width = plot_width)
  }
  
  # Generate plots for each treatment using the processed data
  for (treatment in treatments) {
    
    # Get the plot-title for the current treatment
    title_treatment <- treatments_titles[treatment]
    # Insert treatment text into plot_title surrounded by ( and )
    plot_title <- paste(title_file, "(", title_treatment, ")")
    
    # Generate plot of given treatment
    plot <- generate_plot(int_data, processed_temps, treatment_to_generate = treatment, plot_title = plot_title)
    
    if (print_plots) {
      # Display plot
      print(plot)
    }
    
    if (save_plots) {
      # Optionally, save plot to file
      output_filename <- paste("plot_", tolower(title_file), "_", treatment, ".png", sep = "")
      ggsave(file.path(subfolder, output_filename), plot = plot, height = plot_height, width = plot_width)
    }
  }
}

lc_split <- split(ld_df, ld_df$Run)

for (run_key in names(lc_split)) {
  value <- lc_split[[run_key]]
  value <- split(value, value$Treatment)
  
  for (treatment_key in names(value)) {
    sub_value <- value[[treatment_key]]
    
    write.csv2(sub_value, file = paste(lc_path, "_run-", run_key, "_treatment-", treatment_key, ".csv", sep = ""), row.names = FALSE)
  }
}
