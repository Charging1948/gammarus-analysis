# Required Libraries
library(tidyverse)    # Includes dplyr, ggplot2
library(janitor)      # For clean_names function
library(readxl)       # For reading excel files
library(ggthemes)     # For theme_minimal in plots
library(survival)     # For survival analysis
library(survminer)    # For Kaplan-Meier survival plots
library(flexsurv)     # For flexible survival models
library(drc)          # For dose-response models

# Remove all old stuff from the environment
rm(list = ls())

# Constants
SUBFOLDER <- "plots"
RESULTS_DIR <- "results"
PLOT_WIDTH <- 13
PLOT_HEIGHT <- 8.5
LD_DF <- data.frame()
LC_PATH <- "./results/lc50"
PRINT_PLOTS <- TRUE
SAVE_PLOTS <- TRUE
GGTHEME <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_line(color = "lightgray")
  )

# List of files and treatments
# INFO: Windows
# FILES <- c("C:\\Uni_ffm\\Bachelor\\Protokoll\\Rohdaten_1_neu.xlsx",
#            "C:\\Uni_ffm\\Bachelor\\Protokoll\\Rohdaten_2.xlsx",
#            "C:\\Uni_ffm\\Bachelor\\Protokoll\\Rohdaten_3_neu.xlsx")
FILES <- c("./data/Rohdaten_1_neu.xlsx",
           "./data/Rohdaten_2.xlsx",
           "./data/Rohdaten_3_neu.xlsx")
FILE_TITLES <- setNames(c("Winter", "Spring", "Summer"), FILES)

# All possible treatments
TREATMENTS <- c("ot", "mt")

# Create an object where one can access a title of given treatment using `treatment_titles["mt"]`
TREATMENT_TITLES <- setNames(c("no thiacloprid", "thiacloprid"), TREATMENTS)

# Create necessary directories
create_dir <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
}
create_dir(SUBFOLDER)
create_dir(RESULTS_DIR)

# Function to read all data inside a given Excel sheet
# Cleans header names and converts date column to the correct type
process_data <- function(file_path, sheet) {
  data <- read_excel(file_path, sheet = sheet) %>%
    clean_names() %>%
    mutate(
      date = as.Date(date),
      time_until_death = as.numeric(difftime(.$date, min(.$date), units = "days")) + 1,
      event = 1,
      dose = ifelse(treatment == "mt", 1, 0)
    )
  
  return(data)
}

# Function to process temperatures
process_temps <- function(data) {
  data %>%
    group_by(date, tub) %>%
    summarise(mean_tub = mean(m_temp, na.rm = TRUE)) %>%
    group_by(date) %>%
    summarise(mittelwert_all_tubs = mean(mean_tub, na.rm = TRUE)) %>%
    return()
}

# Function to handle plot saving and printing
handle_plot <- function(plot, filename) {
  if (PRINT_PLOTS) {
    print(plot)
  }
  if (SAVE_PLOTS) {
    ggsave(file.path(SUBFOLDER, filename), plot = plot, height = PLOT_HEIGHT, width = PLOT_WIDTH)
  }
}

# LC50 Calculations and Saving Results for each Run/Treatment
save_lc50 <- function(ld_df, lc_path) {
  lc_split <- split(ld_df, ld_df$Run)
  
  for (run_key in names(lc_split)) {
    value <- lc_split[[run_key]]
    value <- split(value, value$Treatment)
    
    for (treatment_key in names(value)) {
      sub_value <- value[[treatment_key]]
      write.csv2(sub_value, file = paste(lc_path, "_run-", run_key, "_treatment-", treatment_key, ".csv", sep = ""), row.names = FALSE)
    }
  }
}

calculate_lc50 <- function(processed_data, processed_temps, title_file) {
  # Find first and last recorded date of dataset
  first_date <- min(processed_data$date)
  last_date <- max(processed_data$date)
  all_dates <- seq(from = first_date, to = last_date, by = "day")
  all_sites <- sort(unique(processed_data$site))
  
  cumulative_mortality <- calculate_cumulative_mortality(processed_data, processed_temps) %>%
    split(.$dose)
  
  for (dataset in cumulative_mortality) {
    current_dose <- dataset$dose[1]
    
    dr_model <- drm(cumulative_mortality ~ `mittelwert_all_tubs`, site,fct = LL.2(), data = dataset, type = "binomial")
    
    lc50 <- ED(dr_model, 50)
    
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
    
    # Bind row to global dataset
    LD_DF <<- rbind(LD_DF, ld_sub_df)
    
    # Currently unused
    # compParm(dr_model_site, "e", "/")
  }
}

calculate_cumulative_mortality <- function(processed_data, processed_temps, treatment_to_generate = NULL) {
  # Find first and last recorded date of dataset
  first_date <- min(processed_data$date)
  last_date <- max(processed_data$date)
  
  all_sites <- sort(unique(processed_data$site))
  all_dates <- seq(from = first_date, to = last_date, by = "day")
  
  # Calculate total number of individuals per site and treatment
  total_individuals <- processed_data %>%
    group_by(site, treatment) %>%
    summarise(total_ind = n(), .groups = 'drop')
  
  # All Combinations of parameters, if treatment is set, only use set treatment
  all_combinations <- expand.grid(site = unique(processed_data$site),
                                  treatment = unique(ifelse(is.null(treatment_to_generate), processed_data$treatment, treatment_to_generate)),
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
  cumulative_mortality <- full_daily_mortality %>%
    group_by(site, treatment) %>%
    arrange(date) %>%
    left_join(processed_temps, by = "date", suffix = c("", ".new")) %>%
    mutate(
      mittelwert_all_tubs = mittelwert_all_tubs.new,
      cumulative_mortality = cumsum(daily_mortality_rate),
      date_label = paste(format(date, "%d.%m.%Y"), "\n", round(mittelwert_all_tubs, 2), "°C"),
      dose = ifelse(treatment == "mt", 1, 0)
    ) %>%
    return()
}

generate_significance <- function(processed_data, processed_temps, treatment_to_generate = NULL, title_file, site_to_generate = NULL) {
  if (!is.null(site_to_generate)) {
    processed_data <- processed_data %>%
      filter(site == site_to_generate)
  }
  
  if (!is.null(treatment_to_generate)) {
    processed_data <- processed_data %>%
      filter(treatment == treatment_to_generate)
    
    km_fit <- survfit(Surv(time_until_death, event) ~ site, data = processed_data)
  } else {
    km_fit <- survfit(Surv(time_until_death, event) ~ dose, data = processed_data)
  }
  
  g <- ggsurvplot(km_fit,
                  data = processed_data,
                  pval = TRUE,
                  conf.int = TRUE,
                  risk.table = TRUE,
                  ggtheme = GGTHEME) +
    ggtitle(paste("Kaplan-Meier Survival Plot |", title_file,
                  ifelse(is.null(treatment_to_generate), "", paste("(", `TREATMENT_TITLES`[treatment_to_generate], ")")),
                  ifelse(is.null(site_to_generate), "", paste("( Site", site_to_generate, ")"))
    ))
  
  return(g)
}

generate_plot <- function(processed_data, processed_temps, treatment_to_generate, plot_title) {
  
  cumulative_mortality <- calculate_cumulative_mortality(processed_data, processed_temps, treatment_to_generate)
  
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
    GGTHEME
  
  # Return the plot
  return(p)
}

main <- function(files) {
  # Loop over files, process data once, and generate plots for each treatment
  for (file in files) {
    # Get title of current run
    title_file <- FILE_TITLES[file]
    
    # Process the data once per file
    raw_data <- process_data(file, ifelse(file == files[3], "Tabelle1 (2)", "Tabelle1"))
    processed_temps <- process_temps(raw_data)
    
    # Filter out entries, that are not an individual (those are only used to track temperature)
    filtered_data <- raw_data %>%
      filter(!is.na(individual))
    
    calculate_lc50(filtered_data, processed_temps, title_file)
    
    for (site in unique(filtered_data$site)) {
      
      # Generate significance of a specific site, for comparing the treatments
      significance_site <- generate_significance(filtered_data, processed_temps, title_file = title_file, site_to_generate = site)
      
      handle_plot(significance_site$plot, filename = paste0("significance_", tolower(title_file), "_site-", site, ".png"))
    }
    
    # Get the plot that shows the statistical significance for a run, comparing between the two different treatments
    significance <- generate_significance(filtered_data, processed_temps, title_file = title_file)
    
    handle_plot(significance$plot, filename = paste0("significance_", tolower(title_file), ".png"))
    
    # Generate plots for each treatment using the processed data
    for (treatment in TREATMENTS) {
      
      # Get the plot-title for the current treatment
      title_treatment <- TREATMENT_TITLES[treatment]
      # Insert treatment text into plot_title surrounded by ( and )
      plot_title <- paste(title_file, "(", title_treatment, ")")
      
      # Generate plot of given treatment
      plot <- generate_plot(filtered_data, processed_temps, treatment_to_generate = treatment, plot_title = plot_title)
      
      handle_plot(plot, filename = paste0("plot_", tolower(title_file), "_", treatment, ".png"))
    }
  }
  
  save_lc50(LD_DF, LC_PATH)
}

main(FILES)