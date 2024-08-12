# install.packages("readxl")
# install.packages("tidyverse")

library(readxl)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("E:/Desktop_F/calculo_tasas_crecimientos")

#  Read a specific sheet
reactors_data <- read_excel("reactores.xlsx")

head(reactors_data)
str(reactors_data)
colnames(reactors_data)
row.names(reactors_data)

# Calculate the average of the replicates for each reactor and time point
reactors_data <- reactors_data %>%
  mutate(across(contains("Growth_reactor1_r"), mean_growth_reactor1 = ~mean(c_across(starts_with("Growth_reactor1_r")), na.rm = TRUE))) %>%
  mutate(across(contains("Growth_reactor2_r"), mean_growth_reactor2 = ~mean(c_across(starts_with("Growth_reactor2_r")), na.rm = TRUE))) %>%
  # Calculates the standard deviation between reactors for growth.
  mutate(sd_growth = ~sd(c(mean_growth_reactor1, mean_growth_reactor2), na.rm = TRUE)) %>%
  # Calculates the average and standard deviation for the other measurements
  mutate(across(starts_with("Butyrate_reactor"), list(mean = ~ mean(., na.rm = TRUE), 
                                                      sd = ~ sd(., na.rm = TRUE)))) %>%
  # ...do this for all the other measurements....
  select(Time = `Time [h]`, starts_with("mean"), starts_with("sd"))

# See the results
print(reactors_data)

#------------------------------------------------------------------------------------------------------------------------------------------------

# Calculate averages, standard deviations and standard errors for growth
reactors_data_summary <- reactors_data %>%
  pivot_longer(cols = starts_with("Growth"), names_to = "measurement", values_to = "value") %>%
  separate(measurement, into = c("type", "reactor", "replicate"), sep = "_", remove = FALSE) %>%
  group_by(`Time [h]`, type) %>%
  summarize(mean = mean(value, na.rm = TRUE), 
            sd = sd(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE) / sqrt(n()),
            .groups = 'drop') %>%
  ungroup()

# See the results
print(reactors_data_summary)
colnames(reactors_data_summary)
str(reactors_data_summary)

# --------------------------------------
# combined reactors:
reactors_data_2 <- read_excel("exponential_reactors.xlsx")

# first exponential phase
start_exponential_phase <- 1.83
end_exponential_phase <- 5.78

# Calculate the natural logarithm of the DOs for each replicate of both reactors
reactors_data_2 <- reactors_data_2 %>%
  mutate(across(starts_with("Reactor"), log))

# Calculate the mean of the logarithmic DOs combining all replicates of both reactors
reactors_data_2 <- reactors_data_2 %>%
  mutate(ln_OD_mean_all = rowMeans(select(., starts_with("Reactor")), na.rm = TRUE))

# Filter the data for the exponential phase
exponential_data_all <- filter(reactors_data_2, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)

# Check if there are enough and non-NA rows in exponential_data_all
print(head(exponential_data_all))
print(sum(!is.na(exponential_data_all$ln_OD_mean_all)))

# Attempt to perform regression if sufficient data is available
if(nrow(exponential_data_all) > 0 && sum(!is.na(exponential_data_all$ln_OD_mean_all)) > 0) {
  model_mean_all <- lm(ln_OD_mean_all ~ Time_h, data = exponential_data_all)
  summary_mean_all <- summary(model_mean_all)
  
  # Extract and display the desired statistics
} else {
  cat("Not enough data or all NA in the specified range to perform the regression.\n")
}

exponential_data_all

# Perform linear regression in the exponential phase for the combined mean of all reactors and replicates
model_mean_all <- lm(ln_OD_mean_all ~ Time_h, data = exponential_data_all)
summary_mean_all <- summary(model_mean_all)

# Obtain model slope (growth rate), R-squared, p-value, and standard errors
growth_rate_all <- coef(summary_mean_all)["Time_h", "Estimate"]
r_squared_all <- summary_mean_all$r.squared
p_value_all <- coef(summary_mean_all)["Time_h", "Pr(>|t|)"]
std_error_all <- coef(summary_mean_all)["Time_h", "Std. Error"]

# Calculate the standard deviation and standard error of the mean
residuals_all <- residuals(model_mean_all)
sd_all <- sd(residuals_all)
sem_all <- sd_all / sqrt(length(residuals_all))

# Print final results
cat("Final combined growth rate:", growth_rate_all, "\n")
cat("R-squared:", r_squared_all, "\n")
cat("p-value:", p_value_all, "\n")
cat("Standard error:", std_error_all, "\n")
cat("Standard deviation:", sd_all, "\n")
cat("Standard error of the mean:", sem_all, "\n")

###### second exponential phase
reactors_data_2 <- read_excel("exponential_reactors.xlsx")

start_exponential_phase <- 17.73
end_exponential_phase <- 24.23

# Calculate the natural logarithm of the DOs for each replicate of both reactors
reactors_data_2 <- reactors_data_2 %>%
  mutate(across(starts_with("Reactor"), log))

# Calculate the mean of the logarithmic DOs combining all replicates of both reactors
reactors_data_2 <- reactors_data_2 %>%
  mutate(ln_OD_mean_all = rowMeans(select(., starts_with("Reactor")), na.rm = TRUE))

# Filter the data for the exponential phase
exponential_data_all <- filter(reactors_data_2, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)

# Check if there are enough and non-NA rows in exponential_data_all
print(head(exponential_data_all))
print(sum(!is.na(exponential_data_all$ln_OD_mean_all)))

# Attempt to perform regression if sufficient data is available
if(nrow(exponential_data_all) > 0 && sum(!is.na(exponential_data_all$ln_OD_mean_all)) > 0) {
  model_mean_all <- lm(ln_OD_mean_all ~ Time_h, data = exponential_data_all)
  summary_mean_all <- summary(model_mean_all)
  
  # Extract and display the desired statistics
} else {
  cat("No hay datos suficientes o son todos NA en el rango especificado para realizar la regresión.\n")
}

exponential_data_all

# Perform linear regression in the exponential phase for the combined mean of all reactors and replicates
model_mean_all <- lm(ln_OD_mean_all ~ Time_h, data = exponential_data_all)
summary_mean_all <- summary(model_mean_all)

# Obtain model slope (growth rate), R-squared, p-value, and standard errors
growth_rate_all <- coef(summary_mean_all)["Time_h", "Estimate"]
r_squared_all <- summary_mean_all$r.squared
p_value_all <- coef(summary_mean_all)["Time_h", "Pr(>|t|)"]
std_error_all <- coef(summary_mean_all)["Time_h", "Std. Error"]

# Calculate the standard deviation and standard error of the mean
residuals_all <- residuals(model_mean_all)
sd_all <- sd(residuals_all)
sem_all <- sd_all / sqrt(length(residuals_all))

# Print final results
cat("Final combined growth rate:", growth_rate_all, "\n")
cat("R-squared:", r_squared_all, "\n")
cat("p-value:", p_value_all, "\n")
cat("Standard error:", std_error_all, "\n")
cat("Standard deviation:", sd_all, "\n")
cat("Standard error of the mean:", sem_all, "\n")
