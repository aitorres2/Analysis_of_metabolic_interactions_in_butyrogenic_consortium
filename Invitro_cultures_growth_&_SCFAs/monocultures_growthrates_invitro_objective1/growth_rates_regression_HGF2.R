
# install.packages("readxl")
# install.packages("tidyverse")

library(readxl)
library(tidyverse)

setwd("E:/Desktop_F/calculo_tasas_crecimientos")

# Read a specific sheet
HGF2_data <- read_excel("growth_curves_tesis.xlsx", sheet = "HGF2")

head(HGF2_data)

# Calculate the natural logarithm of the DO for each replicate
HGF2_data <- HGF2_data %>%
  mutate(ln_R1 = log(R1),
         ln_R2 = log(R2),
         ln_R3 = log(R3))

head(HGF2_data)

#-----------------------------------------------------------------------------------------------
ggplot(data = HGF2_data, aes(x = Time_h)) +
  geom_line(aes(y = R1, colour = "Replica 1")) +
  geom_line(aes(y = R2, colour = "Replica 2")) +
  geom_line(aes(y = R3, colour = "Replica 3")) +
  theme_minimal() +
  labs(title = "Curvas de Crecimiento Bacteriano",
       x = "Tiempo (horas)",
       y = "OD620 nm",
       colour = "Réplica") +
  scale_colour_manual(values = c("Replica 1" = "blue", "Replica 2" = "red", "Replica 3" = "green"))

# Create dataframe for the chart
HGF2_data_long <- HGF2_data %>%
  gather(key = "Replica", value = "ln_OD", ln_R1, ln_R2, ln_R3)

# Use ggplot2 to create the graph
ggplot(data = HGF2_data_long, aes(x = `Time_h`, y = ln_OD, color = Replica)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "Logaritmo natural de OD vs. Tiempo",
       x = "Tiempo (h)",
       y = "ln(OD)",
       color = "Réplica")
#-----------------------------------------------------------------------------------------------
# phase 1 glucose:
start_exponential_phase <- 4.99994444
end_exponential_phase <- 9.50019444

# Function to perform linear regression and extract the growth rate
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Returns the growth rate
}

# Calculate the growth rate for each replicate
growth_rate_R1 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R3")

# Calculate the average of the growth rates of the three replicates
average_growth_rate <- mean(c(growth_rate_R1, growth_rate_R2, growth_rate_R3))

# Print the average growth rate
print(average_growth_rate)

# Assuming that you have a vector 'growth_rates' containing the growth rates of your replicas
growth_rates <- c(growth_rate_R1, growth_rate_R2, growth_rate_R3)

# Calculate the standard deviation (SD)
sd_growth_rates <- sd(growth_rates)

# Calculate the standard error of the mean (SEM)
sem_growth_rates <- sd_growth_rates / sqrt(length(growth_rates))

# Print the SD and SEM
print(paste("SD:", sd_growth_rates))
print(paste("SEM:", sem_growth_rates))
#------------------------------------------------------------------------------------------------------
# To present results
# First, calculate the average of the logarithmic DOs for each time point
HGF2_data$ln_OD_mean <- rowMeans(HGF2_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Then, perform the linear regression in the exponential phase
exponential_data <- filter(HGF2_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
model_mean <- lm(ln_OD_mean ~ Time_h, data = exponential_data)

# Prints the model summary for the mean of the replicates
summary_mean <- summary(model_mean)
print(summary_mean$r.squared)  # For R-squared
print(summary_mean$coefficients[2,4])  # For the p-value of the slope

# To obtain the slope (growth rate) of the model
slope <- summary_mean$coefficients["Time_h", "Estimate"]
print(slope)

#-----------------------------------------------------------------------------------------------
# phase 2 inulin/FOS:
start_exponential_phase <- 25.5010556
end_exponential_phase <- 33.5015

# Function to perform linear regression and extract the growth rate
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Returns the growth rate
}

# Calculate the growth rate for each replicate
growth_rate_R1 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R3")

# Calculate the average of the growth rates of the three replicates
average_growth_rate <- mean(c(growth_rate_R1, growth_rate_R2, growth_rate_R3))

# Print the average growth rate
print(average_growth_rate)

# Assuming that you have a vector 'growth_rates' containing the growth rates of your replicas
growth_rates <- c(growth_rate_R1, growth_rate_R2, growth_rate_R3)

# Calculate the standard deviation (SD)
sd_growth_rates <- sd(growth_rates)

# Calculate the standard error of the mean (SEM)
sem_growth_rates <- sd_growth_rates / sqrt(length(growth_rates))

# Print the SD and SEM
print(paste("SD:", sd_growth_rates))
print(paste("SEM:", sem_growth_rates))
#------------------------------------------------------------------------------------------------------
# To present results
# First, calculate the average of the logarithmic DOs for each time point
HGF2_data$ln_OD_mean <- rowMeans(HGF2_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Then, perform the linear regression in the exponential phase
exponential_data <- filter(HGF2_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
model_mean <- lm(ln_OD_mean ~ Time_h, data = exponential_data)

# Prints the model summary for the mean of the replicates
summary_mean <- summary(model_mean)
print(summary_mean$r.squared)  # For R-squared
print(summary_mean$coefficients[2,4])  # For the p-value of the slope

# To obtain the slope (growth rate) of the model
slope <- summary_mean$coefficients["Time_h", "Estimate"]
print(slope)

#------------------------------------------------------------------------------------------------------
# OD data provided
OD_values <- c(0.8893,	0.7631,	0.9039)

# Calculate the average of the DOs
average_OD <- mean(OD_values)

# Calculate the standard deviation (SD)
sd_OD <- sd(OD_values)

# Calculate the standard error (SEM)
sem_OD <- sd_OD / sqrt(length(OD_values))

# Print the results
print(paste("OD Final Promedio:", average_OD))
print(paste("SD:", sd_OD))
print(paste("SEM:", sem_OD))




