
# install.packages("readxl")
# install.packages("tidyverse")

library(readxl)
library(tidyverse)

setwd("E:/Desktop_F/calculo_tasas_crecimientos")

# Read a specific sheet
M38_data <- read_excel("growth_curves_tesis.xlsx", sheet = "M38")

head(M38_data)

# Calculate the natural logarithm of the DO for each replicate
M38_data <- M38_data %>%
  mutate(ln_R1 = log(R1),
         ln_R2 = log(R2),
         ln_R3 = log(R3))

head(M38_data)

#-----------------------------------------------------------------------------------------------
ggplot(data = M38_data, aes(x = Time_h)) +
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
M38_data_long <- M38_data %>%
  gather(key = "Replica", value = "ln_OD", ln_R1, ln_R2, ln_R3)

# Use ggplot2 to create the graph
ggplot(data = M38_data_long, aes(x = `Time_h`, y = ln_OD, color = Replica)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "Logaritmo natural de OD vs. Tiempo",
       x = "Tiempo (h)",
       y = "ln(OD)",
       color = "Réplica")
#-----------------------------------------------------------------------------------------------

start_exponential_phase <- 4.49991667
end_exponential_phase <- 9.50019444

# Function to perform linear regression and extract the growth rate
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Devuelve la tasa de crecimiento
}

# Calculate the growth rate for each replicate
growth_rate_R1 <- get_growth_rate(M38_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(M38_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(M38_data, start_exponential_phase, end_exponential_phase, "ln_R3")

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
M38_data$ln_OD_mean <- rowMeans(M38_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Then, perform the linear regression in the exponential phase
exponential_data <- filter(M38_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
model_mean <- lm(ln_OD_mean ~ Time_h, data = exponential_data)

# Prints the model summary for the mean of the replicates
summary_mean <- summary(model_mean)
print(summary_mean$r.squared)  # Para el R cuadrado
print(summary_mean$coefficients[2,4])  # Para el p-value de la pendiente

# To obtain the slope (growth rate) of the model
slope <- summary_mean$coefficients["Time_h", "Estimate"]
print(slope)

#------------------------------------------------------------------------------------------------------
# Obtain the final ODs of each replica (assuming they are in row 95)
OD_final_R1 <- M38_data$R1[95]
OD_final_R2 <- M38_data$R2[95]
OD_final_R3 <- M38_data$R3[95]

# Create a vector with the final ODs
OD_finals <- c(OD_final_R1, OD_final_R2, OD_final_R3)

# Calculate the average of the final DOs
average_OD_final <- mean(OD_finals)

# Calculate the standard deviation (SD)
sd_OD_final <- sd(OD_finals)

# Calculate the standard error (SEM)
sem_OD_final <- sd_OD_final / sqrt(length(OD_finals))

# Print the results
print(paste("OD Final Promedio:", average_OD_final))
print(paste("SD:", sd_OD_final))
print(paste("SEM:", sem_OD_final))
#------------------------------------------------------------------------------------------------------
# OD data provided
OD_values <- c(0.90090003, 0.91590002, 0.94439998)

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



