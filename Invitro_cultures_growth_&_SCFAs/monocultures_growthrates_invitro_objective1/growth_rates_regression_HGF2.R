
# install.packages("readxl")
# install.packages("tidyverse")

library(readxl)
library(tidyverse)

setwd("E:/Desktop_F/calculo_tasas_crecimientos")
#E:/THESIS PHD FINAL/calculo_tasas_crecimientos

# Leer una hoja específica
HGF2_data <- read_excel("growth_curves_tesis.xlsx", sheet = "HGF2")
VPI_5482_data <- read_excel("growth_curves_tesis.xlsx", sheet = "VPI_5482")
WAL14673_data <- read_excel("growth_curves_tesis.xlsx", sheet = "WAL14673")
M38_data <- read_excel("growth_curves_tesis.xlsx", sheet = "M38")
M62_data <- read_excel("growth_curves_tesis.xlsx", sheet = "M62")
reactor_HGF2_HGF2_VPI_data <- read_excel("growth_curves_tesis.xlsx", sheet = "reactor_HGF2_HGF2_VPI")

m38_paper_data <- read_excel("growth_curves_tesis.xlsx", sheet = "m38_paper")

head(HGF2_data)

# Calcular el logaritmo natural de la OD para cada réplica
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

# Crear un dataframe largo para el gráfico
HGF2_data_long <- HGF2_data %>%
  gather(key = "Replica", value = "ln_OD", ln_R1, ln_R2, ln_R3)

# Usar ggplot2 para crear el gráfico
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

# Función para realizar la regresión lineal y extraer la tasa de crecimiento
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Devuelve la tasa de crecimiento
}

# Calcular la tasa de crecimiento para cada réplica
growth_rate_R1 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R3")

# Calcular el promedio de las tasas de crecimiento de las tres réplicas
average_growth_rate <- mean(c(growth_rate_R1, growth_rate_R2, growth_rate_R3))

# Imprimir la tasa de crecimiento promedio
print(average_growth_rate)

# Suponiendo que tienes un vector 'growth_rates' que contiene las tasas de crecimiento de tus réplicas
growth_rates <- c(growth_rate_R1, growth_rate_R2, growth_rate_R3)  # Reemplaza con tus valores reales

# Calcular la desviación estándar (SD)
sd_growth_rates <- sd(growth_rates)

# Calcular el error estándar de la media (SEM)
sem_growth_rates <- sd_growth_rates / sqrt(length(growth_rates))

# Imprimir el SD y SEM
print(paste("SD:", sd_growth_rates))
print(paste("SEM:", sem_growth_rates))
#------------------------------------------------------------------------------------------------------
# Para presentar resultados
# Primero, calcula la media de las ODs logarítmicas para cada punto de tiempo
HGF2_data$ln_OD_mean <- rowMeans(HGF2_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Luego, realiza la regresión lineal en la fase exponencial
exponential_data <- filter(HGF2_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
model_mean <- lm(ln_OD_mean ~ Time_h, data = exponential_data)

# Imprime el resumen del modelo para la media de las réplicas
summary_mean <- summary(model_mean)
print(summary_mean$r.squared)  # Para el R cuadrado
print(summary_mean$coefficients[2,4])  # Para el p-value de la pendiente

# Para obtener la pendiente (growth rate) del modelo
slope <- summary_mean$coefficients["Time_h", "Estimate"]
print(slope)

#-----------------------------------------------------------------------------------------------
# phase 2 inulin/FOS:
start_exponential_phase <- 25.5010556
end_exponential_phase <- 33.5015

# Función para realizar la regresión lineal y extraer la tasa de crecimiento
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Devuelve la tasa de crecimiento
}

# Calcular la tasa de crecimiento para cada réplica
growth_rate_R1 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(HGF2_data, start_exponential_phase, end_exponential_phase, "ln_R3")

# Calcular el promedio de las tasas de crecimiento de las tres réplicas
average_growth_rate <- mean(c(growth_rate_R1, growth_rate_R2, growth_rate_R3))

# Imprimir la tasa de crecimiento promedio
print(average_growth_rate)

# Suponiendo que tienes un vector 'growth_rates' que contiene las tasas de crecimiento de tus réplicas
growth_rates <- c(growth_rate_R1, growth_rate_R2, growth_rate_R3)  # Reemplaza con tus valores reales

# Calcular la desviación estándar (SD)
sd_growth_rates <- sd(growth_rates)

# Calcular el error estándar de la media (SEM)
sem_growth_rates <- sd_growth_rates / sqrt(length(growth_rates))

# Imprimir el SD y SEM
print(paste("SD:", sd_growth_rates))
print(paste("SEM:", sem_growth_rates))
#------------------------------------------------------------------------------------------------------
# Para presentar resultados
# Primero, calcula la media de las ODs logarítmicas para cada punto de tiempo
HGF2_data$ln_OD_mean <- rowMeans(HGF2_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Luego, realiza la regresión lineal en la fase exponencial
exponential_data <- filter(HGF2_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
model_mean <- lm(ln_OD_mean ~ Time_h, data = exponential_data)

# Imprime el resumen del modelo para la media de las réplicas
summary_mean <- summary(model_mean)
print(summary_mean$r.squared)  # Para el R cuadrado
print(summary_mean$coefficients[2,4])  # Para el p-value de la pendiente

# Para obtener la pendiente (growth rate) del modelo
slope <- summary_mean$coefficients["Time_h", "Estimate"]
print(slope)

#------------------------------------------------------------------------------------------------------
# Datos de OD proporcionados
OD_values <- c(0.8893,	0.7631,	0.9039)

# Calcular el promedio de las OD
average_OD <- mean(OD_values)

# Calcular la desviación estándar (SD)
sd_OD <- sd(OD_values)

# Calcular el error estándar (SEM)
sem_OD <- sd_OD / sqrt(length(OD_values))

# Imprimir los resultados
print(paste("OD Final Promedio:", average_OD))
print(paste("SD:", sd_OD))
print(paste("SEM:", sem_OD))




