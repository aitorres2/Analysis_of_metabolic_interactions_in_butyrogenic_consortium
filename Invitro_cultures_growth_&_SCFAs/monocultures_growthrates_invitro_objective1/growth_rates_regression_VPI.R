
# install.packages("readxl")
# install.packages("tidyverse")

library(readxl)
library(tidyverse)

setwd("E:/Desktop_F/calculo_tasas_crecimientos")


# Leer una hoja específica
VPI_5482_data <- read_excel("growth_curves_tesis.xlsx", sheet = "VPI_5482")


head(VPI_5482_data)

# Calcular el logaritmo natural de la OD para cada réplica
VPI_5482_data <- VPI_5482_data %>%
  mutate(ln_R1 = log(R1),
         ln_R2 = log(R2),
         ln_R3 = log(R3))

head(VPI_5482_data)

#-----------------------------------------------------------------------------------------------
ggplot(data = VPI_5482_data, aes(x = Time_h)) +
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
VPI_5482_data_long <- VPI_5482_data %>%
  gather(key = "Replica", value = "ln_OD", ln_R1, ln_R2, ln_R3)

# Usar ggplot2 para crear el gráfico
ggplot(data = VPI_5482_data_long, aes(x = `Time_h`, y = ln_OD, color = Replica)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "Logaritmo natural de OD vs. Tiempo",
       x = "Tiempo (h)",
       y = "ln(OD)",
       color = "Réplica")
#-----------------------------------------------------------------------------------------------

start_exponential_phase <- 2.99983333
end_exponential_phase <- 8.00008333

# Función para realizar la regresión lineal y extraer la tasa de crecimiento
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Devuelve la tasa de crecimiento
}

# Calcular la tasa de crecimiento para cada réplica
growth_rate_R1 <- get_growth_rate(VPI_5482_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(VPI_5482_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(VPI_5482_data, start_exponential_phase, end_exponential_phase, "ln_R3")

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
VPI_5482_data$ln_OD_mean <- rowMeans(VPI_5482_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Luego, realiza la regresión lineal en la fase exponencial
exponential_data <- filter(VPI_5482_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
model_mean <- lm(ln_OD_mean ~ Time_h, data = exponential_data)

# Imprime el resumen del modelo para la media de las réplicas
summary_mean <- summary(model_mean)
print(summary_mean$r.squared)  # Para el R cuadrado
print(summary_mean$coefficients[2,4])  # Para el p-value de la pendiente

# Para obtener la pendiente (growth rate) del modelo
slope <- summary_mean$coefficients["Time_h", "Estimate"]
print(slope)

#------------------------------------------------------------------------------------------------------
# Obtener las OD finales de cada réplica (asumiendo que están en la fila 95)
OD_final_R1 <- VPI_5482_data$R1[95]
OD_final_R2 <- VPI_5482_data$R2[95]
OD_final_R3 <- VPI_5482_data$R3[95]

# Crear un vector con las OD finales
OD_finals <- c(OD_final_R1, OD_final_R2, OD_final_R3)

# Calcular el promedio de las OD finales
average_OD_final <- mean(OD_finals)

# Calcular la desviación estándar (SD)
sd_OD_final <- sd(OD_finals)

# Calcular el error estándar (SEM)
sem_OD_final <- sd_OD_final / sqrt(length(OD_finals))

# Imprimir los resultados
print(paste("OD Final Promedio:", average_OD_final))
print(paste("SD:", sd_OD_final))
print(paste("SEM:", sem_OD_final))
#------------------------------------------------------------------------------------------------------
# Datos de OD proporcionados
OD_values <- c(0.44684997, 0.46045002, 0.48074999)

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



