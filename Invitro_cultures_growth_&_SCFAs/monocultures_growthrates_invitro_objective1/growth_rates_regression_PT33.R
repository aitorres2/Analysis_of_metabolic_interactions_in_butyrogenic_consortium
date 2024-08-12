
# install.packages("readxl")
# install.packages("tidyverse")

library(readxl)
library(tidyverse)

setwd("E:/Desktop_F/calculo_tasas_crecimientos")

  
# Leer una hoja específica
PT33_data <- read_excel("growth_curves_tesis.xlsx", sheet = "PT33")
HGF2_data <- read_excel("growth_curves_tesis.xlsx", sheet = "HGF2")
VPI_5482_data <- read_excel("growth_curves_tesis.xlsx", sheet = "VPI_5482")
WAL14673_data <- read_excel("growth_curves_tesis.xlsx", sheet = "WAL14673")
M38_data <- read_excel("growth_curves_tesis.xlsx", sheet = "M38")
M62_data <- read_excel("growth_curves_tesis.xlsx", sheet = "M62")
reactor_HGF2_PT33_VPI_data <- read_excel("growth_curves_tesis.xlsx", sheet = "reactor_HGF2_PT33_VPI")

head(PT33_data)

# Calcular el logaritmo natural de la OD para cada réplica
PT33_data <- PT33_data %>%
  mutate(ln_R1 = log(R1),
         ln_R2 = log(R2),
         ln_R3 = log(R3))

head(PT33_data)

#-----------------------------------------------------------------------------------------------
ggplot(data = PT33_data, aes(x = Time_h)) +
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
PT33_data_long <- PT33_data %>%
  gather(key = "Replica", value = "ln_OD", ln_R1, ln_R2, ln_R3)

# Usar ggplot2 para crear el gráfico
ggplot(data = PT33_data_long, aes(x = `Time_h`, y = ln_OD, color = Replica)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "Logaritmo natural de OD vs. Tiempo",
       x = "Tiempo (h)",
       y = "ln(OD)",
       color = "Réplica")
#-----------------------------------------------------------------------------------------------

start_exponential_phase <- 0.4996944 # 1.9997778
end_exponential_phase <- 6.0000000 # 1.9997778

# Función para realizar la regresión lineal y extraer la tasa de crecimiento
get_growth_rate <- function(data, start_time, end_time, response_var) {
  exponential_data <- filter(data, Time_h >= start_time & Time_h <= end_time)
  model <- lm(reformulate("Time_h", response = response_var), data = exponential_data)
  return(coef(model)["Time_h"]) # Devuelve la tasa de crecimiento
}

# Calcular la tasa de crecimiento para cada réplica
growth_rate_R1 <- get_growth_rate(PT33_data, start_exponential_phase, end_exponential_phase, "ln_R1")
growth_rate_R2 <- get_growth_rate(PT33_data, start_exponential_phase, end_exponential_phase, "ln_R2")
growth_rate_R3 <- get_growth_rate(PT33_data, start_exponential_phase, end_exponential_phase, "ln_R3")

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
PT33_data$ln_OD_mean <- rowMeans(PT33_data[, c("ln_R1", "ln_R2", "ln_R3")])

# Luego, realiza la regresión lineal en la fase exponencial
exponential_data <- filter(PT33_data, Time_h >= start_exponential_phase & Time_h <= end_exponential_phase)
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
OD_final_R1 <- PT33_data$R1[95]
OD_final_R2 <- PT33_data$R2[95]
OD_final_R3 <- PT33_data$R3[95]

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
OD_values <- c(0.3128, 0.3804, 0.3333)

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



