library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

setwd("/home/alexis/Desktop/RESULTADOS_FINALES_ERGO/kallisto_TRABAJADO/PT33/descargados_ERGO_deseq2/GSEA_gage/3_PT33Con_late_vs_PT33_late")

# Leer los datos desde los archivos CSV
upregulated_ergo <- read_csv("gage (ergo): Consortium_Late - PT8_Late_upregulated.csv")
downregulated_ergo <- read_csv("gage (ergo): Consortium_Late - PT8_Late_downregulated.csv")
colnames(upregulated_ergo)

# Filtrar los datos para mostrar solo elementos significativos
umbral_p_valor <- 0.05
umbral_q_valor <- 0.05
upregulated_filtered <- upregulated_ergo %>%
  filter(`P Value` <= umbral_p_valor, `Q Value` <= umbral_q_valor)

upregulated_ergo_plot <- ggplot(upregulated_filtered, aes(x=reorder(Name, -`Stat Mean`), y=`Stat Mean`, fill=`Q Value`)) +
  geom_bar(stat="identity") +
  coord_flip() + # Las barras horizontalmente
  scale_fill_gradient(low="red", high="blue") + # Los valores Q más bajos en rojo, los altos en azul
  labs(title = "PT33 Consortium Late vs PT33 Single Late - Upregulated genes", x="ERGO database", y="Stat Mean", fill="Q Value") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    text = element_text(face = "bold"), # Hace todo el texto negrita
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Reduce aún más el tamaño y rota las etiquetas del eje x
    axis.text.y = element_text(size = 9), # Reduce aún más el tamaño de las etiquetas del eje y
    axis.title = element_text(size = 12), # Reduce el tamaño del título de los ejes
    plot.title = element_text(size = 10), # Reduce el tamaño del título del gráfico
    legend.title = element_text(size = 14), # Reduce el tamaño del título de la leyenda
    legend.text = element_text(size = 12), # Reduce el tamaño del texto de la leyenda
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"), # Ajusta los márgenes del gráfico
    panel.grid.major.y = element_line(color = "grey", size = 0.5), # Líneas de referencia principales en gris
    panel.grid.minor.y = element_blank(), # Sin líneas de referencia menores
    panel.grid.major.x = element_blank(), # Sin líneas de cuadrícula en el eje x
    panel.grid.minor.x = element_blank() # Sin líneas de cuadrícula menores en el eje x
  )

# Visualizar el gráfico 'upregulated' con las mejoras
upregulated_ergo_plot

# Guardar la imagen con fondo blanco
ggsave("upregulated_ergo_plot.png", upregulated_ergo_plot, width = 20, height = 15, dpi = 600, bg = "white")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
downregulated_filtered <- downregulated_ergo %>%
  filter(`P Value` <= umbral_p_valor, `Q Value` <= umbral_q_valor)

downregulated_ergo_plot <- ggplot(downregulated_filtered, aes(x=reorder(Name, -`Stat Mean`), y=`Stat Mean`, fill=`Q Value`)) +
  geom_bar(stat="identity") +
  coord_flip() + # Las barras horizontalmente
  scale_fill_gradient(low="red", high="blue") + # Los valores Q más bajos en rojo, los altos en azul
  labs(title = "PT33 Consortium Late vs PT33 Single Late - Downregulated genes", x="ERGO database", y="Stat Mean", fill="Q Value") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    text = element_text(face = "bold"), # Hace todo el texto negrita
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), # Reduce aún más el tamaño y rota las etiquetas del eje x
    axis.text.y = element_text(size = 9), # Reduce aún más el tamaño de las etiquetas del eje y
    axis.title = element_text(size = 12), # Reduce el tamaño del título de los ejes
    plot.title = element_text(size = 10), # Reduce el tamaño del título del gráfico
    legend.title = element_text(size = 14), # Reduce el tamaño del título de la leyenda
    legend.text = element_text(size = 12), # Reduce el tamaño del texto de la leyenda
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"), # Ajusta los márgenes del gráfico
    panel.grid.major.y = element_line(color = "grey", size = 0.5), # Líneas de referencia principales en gris
    panel.grid.minor.y = element_blank(), # Sin líneas de referencia menores
    panel.grid.major.x = element_blank(), # Sin líneas de cuadrícula en el eje x
    panel.grid.minor.x = element_blank() # Sin líneas de cuadrícula menores en el eje x
  )

# Visualizar el gráfico 'upregulated' con las mejoras
downregulated_ergo_plot

# Guardar la imagen con fondo blanco
ggsave("downregulated_ergo_plot.png", downregulated_ergo_plot, width = 20, height = 15, dpi = 600, bg = "white")




#-------------------------------------------------------------------------------------------------------------------------------------------------------------v
upregulated_filtered

# Gense set size is not important data to show
# Crear el gráfico completo sin geom_text_repel
dot_plot <- ggplot(upregulated_filtered, aes(x=`Stat Mean`, y=reorder(Name, `Stat Mean`))) +
  geom_point(aes(size=`Set Size`, color=`Q Value`)) +
  scale_size(range = c(3, 10)) + # Ajusta el rango de tamaño de los puntos
  scale_color_gradient(low="red", high="blue") + # Escala de color de Q Value
  labs(
    title = "PT33 Consortium Late vs PT33 Late - Upregulated genes",
    x = "Stat Mean", 
    y = "ERGO database", 
    color = "Q Value", 
    size = "Gene Set Size"
  ) +
  theme_minimal() +
  theme(
    text = element_text(face = "bold"), # Hace todo el texto negrita
    legend.position = "right",
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(color = "grey", linewidth = 0.5), # Líneas de referencia
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Mostrar el gráfico completo
print(dot_plot)

# Save the plot ensuring the device has enough space to handle the dimensions
ggsave("upregulated_filtered_dot_plot_complete.png", plot = dot_plot, width = 14, height = 10, dpi = 600, bg = "white")



