library(readxl)
library(tidyverse)
library(stringr)
library("pheatmap")
library("grid")
library(data.table)
library("writexl")
library(ggplot2)
library(matrixStats)
library(reshape2)
library(grid)
library(gridExtra)

# loading flux sampling output from matlab
setwd("/media/alexis/hdd2/objetivo_1_tesis_doctoral/matlab_scripts_genomas")
#sampling <- fread("samples_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_200000_modelseed_ids.csv") TESIS ORIGINAL
sampling <- fread("samples_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_200000_modelseed_Asegurando_que_PT33_cresca_en_monocultivo.csv") # USANDO RXNS DE VMH EN MODELOS GAPSEQ PARA INULINA EXTRACELULAR

#sampling_rxn_names <- fread("samples_EcComRXNS_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_modelseed_ids.csv", header = FALSE)
sampling_rxn_names <- fread("samples_EcComRXNS_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_modelseed_Asegurando_que_PT33_cresca_en_monocultivo.csv", header = FALSE)

modelseed_names <- fread("modelseed_compounds.tsv", header = TRUE)
colnames(modelseed_names)
str(modelseed_names)

# Exploring tables to understand merge
colnames(sampling)[1:2]
colnames(sampling_rxn_names)[1:2]
rownames(sampling)[1:2]
rownames(sampling_rxn_names)[1:2]

# Set the row names of sampling using the first column of sampling_rxn_names
row.names(sampling) <- sampling_rxn_names[[1]]
rownames(sampling)[1:2]

# Identify rows containing "EX_" in their name
ex_rows <- grepl("EX_", row.names(sampling))

# Identify rows containing "biomass" in their name
biomass_rows <- grepl("biomass", row.names(sampling))

# Combine both conditions using the logical operator OR (|)
combined_rows <- ex_rows | biomass_rows

# Filter the original "sampling" table to obtain the rows that meet the following conditions
sampling_filtered <- sampling[combined_rows, ]

# Set the row names of the filtered dataset to the original row names that meet the condition
row.names(sampling_filtered) <- row.names(sampling)[combined_rows]

# Print selected rows
rownames(sampling_filtered)

# Convert row names to character strings
rownames(sampling_filtered) <- as.character(rownames(sampling_filtered))
rownames(sampling_filtered)[1:2]
rownames(sampling_filtered)

rm(sampling) # Use Restart R from Session menu to recover some RAM
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Transform modelseed id names to metabolite names

library(stringr)

# Step 1: Ensure that the row names are strings
rownames(sampling_filtered) <- as.character(rownames(sampling_filtered))
modelseed_names$id <- as.character(modelseed_names$id)

# Step 2: Create a function to look up and get the corresponding names in modelseed_names
get_metabolite_name <- function(metabolite_id) {
  index <- match(metabolite_id, modelseed_names$id)
  if (!is.na(index)) {
    return(modelseed_names$name[index])
  } else {
    return(NA) # Returns NA if the identifier is not found in modelseed_names
  }
}

# Step 3: Replace the row names in sampling_filtered with the corresponding names
unique_metabolite_names <- sapply(rownames(sampling_filtered), function(row_name) {
  bacterium_name <- str_extract(row_name, "^[^_]+") # Extract the name of the bacterium (part before the first underscore)
  metabolite_id <- str_extract(row_name, "(?<=cpd)\\d+") # # Extract the identifier only (cpd...)
  
  if (!is.na(metabolite_id)) {
    metabolite_name <- get_metabolite_name(paste0("cpd", metabolite_id))
    if (!is.na(metabolite_name)) {
      return(paste0(bacterium_name, "_", metabolite_name))
    }
  }
  
  return(row_name) # Return the original row name if the identifier is not found
})

# Resolve duplicates and make row names unique
unique_metabolite_names <- make.unique(unique_metabolite_names, sep = "_")

# Replace the row names in sampling_filtered with the corresponding unique names
sampling_filtered$metabolite_name <- unique_metabolite_names
rownames(sampling_filtered) <- sampling_filtered$metabolite_name
sampling_filtered$metabolite_name <- NULL # We remove the temporary column from metabolite_name

# Verification: Shows the first sampling_filtered rows after replacement
rownames(sampling_filtered)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Making distribution plots

# Filtering only to metabolites of interest
# Vector with all search patterns (search_strings)
search_strings <- c("biomass", "H2O", "Cholate", "Glycine", "Phosphate", "L-Glutamate", "D-Glucose", "L-Lysine", "L-Aspartate",
                    "L-Arginine", "L-Glutamine", "Ornithine", "L-Leucine", "Putrescine", "L-Histidine", "L-Proline",
                    "Xylose", "L-Valine", "L-Threonine", "L-Arabinose", "H+", "Uracil", "Guanosine", "Citrate",
                    "L-Alanine", "L-Methionine", "L-Cysteine", "L-Isoleucine", "L-Tyrosine", "L-Phenylalanine",
                    "L-Serine", "NH3", "Pyruvate", "Acetate", "L-Asparagine", "Fumarate", "GABA", "L-Lactate", "Riboflavin",
                    "L-Tryptophan", "Xanthosine", "Deoxyguanosine", "D-Lactate", "Sucrose", "Succinate", "NADP", "EX_Folate",
                    "D-Mannose", "D-Mannitol", "3-dehydrocholate", "CO2", "Acetaldehyde", "Ethanol", "Formaldehyde", "L-Fucose",
                    "D-Glyceraldehyde", "H2S", "Inulin", "Butyrate", "Methanol", "D-Fructose", "D-Ribose", "Lacto-N-biose"
                    )

# Function for escaping special characters in search patterns
escape_special_characters <- function(string) {
  gsub("([\\{\\}\\[\\]\\?\\*\\+\\(\\)\\|\\^\\$\\#\\:\\\\\\.])", "\\\\\\1", string)
}

# Search patterns with escaped special characters
escaped_search_strings <- escape_special_characters(search_strings)

# Filter metabolites that contain any of the search patterns in their name
matching_metabolites <- rownames(sampling_filtered)[grep(paste(escaped_search_strings, collapse = "|"), rownames(sampling_filtered))]

# Print matching metabolites
matching_metabolites

# Filter out metabolites that do NOT contain "EX_" in their name
matching_metabolites_cleaned <- matching_metabolites[!grepl("EX_", matching_metabolites)]
matching_metabolites_cleaned

# Use grepl to find the rows that contain the search strings in their names
row_indices <- grepl(paste(matching_metabolites_cleaned, collapse="|"), rownames(sampling_filtered))

# Subset the filtered data to keep only the rows that contain the search strings in their names
sampling_final <- sampling_filtered[row_indices, , drop = FALSE]

# Set the row names of the filtered data to the original row names
row.names(sampling_final) <- row.names(sampling_filtered)[row_indices]
rownames(sampling_final)[1:2]

rownames(sampling_final)

# Remove the row with the name "EX_acgam[u]tr" # obtenidos al visualizar el anterior "rownames(sampling_final)"
patterns_to_remove <- c("Heme",
                        "Putrescine",
                        "Xylose",
                        "L-Arabinose",
                        "Uracil",
                        "Guanosine",
                        "Citrate",
                        "Xanthosine",
                        "Deoxyguanosine",
                        "Hexanoate",
                        "MTTL",
                        "HCN",
                        "NADP",
                        "Phosphate",
                        "TRHL",
                        "GSH",
                        "indol",
                        "NH3",
                        "CO2",
                        "D-Ribose"
                        )

for(pattern in patterns_to_remove){
  sampling_final <- sampling_final %>%
    filter(!grepl(pattern, rownames(.)))
}

class(sampling_final)
rownames(sampling_final)

# Adjust the names of the biomasses in sampling_final
sampling_final$metabolite <- gsub("Bifidobacterium_animalis_subsp_lactis_PT33bio1_biomass", "Biomass", sampling_final$metabolite)
sampling_final$metabolite <- gsub("Clostridium_innocuum_HFG2bio1_biomass", "Biomass", sampling_final$metabolite)

# Define the function to extract data by metabolite name
extract_data_by_metabolite <- function(data_frame, metabolites) {
  all_data <- data.frame()  # Create an empty data frame to hold the data for all organisms
  
  for (metabolite in metabolites) {
    # Initialize data1 and data2 for each iteration
    data1 <- data2 <- numeric(0)
    
    # Extract data for Biomass specifically
    if (metabolite == "Biomass") {
      data1 <- na.omit(as.numeric(data_frame[rownames(data_frame) %in% c("Bifidobacterium_animalis_subsp_lactis_PT33bio1_biomass"), 2:ncol(data_frame)]))
      data2 <- na.omit(as.numeric(data_frame[rownames(data_frame) %in% c("Clostridium_innocuum_HFG2bio1_biomass"), 2:ncol(data_frame)]))
    } else {
      # Extract data for other metabolites
      bifido_row <- paste("Bifidobacterium", metabolite, sep = "_")
      clostridium_row <- paste("Clostridium", metabolite, sep = "_")
      
      data1 <- na.omit(as.numeric(data_frame[rownames(data_frame) == bifido_row, 2:ncol(data_frame)]))
      data2 <- na.omit(as.numeric(data_frame[rownames(data_frame) == clostridium_row, 2:ncol(data_frame)]))
    }
    
    # Combining data from both agencies in one dataframe
    df <- data.frame(value = c(data1, data2),
                     organism = c(rep("Bifidobacterium", length(data1)), 
                                  rep("Clostridium", length(data2))),
                     metabolite = factor(rep(metabolite, length(data1) + length(data2)), 
                                         levels = metabolites))
    
    all_data <- rbind(all_data, df)
  }
  
  return(all_data)
}

rownames(sampling_final)

# List of metabolites of interest for final distribution plot
metabolites <- c("Inulin", "Butyrate", "L-Aspartate", "L-Alanine", 
                 "L-Glutamate", "L-Lysine", "L-Isoleucine", "L-Phenylalanine",
                 "L-Proline", "L-Valine", "D-Glucose", "D-Fructose",
                 "Acetate", "L-Lactate", "D-Lactate", "Succinate", 
                 "D-Glyceraldehyde", "Acetaldehyde", "Formaldehyde", "3-dehydrocholate")

# Function call with the data "sampling_final" and the list of metabolites
result_data <- extract_data_by_metabolite(sampling_final, metabolites)

# Generate the graph with ggplot
final_plot <- ggplot(result_data, aes(x = value, fill = organism)) +
  geom_histogram(position = "identity", bins = 120, alpha = 0.5) +
  labs(x = "Range of Values [mmol gDW-1 h-1]", y = "Frequency of Occurrence") +
  theme_bw() +
  scale_fill_manual(name = "Species", 
                    values = c("Bifidobacterium" = "#1b9e77", "Clostridium" = "#d95f02"),
                    labels = c("Bifidobacterium" = "Bifidobacterium animalis lactis PT33", 
                               "Clostridium" = "Clostridium innocuum HGF2")) +
  facet_wrap(~metabolite, ncol = 4, scales = "free", strip.position = "top") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"), 
        panel.spacing = unit(1, "lines"),
        axis.text = element_text(size = 10), 
        legend.text = element_text(size = 18), 
        text = element_text(size = 17))

final_plot 

# Save figure
ggsave("sampling_final_plot_HGF2_lactis_modelseed_final_2.pdf", plot = final_plot, width = 20, height = 12, units = "in", dpi = 320, bg = "white")
ggsave("sampling_final_plot_HGF2_lactis_modelseed_final_2.png", plot = final_plot, width = 20, height = 12, units = "in", dpi = 320, bg = "white")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate basic statistics per metabolite from the distribution of flux sampling

# Filter to obtain only Biomass data
biomass_data <- result_data %>% 
  filter(metabolite == "Biomass")

# See the first Biomass records
head(biomass_data)

# Summary statistics for Biomass values
summary(biomass_data$value) # biomass of clostridum: 0.04115

# Same for Butyrate
butyrate_data <- result_data %>% 
  filter(metabolite == "Butyrate")

# See the first records of Butyrate
head(butyrate_data)

# Summary statistics for Butyrate values
summary(butyrate_data$value)

# Exploring the distribution with a boxplot to identify outliers
#boxplot(biomass_data$value, main = "Biomass Values")
boxplot(butyrate_data$value, main = "Butyrate Values")

# Another option is to explore the density of the data
#plot(density(na.omit(biomass_data$value)), main = "Density of Biomass Values")
plot(density(na.omit(butyrate_data$value)), main = "Density of Butyrate Values")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Making a correlation matrix from external metabolite between bifidobacterium and clostridium from flux sampling

# Select only Clostridium and Bifidobacterium related fluxes
nombres_flujos <- rownames(sampling_filtered)
# Create an index based on whether the stream name contains "Clostridium" or "Bifidobacterium"
indice_flujos <- grep("Clostridium|Bifidobacterium", nombres_flujos)

# Filtering the data and keeping the original names of each flux
flujos_bacterias <- sampling_filtered[indice_flujos, ]
rownames(flujos_bacterias) <- nombres_flujos[indice_flujos]

rownames(flujos_bacterias)

# Transpose the data set
flujos_transpuestos <- t(flujos_bacterias)

# Separate Clostridium and Bifidobacterium streams
flujos_clostridium <- flujos_bacterias[, grepl("Clostridium_", colnames(flujos_bacterias))]
flujos_bifidobacterium <- flujos_bacterias[, grepl("Bifidobacterium_", colnames(flujos_bacterias))]

# debuggin
# Make sure that the transaction_flows have column names
colnames(flujos_transpuestos) <- rownames(flujos_bacterias)

# Separate data for Clostridium and Bifidobacterium using the correct names
flujos_clostridium <- flujos_transpuestos[, grepl("Clostridium_", colnames(flujos_transpuestos))]
flujos_bifidobacterium <- flujos_transpuestos[, grepl("Bifidobacterium_", colnames(flujos_transpuestos))]

# Verify that there is data before proceeding
if (ncol(flujos_clostridium) > 0 && ncol(flujos_bifidobacterium) > 0) {
  # Calculate standard deviations to verify the data
  sds_clostridium <- apply(flujos_clostridium, 2, sd)
  sds_bifidobacterium <- apply(flujos_bifidobacterium, 2, sd)
  print(sds_clostridium)
  print(sds_bifidobacterium)
  
  # If there are no standard deviations from zero, calculate the correlation
  if (!any(sds_clostridium == 0) && !any(sds_bifidobacterium == 0)) {
    matriz_correlacion_cruzada <- cor(flujos_clostridium, flujos_bifidobacterium, method = "spearman")
    print(matriz_correlacion_cruzada)
  } else {
    print("There are columns with zero standard deviation. Correlation cannot be calculated.")
  }
} else {
  print("Insufficient data in one or both data sets to calculate the correlation.")
}

# Calculates the standard deviation for the Clostridium and Bifidobacterium assemblies
sds_clostridium <- apply(flujos_clostridium, 2, sd)
sds_bifidobacterium <- apply(flujos_bifidobacterium, 2, sd)

# Filters columns where the standard deviation is non-zero
flujos_clostridium <- flujos_clostridium[, sds_clostridium > 0]
flujos_bifidobacterium <- flujos_bifidobacterium[, sds_bifidobacterium > 0]

# Check again if there is enough data to calculate the correlation
if (ncol(flujos_clostridium) > 0 && ncol(flujos_bifidobacterium) > 0) {
  matriz_correlacion_cruzada <- cor(flujos_clostridium, flujos_bifidobacterium, method = "spearman")
  print(matriz_correlacion_cruzada)
} else {
  print("There are not enough columns with variability to calculate the correlation.")
}

matriz_correlacion_cruzada

# Convert the correlation matrix to long format
matriz_larga <- melt(matriz_correlacion_cruzada)

library(tidyverse)
library(ggplot2)

# Plotting the heatmap
ggplot(data = matriz_larga, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + # Usar geom_tile para crear el heatmap
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 45)) +
  labs(fill = "Correlation") # Cambiar la etiqueta de la leyenda de color

# improve order for better heatmap visualization

# Extract the names of the metabolites of Bifidobacterium and Clostridium
nombres_metabolitos_bifido <- gsub("Bifidobacterium_", "", colnames(matriz_correlacion_cruzada))
nombres_metabolitos_clostridium <- gsub("Clostridium_", "", rownames(matriz_correlacion_cruzada))

# Identify common names
nombres_comunes <- intersect(nombres_metabolitos_bifido, nombres_metabolitos_clostridium)

# Create an order vector for rows and columns based on common names
orden_filas <- match(nombres_comunes, nombres_metabolitos_clostridium)
orden_columnas <- match(nombres_comunes, nombres_metabolitos_bifido)

# Create the new ordered correlation matrix
matriz_correlacion_ordenada <- matriz_correlacion_cruzada[orden_filas, orden_columnas]

# Make sure that row and column names are correct
rownames(matriz_correlacion_ordenada) <- paste("Clostridium", nombres_comunes, sep="_")
colnames(matriz_correlacion_ordenada) <- paste("Bifidobacterium", nombres_comunes, sep="_")

# Convert the correlation matrix to long format
matriz_larga_ordenada <- melt(matriz_correlacion_ordenada)

# Plotting the heatmap
p4 <- ggplot(data = matriz_larga_ordenada, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black", size = 0.1) +  # Define el color y el tamaño de las líneas entre los cuadros
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", name="Spearman Correlation") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", colour = "white"), # Fondo blanco para el gráfico
    panel.background = element_rect(fill = "white", colour = "white"), # Fondo blanco para el panel
    axis.text.x = element_text(angle = 70, vjust = 1, hjust=1, size = 6, color = "black", face = "bold"), # Hacer texto del eje X más grande y negro
    axis.text.y = element_text(angle = 0, vjust = 1, hjust=1, size = 5, color = "black", face = "bold"), # Hacer texto del eje Y más grande y negro
    axis.title.x = element_blank(), # Quitar "Var1"
    axis.title.y = element_blank(), # Quitar "Var2"
    legend.position = "right" # Posición de la leyenda
  ) +
  labs(x = NULL, y = NULL) # Quitar etiquetas de los ejes

p4

ggsave("matriz_correlacion_sampling_spearman_200000_final.png", p4, width = 10, height = 8, dpi = 600)
