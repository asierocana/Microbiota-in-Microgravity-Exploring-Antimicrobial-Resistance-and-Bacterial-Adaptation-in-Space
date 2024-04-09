library(ggVennDiagram)
library(dplyr)

# RESISTENCIAS ANTIMICROBIANAS
# Cargar datos abricate del grupo_F y grupo_GC
# Cargar los resultados para cada grupo
setwd("C:/Users/asier/OneDrive/Escritorio/kraken_biom_R//")
grupo_F_abricate <- read.csv("grupo_F_abricate.csv")
grupo_GC_abricate <- read.csv("grupo_GC_abricate.csv")
head(grupo_F_abricate)

# Cargar los datos con delimitador de tabulador
grupo_F_abricate <- read.table("grupo_F_abricate.csv", sep = "\t", header = TRUE, quote = "")

# Visualizar la estructura de los datos
str(grupo_F_abricate)
# Realizar un análisis exploratorio
summary(grupo_F_abricate)
colnames(grupo_F_abricate) <- c("Muestra", "Sequence", "Start", "End", "Gene", "Coverage", "Coverage_map", "Gaps", "%Coverage", "%Identity", "Database", "Accession")
View(grupo_F_abricate)

# Cargar los datos con delimitador de tabulador
grupo_GC_abricate <- read.table("grupo_GC_abricate.csv", sep = "\t", header = TRUE, quote = "")
# Visualizar la estructura de los datos
str(grupo_GC_abricate)
# Realizar un análisis exploratorio
summary(grupo_GC_abricate)
colnames(grupo_GC_abricate) <- c("Muestra", "Sequence", "Start", "End", "Gene", "Coverage", "Coverage_map", "Gaps", "%Coverage", "%Identity", "Database", "Accession")
View(grupo_GC_abricate)

# Realizar un análisis exploratorio
summary(grupo_F_abricate)
summary(grupo_GC_abricate)

# Calcular la frecuencia de cada gen de resistencia en cada grupo
freq_grupo_F <- table(grupo_F_abricate$Gene)
orden_grupo_F <- sort(freq_grupo_F, decreasing = TRUE)
View(orden_grupo_F)
freq_grupo_GC <- table(grupo_GC_abricate$Gene)
orden_grupo_GC <- sort(freq_grupo_GC, decreasing = TRUE)
orden_grupo_GC

tabla_F_data<- as.data.frame(orden_grupo_F)
tabla_F_data
tabla_GC_data<- as.data.frame(orden_grupo_GC)
tabla_GC_data

datos_combinados <- bind_rows(
  mutate(tabla_F_data, Grupo = "Grupo_F"),
  mutate(tabla_GC_data, Grupo = "Grupo_GC")
)

View(datos_combinados)

library(tidyr)

# Usar spread para separar los grupos en dos columnas
datos_separados <- spread(datos_combinados, key = Grupo, value = Freq)

# Visualizar los datos con los grupos separados en columnas
View(datos_separados)

# Reemplazar NA por 0 en el conjunto de datos
datos_separados <- replace(datos_separados, is.na(datos_separados), 0)

# Visualizar el conjunto de datos actualizado
View(datos_separados)

library(dplyr)

# Agregar columnas para indicar presencia/ausencia de genes para cada grupo
datos_separados <- datos_separados %>%
  mutate(presencia_genes_F = ifelse(Grupo_F > 0, 1, 0),
         presencia_genes_GC = ifelse(Grupo_GC > 0, 1, 0))

# Visualizar el conjunto de datos actualizado
View(datos_separados)
# Cambiar el nombre de la columna Var1 por Genes
datos_separados <- datos_separados %>%
  rename(Genes = Var1)

# Visualizar el conjunto de datos actualizado
View(datos_separados)



df <- read.delim("Genes_R_por_grupos")
vector_F <- Genes_R_por_grupos %>% filter(presencia_genes_F == 1) %>% pull(Genes)
vector_GC <- Genes_R_por_grupos %>% filter(presencia_genes_GC == 1) %>% pull(Genes)


lista <- list(Casos = vector_F, Controles = vector_GC)

ggVennDiagram(lista)

no_compartidos <- Genes_R_por_grupos %>% filter(presencia_genes_F == 1 & presencia_genes_GC == 0)
View(no_compartidos)

ggVennDiagram(lista, color = "black", lwd = 0.8, lty = 1) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  ggtitle("Diagrama de Venn con los genes en común y no compartidos")







