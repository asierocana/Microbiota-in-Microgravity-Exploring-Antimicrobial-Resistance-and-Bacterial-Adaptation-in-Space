library("phyloseq")
library(pheatmap)
library(tidyr)
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("biomformat")
library("ampvis2")
library("BiocManager")
library(vegan)
library(ggbiplot)
library(dplyr)
library(stats)

# 27/02. PREPARACIÓN DE LOS DATOS
# Importar datos kraken
setwd("C:/Users/asier/OneDrive/Escritorio/kraken_biom_R//")
muestras_biom <- import_biom ("table_por_muestra.biom")
View(muestras_biom@otu_table)

# Ver la estructura del objeto phyloseq
str(muestras_biom)
class(muestras_biom)
View(muestras_biom@tax_table@.Data)
muestras_biom@tax_table@.Data <- substring(muestras_biom@tax_table@.Data, 4)
colnames(muestras_biom@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
View(muestras_biom@tax_table)

# Realizar un gráfico de barras para la familia
muestras_biom_sin_huecos <- subset_taxa(muestras_biom, Kingdom != "" & Phylum != "" & Class != "" & Order != "" & Family != "" & Genus != "" & Species != "")
View(muestras_biom_sin_huecos@tax_table)

# Acceder a la tabla de OTU
otu_table <- otu_table(muestras_biom)
View(otu_table)
otu_table_data <- as.data.frame(otu_table)
# Acceder a los datos de muestra
sample_data <- sample_data(muestras_biom)
View(sample_data)
# Acceder a la tabla de taxonomía
tax_table <- tax_table(muestras_biom)
View(tax_table)

# Familias más abundantes
familias_abundantes_limpias <- names(sort(rowSums(otu_table(muestras_biom_sin_huecos)), decreasing = TRUE)[1:20])
familias_abundantes_limpias
muestras_filtradas_limpias <- prune_taxa(familias_abundantes_limpias, muestras_biom_sin_huecos)
muestras_filtradas_limpias
View(muestras_filtradas_limpias@sam_data)

head(muestras_filtradas_limpias@otu_table)
head(muestras_filtradas_limpias@tax_table)


# Sumar las abundancias de cada familia en todas las muestras
sumas_familias_grupo_F <- rowSums(muestras_filtradas_limpias@otu_table[, 1:10])
sumas_familias_grupo_F
# Crear un data frame con la información
data_familias_grupo_F <- data.frame(muestras_filtradas_limpias@tax_table, Abundancia = sumas_familias_grupo_F)
View(data_familias_grupo_F)
# Ordenar el data frame por abundancia de mayor a menor
data_familias_ordenado_grupo_F <- data_familias_grupo_F[order(data_familias_grupo_F$Abundancia, decreasing = TRUE), ]
# Visualizar las familias más abundantes
View(data_familias_ordenado_grupo_F)

ggplot(data_familias_ordenado_grupo_F, aes(x = Family, y = Abundancia)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Abundancia de las familias más abundantes en muestras del grupo_F",
       x = "Familia",
       y = "Abundancia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas del eje x para mejor legibilidad


# Sumar las abundancias de cada familia en todas las muestras
sumas_familias_grupo_GC <- rowSums(muestras_filtradas_limpias@otu_table[, 10:20])
sumas_familias_grupo_GC
# Crear un data frame con la información
data_familias_grupo_GC <- data.frame(muestras_filtradas_limpias@tax_table, Abundancia = sumas_familias_grupo_GC)
View(data_familias_grupo_GC)
# Ordenar el data frame por abundancia de mayor a menor
data_familias_ordenado_grupo_GC <- data_familias_grupo_GC[order(data_familias_grupo_GC$Abundancia, decreasing = TRUE), ]
# Visualizar las familias más abundantes
View(data_familias_ordenado_grupo_GC)

ggplot(data_familias_ordenado_grupo_GC, aes(x = Family, y = Abundancia)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Abundancia de las familias más abundantes en muestras del grupo_F",
       x = "Familia",
       y = "Abundancia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotar etiquetas del eje x para mejor legibilidad


# FAMILIAS MÁS ABUNDANTES
# Añadir una columna "Grupo" a los datos
data_familias_ordenado_grupo_F$Grupo <- "Grupo_F"
data_familias_ordenado_grupo_F
data_familias_ordenado_grupo_GC$Grupo <- "Grupo_GC"

# Combinar los datos de ambos grupos
data_combinada <- rbind(data_familias_ordenado_grupo_F, data_familias_ordenado_grupo_GC)
data_combinada
# Crear la gráfica combinada
plot_combinado <- ggplot(data_combinada, aes(x = Family, y = Abundancia, fill = Grupo)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Familias más abundantes en las muestras de los grupos F y GC",
       x = "Familia",
       y = "Abundancia") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar etiquetas del eje x para mejor legibilidad
  scale_fill_manual(values = c("Grupo_F" = "darkblue", "Grupo_GC" = "orange")) +  # Colores personalizados para los grupos
  guides(fill = guide_legend(title = "Grupo"))  # Añadir leyenda

# Visualizar la gráfica combinada
print(plot_combinado)


# A. DIVERSIDAD ALFA.
# Extraer las matrices de abundancias de los grupos F y GC
View(muestras_filtradas_limpias@otu_table)
matriz_abundancias_grupo_F <- as.matrix(muestras_filtradas_limpias@otu_table[, 1:10])
View(matriz_abundancias_grupo_F)
matriz_F_invertida <- t(matriz_abundancias_grupo_F)
matriz_F_invertida
matriz_abundancias_grupo_GC <- as.matrix(muestras_filtradas_limpias@otu_table[, 11:20])
matriz_GC_invertida <- t(matriz_abundancias_grupo_GC)
matriz_GC_invertida

# Análisis de diversidad alpha. SHANNON (POR MUESTRA), TEST WILCOXON
summary(muestras_biom@tax_table@.Data== "")
muestras_biom_sin_huecos <- subset_taxa(muestras_biom, Kingdom != "" & Phylum != "" & Class != "" & Order != "" & Family != "" & Genus != "" & Species != "")
summary(muestras_biom_sin_huecos@tax_table@.Data== "")
head(muestras_biom_sin_huecos@otu_table@.Data)
percentages <- transform_sample_counts(muestras_biom_sin_huecos, function(x) x*100 / sum(x) )
head(percentages@otu_table@.Data)
View(percentages@otu_table)


# Extraer la tabla de abundancias
abundancias <- as.matrix(otu_table(muestras_biom_sin_huecos))
# Calcular la diversidad de Shannon
shannon_diversity <- diversity(abundancias, index = "shannon")
shannon_diversity
# Agregar la diversidad de Shannon al objeto phyloseq
physeq_with_shannon <- merge_phyloseq(muestras_biom_sin_huecos, shannon_diversity)
# Ver las primeras filas del resultado
head(otu_table(physeq_with_shannon))
# Plot de la diversidad de Shannon
plot_richness(muestras_biom_sin_huecos, measures = "Shannon") +
  ggtitle("Diversidad alfa. SHANNON")

# Calcular el índice de Shannon para cada grupo
diversidad_alfa_grupo_F <- diversity(matriz_F_invertida, index = "shannon")
diversidad_alfa_grupo_F
diversidad_alfa_grupo_GC <- diversity(matriz_GC_invertida, index = "shannon")
diversidad_alfa_grupo_GC

grupo <- c(rep("Grupo F", length(diversidad_alfa_grupo_F)), rep("Grupo GC", length(diversidad_alfa_grupo_GC)))
shannon <- c(diversidad_alfa_grupo_F, diversidad_alfa_grupo_GC)
df <- data.frame(grupo, shannon)

# Grafica de los valores de diversidad de Shannon para cada grupo
ggplot(df, aes(x = grupo, y = shannon, fill = grupo)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5, alpha = 0.8) +
  labs(x = "Grupo", y = "Diversidad de Shannon", title = "Diversidad de Shannon por Grupo") +
  theme_minimal()



# USO DE AMPVIS2 (POR MUESTRA)
# Extraer la tabla de OTUs 
muestras_biom_sin_huecos@otu_table
muestras_biom_grupo <- muestras_biom_sin_huecos@sam_data
muestras_biom_grupo
# Crear un vector de grupos
grupos <- c(rep("Grupo_F", 10), rep("Grupo_GC", 10))
grupos
muestras_biom_grupo$Grupo <- grupos
View(muestras_biom_grupo)

# TENGO QUE SEPARAR LOS DATOS POR LOS DOS GRUPOS, JUNTARLOS Y LUEGO HACER EL DATOS_amp
grupo_F_amp <- muestras_biom_sin_huecos@otu_table [, 0:10]
View(grupo_F_amp)
grupo_F_amp <- rowSums(grupo_F_amp)
grupo_F_amp <- as.data.frame(grupo_F_amp)

grupo_GC_amp <- muestras_biom_sin_huecos@otu_table [, 11:20]
View(grupo_GC_amp)
grupo_GC_amp <- rowSums(grupo_GC_amp)
grupo_GC_amp <- as.data.frame(grupo_GC_amp)

grupos_amp <- cbind(grupo_F_amp, grupo_GC_amp)
View(grupos_amp)

datos_amp <- amp_load(
  otutable = grupos_amp,
  taxonomy = muestras_biom_sin_huecos@tax_table
)
datos_amp

# HEATMAP
amp_heatmap(datos_amp) +
  ggtitle("Heatmap de los phylum más abundantes en cada grupo") +
  theme(plot.title = element_text(size = 15))

# CURVA DE RAREFACCIÓN
amp_rarefaction_curve(
  datos_amp,
  stepsize = 1000,
  color_by = NULL,
  facet_by = NULL,
  facet_scales = "fixed"
) +
  ggtitle("Curva de Rarefacción observada")

# BOXPLOT
amp_boxplot(datos_amp) +
  ggtitle("BoxPlot de los géneros bacterianos más abundantes") +
  theme(plot.title = element_text(size = 20))






# B. PRUEBAS ESTADÍSTICAS. WILCOXON Y T TEST
# T TEST
# Realiza el t-test
head(grupo_F_amp)
head(grupo_GC_amp)
t_test_result <- t.test(grupo_F_amp, grupo_GC_amp)

# Mostrar los resultados
print(t_test_result)

# No hay una diferencia significativa entre las medias de los grupos: 
# El p-valor obtenido (0.2205) es mayor que 0.05. No hay suficiente 
# evidencia para afirmar que las medias de los 2 grupos son diferentes.
# Intervalo de confianza para la diferencia en las medias: 
# intervalo confianza 95% para la diferencia en las medias desde -2.076093 
# a 8.994948. Con nivel confianza 95%, es probable que la verdadera 
# diferencia en las medias esté en este rango.

# Convertir las columnas a datos numéricos
grupo_F <- as.numeric(unlist(grupo_F_amp))
head(grupo_F)
grupo_GC <- as.numeric(unlist(grupo_GC_amp))

# Realizar el test de Wilcoxon
wilcox_result <- wilcox.test(grupo_F, grupo_GC)
print(wilcox_result)

# Diferencia significativa en la distribución de los datos entre grupos: 
# p-valor muy pequeño (0.00000000373) 
# Rechazo hipótesis nula: p-valor muy bajo, rechaza hipótesis nula 
# de que no hay diferencia en la distribución entre los grupos. 
# Hay diferencia real y significativa entre las medianas de los dos grupos.
# Esta diferencia de los grupos es lo suficientemente grande, lo que sugiere
# que esta diferencia no es el resultado de la variabilidad aleatoria en los datos.



# DIVERISAD BETA
# Convierte el objeto phyloseq en una matriz de datos
matriz_abundancia <- as.matrix(muestras_biom_sin_huecos@otu_table)
print(matriz_abundancia)
comunidad_beta<- as.data.frame(matriz_abundancia)
comunidad_beta

head(grupo_F_amp)
head(grupo_GC_amp)

data_for_vegdist <- cbind(grupo_F_amp, grupo_GC_amp)
head(data_for_vegdist)
# Calcular la matriz de disimilaridad usando la distancia de Bray-Curtis
disimilarity_matrix <- vegdist(data_for_vegdist, method = "bray")
head(disimilarity_matrix)

# NMDS BRAY
set.seed(0)
nmds_result <- metaMDS(disimilarity_matrix, k = 2)
plot(nmds_result)

# NMDS BRAY BIEN 
data_for_vegdist_invertido <- t (data_for_vegdist)
View(data_for_vegdist_invertido)
nmds_result_bray_prueba <- metaMDS(data_for_vegdist_invertido)
plot(nmds_result_bray_prueba)
plot(nmds_result_bray_prueba$points)

sample_data
ID<- muestras_biom_sin_huecos@sam_data
grupos <- c(rep("grupo_F", 10), rep("grupo_GC", 10))
ID$grupo <- grupos
print(ID)
grupos <- factor(grupos)
table(grupos)
# Crear dos data frames separados por grupo
data_frame_F <- ID[ID$grupo == "grupo_F", ]
data_frame_F
data_frame_GC <- ID[ID$grupo == "grupo_GC", ]
data_frame_GC

# Crear un data frame combinando las coordenadas NMDS y los grupos
plot_data <- data.frame(
  MDS1 = nmds_result_bray_prueba$points[, "MDS1"],
  MDS2 = nmds_result_bray_prueba$points[, "MDS2"],
  grupo = ID$grupo
)
plot_data
# Crear el gráfico de dispersión
ggplot(plot_data, aes(x = MDS1, y = MDS2, color = grupo)) +
  geom_point() +
  labs(title = "NMDS Plot",
       x = "NMDS1",
       y = "NMDS2",
       color = "Grupo") +
  stat_ellipse(geom = "polygon", level = 0.95, alpha=0.2)





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

# Visualizar la estructura de los datos
str(grupo_F_abricate)
str(grupo_GC_abricate)

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


# Filtrar el grupo_F_abricate para incluir solo las muestras con genes de resistencia
grupo_F_con_resistencia <- grupo_F_abricate[grupo_F_abricate$Gene != "", ]
grupo_F_con_resistencia
# Filtrar el grupo_GC_abricate de manera similar
grupo_GC_con_resistencia <- grupo_GC_abricate[grupo_GC_abricate$Gene != "", ]

# Verificar las longitudes de los vectores
length(grupo_F_con_resistencia$Gene)
length(grupo_GC_con_resistencia$Gene)
# Como tienen diferentes longitudes, se descarta la posibilidad de realizar un Chi cuadrado
str(grupo_F_con_resistencia$Gene)
str(grupo_GC_con_resistencia$Gene)

View(grupo_F_con_resistencia)

# GENES APARICIÓN GRUPO F
grupo_F_con_resistencia_genes <- grupo_F_con_resistencia %>%
  select(starts_with("Muestra"), starts_with("Gene"))
# Conteo de cada uno de los genes
conteo_genes_F <- grupo_F_con_resistencia_genes %>%
  count(Muestra, Gene) %>%
  arrange(Muestra, Gene)
# Visualizar el conteo
print(conteo_genes_F)
# Pivoteo de los datos para tener los genes en la primera fila
conteo_genes_F_reshaped <- conteo_genes_F %>%
  pivot_wider(names_from = Gene, values_from = n, values_fill = 0)
View(conteo_genes_F_reshaped)

# GENES APARICIÓN GRUPO GC
grupo_GC_con_resistencia_genes <- grupo_GC_con_resistencia %>%
  select(starts_with("Muestra"), starts_with("Gene"))
# Conteo de cada uno de los genes
conteo_genes_GC <- grupo_GC_con_resistencia_genes %>%
  count(Muestra, Gene) %>%
  arrange(Muestra, Gene)
# Visualizar el conteo
print(conteo_genes_GC)
# Pivoteo de los datos para tener los genes en la primera fila
conteo_genes_GC_reshaped <- conteo_genes_GC %>%
  pivot_wider(names_from = Gene, values_from = n, values_fill = 0)
View(conteo_genes_GC_reshaped)

# GENES APARICIÓN TOTAL
# Combinar los datos anteriores
head(conteo_genes_F_reshaped)
head(conteo_genes_F_reshaped)
library(dplyr)
# Suponiendo que tienes los DataFrames conteo_genes_F_reshaped y conteo_genes_GC_reshaped
conteo_genes_combinado <- bind_rows(
  mutate(conteo_genes_F_reshaped, Grupo = "F"),
  mutate(conteo_genes_GC_reshaped, Grupo = "GC")
)
View(conteo_genes_combinado)
conteo_genes_combinado[is.na(conteo_genes_combinado)] <- 0
View(conteo_genes_combinado)

# Guardar el resultado 
write.csv(conteo_genes_combinado, file = "C:/Users/asier/OneDrive/Escritorio/kraken_biom_R/conteo_genes_combinado.csv", row.names = FALSE)

# Sumar el número de veces que aparece cada gen por columna GRUPO_F
suma_por_gen_F <- colSums(conteo_genes_F_reshaped[, -1])
# Mostrar los resultados
print(suma_por_gen_F)
# Ordenar los genes de mayor a menor abundancia
orden <- order(-suma_por_gen_F)
# Aplicar el orden a los datos
suma_por_gen_ordenado_F <- suma_por_gen_F[orden]
# Crear un gráfico de barras ordenado con límites de eje Y ajustados
barplot(suma_por_gen_ordenado_F, main = "Frecuencia de Genes de Resistencia en Grupo_F", 
        xlab = "Genes", ylab = "Frecuencia", las = 2, cex.names = 0.5, col = "skyblue",
        ylim = c(0, max(suma_por_gen_ordenado_F) + 2))


# Sumar el número de veces que aparece cada gen por columna GRUPO_F
suma_por_gen_GC <- colSums(conteo_genes_GC_reshaped[, -1])
# Mostrar los resultados
print(suma_por_gen_GC)
# Ordenar los genes de mayor a menor abundancia
orden <- order(-suma_por_gen_GC)
# Aplicar el orden a los datos
suma_por_gen_ordenado_GC <- suma_por_gen_GC[orden]
# Crear un gráfico de barras ordenado
barplot(suma_por_gen_ordenado_GC, main = "Frecuencia de Genes de Resistencia en Grupo_GC", 
        xlab = "Genes", ylab = "Frecuencia", las = 2, cex.names = 0.5, col = "skyblue",
        ylim = c(0, max(suma_por_gen_ordenado_GC) + 2))

# Test de Fisher 
library(stats)

head(conteo_genes_F_reshaped)
head(conteo_genes_GC_reshaped)
View(conteo_genes_F_reshaped)
View(conteo_genes_GC_reshaped)
# Paso 1: Extraer las columnas relevantes
genes_F <- conteo_genes_F_reshaped[, -1]  # Excluir la columna 'Muestra'
genes_GC <- conteo_genes_GC_reshaped[, -1]  # Excluir la columna 'Muestra'
View(genes_F)
View(genes_GC)

# Calcula las sumas marginales por fila (suma de genes para cada muestra)
sums_F <- apply(genes_F, 1, sum)
sums_GC <- apply(genes_GC, 1, sum)
head(sums_F)
head(sums_GC)

# Construir la tabla de contingencia
# Crea una matriz con los conteos de genes para cada muestra en filas
# y dos columnas para los dos conjuntos de datos
tabla_contingencia <- cbind(sums_F, sums_GC)
tabla_contingencia
# Calcular el estadístico de prueba de Fisher
fisher_test_result <- fisher.test(tabla_contingencia, simulate.p.value = TRUE)
print(fisher_test_result)

# El resultado del test indica p-valor de 0.001499, lo que sugiere diferencia significativa
# en la distribución de genes entre los dos grupos (se rechaza la hipótesis nula)


























 