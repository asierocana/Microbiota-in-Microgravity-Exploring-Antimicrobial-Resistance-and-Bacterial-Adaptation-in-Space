# glm
# Modelo de regresión LOGÍSTICA
# Crear un dataframe con los datos de frecuencia de genes de resistencia
data_genes_completos <- data.frame(
  row.names = c("GLDS-417_F13-ISST_359604_S18_1_contigs.fasta", "GLDS-417_F14-ISST_359599_S13_2_contigs.fasta",
                "GLDS-417_F15-ISST_359600_S14_3_contigs.fasta", "GLDS-417_F16-ISST_359597_S11_4_contigs.fasta",
                "GLDS-417_F17-ISST_359603_S17_5_contigs.fasta", "GLDS-417_F3-ISST_359605_S19_6_contigs.fasta",   
                "GLDS-417_F4-ISST_359602_S16_7_contigs.fasta", "GLDS-417_F5-ISST_359598_S12_8_contigs.fasta",
                "GLDS-417_F6-ISST_359601_S15_9_contigs.fasta", "GLDS-417_F7-ISST_359606_S20_10_contigs.fasta", 
                "GLDS-417_GC13-ISST_359592_S6_11_contigs.fasta", "GLDS-417_GC14-ISST_359588_S2_12_contigs.fasta", 
                "GLDS-417_GC15-ISST_359587_S1_13_contigs.fasta", "GLDS-417_GC16-ISST_359593_S7_14_contigs.fasta", 
                "GLDS-417_GC17-ISST_359596_S10_15_contigs.fasta", "GLDS-417_GC3-ISST_359595_S9_16_contigs.fasta",  
                "GLDS-417_GC4-ISST_359594_S8_17_contigs.fasta", "GLDS-417_GC5-ISST_359590_S4_18_contigs.fasta",  
                "GLDS-417_GC6-ISST_359591_S5_19_contigs.fasta", "GLDS-417_GC7-ISST_359589_S3_20_contigs.fasta"),
  Grupo = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanG_1 = c(2, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 3, 2, 2, 0, 3, 2, 0, 1),
  VanH_D_1 = c(1, 1, 2, 1, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 1, 0, 0, 0, 1, 1),
  VanH_D_4 = c(1, 1, 0, 0, 4, 1, 0, 2, 2, 0, 0, 3, 3, 0, 1, 3, 0, 3, 2, 3),
  VanR_D_2 = c(1, 1, 1, 0, 2, 2, 1, 0, 0, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0),
  VanR_G_1 = c(3, 2, 1, 0, 4, 2, 0, 1, 1, 2, 1, 3, 4, 2, 1, 1, 6, 4, 1, 4),
  VanR_Pt_3 = c(1, 1, 1, 0, 0, 1, 0, 2, 1, 1, 0, 1, 0, 2, 0, 1, 2, 0, 0, 1),
  VanR_Pt_5 = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
  VanS_D_1 = c(2, 2, 2, 1, 2, 1, 3, 1, 2, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 1),
  VanS_G_1 = c(1, 1, 0, 0, 2, 3, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1),
  VanX_D_1 = c(1, 1, 0, 0, 0, 2, 1, 0, 0, 0, 2, 1, 0, 0, 0, 2, 0, 2, 4, 3),
  VanX_D_4 = c(1, 2, 2, 1, 0, 1, 0, 0, 0, 2, 1, 0, 1, 1, 0, 1, 1, 4, 2, 2),
  VanY_D_1 = c(2, 4, 1, 1, 0, 2, 1, 1, 3, 2, 2, 0, 3, 1, 2, 1, 5, 1, 2, 4),
  aadE_1 = c(1, 3, 2, 1, 2, 2, 3, 3, 1, 1, 2, 3, 2, 2, 0, 1, 1, 2, 1, 1),
  optrA_3 = c(1, 0, 0, 0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1),
  tet_32_1 = c(1, 0, 0, 1, 0, 0, 0, 1, 0, 2, 1, 2, 1, 0, 1, 1, 0, 1, 0, 1),
  tet_32_2 = c(2, 3, 2, 1, 6, 1, 3, 3, 2, 1, 5, 3, 3, 2, 2, 2, 4, 2, 5, 2),
  tet_40_2 = c(1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 0, 2, 0, 2, 2, 0, 2, 1, 0, 2),
  tet44_1 = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  tetO_1 = c(2, 2, 2, 3, 6, 5, 2, 2, 3, 2, 5, 2, 1, 4, 1, 1, 2, 0, 2, 3),
  tetO_2 = c(1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1, 2, 1, 2, 0, 1, 1),
  tetW_1 = c(2, 1, 1, 2, 10, 4, 1, 2, 1, 3, 1, 2, 1, 3, 2, 2, 4, 2, 1, 3),
  VanD_1 = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanD_2 = c(0, 3, 1, 0, 4, 1, 1, 0, 0, 0, 2, 2, 0, 1, 1, 4, 0, 0, 1, 2),
  VanD_3 = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
  VanR_A_1 = c(0, 1, 0, 0, 10, 0, 1, 0, 0, 1, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0),
  VanR_D_1 = c(0, 1, 1, 2, 2, 0, 0, 1, 0, 0, 1, 2, 1, 1, 0, 2, 0, 1, 0, 1),
  VanR_D_4 = c(0, 3, 2, 1, 2, 2, 0, 1, 3, 1, 2, 0, 2, 2, 1, 1, 2, 1, 1, 3),
  VanS_D_3 = c(0, 1, 2, 0, 0, 2, 0, 1, 1, 2, 0, 1, 2, 0, 0, 0, 1, 0, 2, 1),
  VanS_Pt_3 = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanU_G_1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
  VanXY_G_1 = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1),
  VanY_D_2 = c(0, 2, 0, 1, 4, 0, 1, 0, 2, 0, 0, 2, 1, 1, 0, 1, 1, 0, 1, 2),
  VanY_G_1 = c(0, 1, 2, 2, 2, 3, 3, 4, 1, 8, 1, 4, 1, 0, 1, 1, 0, 3, 3, 3),
  carA_1 = c(0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
  tet40_1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0),
  tetO_3 = c(0, 1, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 2, 1, 0, 0, 2),
  tetW_6 = c(0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 1),
  vatB_1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanD_4 = c(0, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 2, 0, 1, 3, 1, 1, 0, 1),
  VanS_D_2 = c(0, 0, 2, 1, 0, 1, 1, 1, 1, 0, 2, 0, 1, 0, 1, 3, 1, 2, 2, 0),
  aph3_VIIa_1 = c(0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanS_Pt_5 = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanT_G_1 = c(0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 1, 1, 0),
  lsaB_1 = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0),
  otrA_1 = c(0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
  tetB_3 = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 1, 1, 0, 1, 0, 0),
  VanX_D_2 = c(0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0),
  VanR_C_3 = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  msrD_2 = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)


modelo_glm_prueba <- glm(Grupo ~ ., family = "binomial", data = data_genes_completos)
summary (modelo_glm_prueba)

prueba_poisson_glm <- glm(Grupo ~., data = data_genes_completos, family = poisson(link = "log"))
summary (prueba_poisson_glm)

library(MASS)
prueba_glm_nb <- glm.nb(Grupo ~ ., data = data_genes_completos)
summary (prueba_glm_nb)

# Calcular las sumas por fila para obtener el total de lecturas por muestra
total_lecturas <- rowSums(data_genes_completos[, 2:ncol(data_genes_completos)])

# Dividir cada valor en las filas por el total de lecturas correspondiente para obtener las abundancias relativas
data_abundancias_relativas <- data_genes_completos[, 2:ncol(data_genes_completos)] / total_lecturas

# Agregar la columna de Grupo al dataframe de abundancias relativas
data_abundancias_relativas$Grupo <- data_genes_completos$Grupo

# Ajustar el modelo glm con distribución binomial usando las abundancias relativas
modelo_glm_abundancias <- glm(Grupo ~ ., family = "binomial", data = data_abundancias_relativas)

# Ver el resumen del modelo
summary(modelo_glm_abundancias)


# Crear un dataframe con los datos proporcionados
data_genes <- data.frame(
  VanG_1 = c(2, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 3, 2, 2, 0, 3, 2, 0, 1),
  VanH_D_1 = c(1, 1, 2, 1, 0, 1, 0, 1, 0, 2, 1, 2, 0, 1, 1, 0, 0, 0, 1, 1),
  VanH_D_4 = c(1, 1, 0, 0, 4, 1, 0, 2, 2, 0, 0, 3, 3, 0, 1, 3, 0, 3, 2, 3),
  VanR_D_2 = c(1, 1, 1, 0, 2, 2, 1, 0, 0, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0),
  VanR_G_1 = c(3, 2, 1, 0, 4, 2, 0, 1, 1, 2, 1, 3, 4, 2, 1, 1, 6, 4, 1, 4),
  VanR_Pt_3 = c(1, 1, 1, 0, 0, 1, 0, 2, 1, 1, 0, 1, 0, 2, 0, 1, 2, 0, 0, 1),
  VanR_Pt_5 = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0),
  VanS_D_1 = c(2, 2, 2, 1, 2, 1, 3, 1, 2, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 1),
  VanS_G_1 = c(1, 1, 0, 0, 2, 3, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1),
  VanX_D_1 = c(1, 1, 0, 0, 0, 2, 1, 0, 0, 0, 2, 1, 0, 0, 0, 2, 0, 2, 4, 3),
  VanX_D_4 = c(1, 2, 2, 1, 0, 1, 0, 0, 0, 2, 1, 0, 1, 1, 0, 1, 1, 4, 2, 2),
  VanY_D_1 = c(2, 4, 1, 1, 0, 2, 1, 1, 3, 2, 2, 0, 3, 1, 2, 1, 5, 1, 2, 4),
  aadE_1 = c(1, 3, 2, 1, 2, 2, 3, 3, 1, 1, 2, 3, 2, 2, 0, 1, 1, 2, 1, 1),
  optrA_3 = c(1, 0, 0, 0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1),
  tet_32_1 = c(1, 0, 0, 1, 0, 0, 0, 1, 0, 2, 1, 2, 1, 0, 1, 1, 0, 1, 0, 1),
  tet_32_2 = c(2, 3, 2, 1, 6, 1, 3, 3, 2, 1, 5, 3, 3, 2, 2, 2, 4, 2, 5, 2),
  tet_40_2 = c(1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 0, 2, 0, 2, 2, 0, 2, 1, 0, 2),
  tet44_1 = c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  tetO_1 = c(2, 2, 2, 3, 6, 5, 2, 2, 3, 2, 5, 2, 1, 4, 1, 1, 2, 0, 2, 3),
  tetO_2 = c(1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 2, 0, 0, 1, 2, 1, 2, 0, 1, 1),
  tetW_1 = c(2, 1, 1, 2, 10, 4, 1, 2, 1, 3, 1, 2, 1, 3, 2, 2, 4, 2, 1, 3),
  VanD_1 = c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanD_2 = c(0, 3, 1, 0, 4, 1, 1, 0, 0, 0, 2, 2, 0, 1, 1, 4, 0, 0, 1, 2),
  VanD_3 = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0),
  VanR_A_1 = c(0, 1, 0, 0, 10, 0, 1, 0, 0, 1, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0),
  VanR_D_1 = c(0, 1, 1, 2, 2, 0, 0, 1, 0, 0, 1, 2, 1, 1, 0, 2, 0, 1, 0, 1),
  VanR_D_4 = c(0, 3, 2, 1, 2, 2, 0, 1, 3, 1, 2, 0, 2, 2, 1, 1, 2, 1, 1, 3),
  VanS_D_3 = c(0, 1, 2, 0, 0, 2, 0, 1, 1, 2, 0, 1, 2, 0, 0, 0, 1, 0, 2, 1),
  VanS_Pt_3 = c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanU_G_1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0),
  VanXY_G_1 = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1),
  VanY_D_2 = c(0, 2, 0, 1, 4, 0, 1, 0, 2, 0, 0, 2, 1, 1, 0, 1, 1, 0, 1, 2),
  VanY_G_1 = c(0, 1, 2, 2, 2, 3, 3, 4, 1, 8, 1, 4, 1, 0, 1, 1, 0, 3, 3, 3),
  carA_1 = c(0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
  tet40_1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0),
  tetO_3 = c(0, 1, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 1, 2, 1, 0, 0, 2),
  tetW_6 = c(0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 1, 1),
  vatB_1 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanD_4 = c(0, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 0, 2, 0, 1, 3, 1, 1, 0, 1),
  VanS_D_2 = c(0, 0, 2, 1, 0, 1, 1, 1, 1, 0, 2, 0, 1, 0, 1, 3, 1, 2, 2, 0),
  aph3_VIIa_1 = c(0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanS_Pt_5 = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  VanT_G_1 = c(0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 1, 1, 0),
  lsaB_1 = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0),
  otrA_1 = c(0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0),
  tetB_3 = c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 1, 1, 0, 1, 0, 0),
  VanX_D_2 = c(0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0),
  VanR_C_3 = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  msrD_2 = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)

# Sumar todas las filas
sum_row <- colSums(data_genes)

# Ordenar por abundancia descendente
sum_row <- sort(sum_row, decreasing = TRUE)

# Obtener los 11 genes más abundantes
top_11_genes <- names(sum_row)[1:11]

# Mostrar los 11 genes más abundantes
print(top_11_genes)


# CON LOS 11 MÁS ABUNDANTES
data_genes_11_abundantes <- data.frame(
  row.names = c("GLDS-417_F13-ISST_359604_S18_1_contigs.fasta", "GLDS-417_F14-ISST_359599_S13_2_contigs.fasta",
                "GLDS-417_F15-ISST_359600_S14_3_contigs.fasta", "GLDS-417_F16-ISST_359597_S11_4_contigs.fasta",
                "GLDS-417_F17-ISST_359603_S17_5_contigs.fasta", "GLDS-417_F3-ISST_359605_S19_6_contigs.fasta",   
                "GLDS-417_F4-ISST_359602_S16_7_contigs.fasta", "GLDS-417_F5-ISST_359598_S12_8_contigs.fasta",
                "GLDS-417_F6-ISST_359601_S15_9_contigs.fasta", "GLDS-417_F7-ISST_359606_S20_10_contigs.fasta", 
                "GLDS-417_GC13-ISST_359592_S6_11_contigs.fasta", "GLDS-417_GC14-ISST_359588_S2_12_contigs.fasta", 
                "GLDS-417_GC15-ISST_359587_S1_13_contigs.fasta", "GLDS-417_GC16-ISST_359593_S7_14_contigs.fasta", 
                "GLDS-417_GC17-ISST_359596_S10_15_contigs.fasta", "GLDS-417_GC3-ISST_359595_S9_16_contigs.fasta",  
                "GLDS-417_GC4-ISST_359594_S8_17_contigs.fasta", "GLDS-417_GC5-ISST_359590_S4_18_contigs.fasta",  
                "GLDS-417_GC6-ISST_359591_S5_19_contigs.fasta", "GLDS-417_GC7-ISST_359589_S3_20_contigs.fasta"),
  Grupo = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  VanG_1 = c(2, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 3, 2, 2, 0, 3, 2, 0, 1),
  VanR_G_1 = c(3, 2, 1, 0, 4, 2, 0, 1, 1, 2, 1, 3, 4, 2, 1, 1, 6, 4, 1, 4),
  VanH_D_4 = c(1, 1, 0, 0, 4, 1, 0, 2, 2, 0, 0, 3, 3, 0, 1, 3, 0, 3, 2, 3),
  VanS_D_1 = c(2, 2, 2, 1, 2, 1, 3, 1, 2, 3, 1, 1, 2, 1, 1, 1, 1, 2, 4, 1),
  VanY_D_1 = c(2, 4, 1, 1, 0, 2, 1, 1, 3, 2, 2, 0, 3, 1, 2, 1, 5, 1, 2, 4),
  aadE_1 = c(1, 3, 2, 1, 2, 2, 3, 3, 1, 1, 2, 3, 2, 2, 0, 1, 1, 2, 1, 1),
  tet_32_2 = c(2, 3, 2, 1, 6, 1, 3, 3, 2, 1, 5, 3, 3, 2, 2, 2, 4, 2, 5, 2),
  tetO_1 = c(2, 2, 2, 3, 6, 5, 2, 2, 3, 2, 5, 2, 1, 4, 1, 1, 2, 0, 2, 3),
  tetW_1 = c(2, 1, 1, 2, 10, 4, 1, 2, 1, 3, 1, 2, 1, 3, 2, 2, 4, 2, 1, 3),
  VanR_D_4 = c(0, 3, 2, 1, 2, 2, 0, 1, 3, 1, 2, 0, 2, 2, 1, 1, 2, 1, 1, 3),
  VanY_G_1 = c(0, 1, 2, 2, 2, 3, 3, 4, 1, 8, 1, 4, 1, 0, 1, 1, 0, 3, 3, 3)
)


modelo_glm_prueba <- glm(Grupo ~ ., family = "binomial", data = data_genes_11_abundantes)
summary (modelo_glm_prueba)

# Poisson
prueba_poisson_glm <- glm(Grupo ~., data = data_genes_11_abundantes, family = poisson(link = "log"))
summary (prueba_poisson_glm)
# Ninguna significativa

# Modelo Regresión negativa binomial
library(MASS)
prueba_glm_nb <- glm.nb(Grupo ~ ., data = data_genes_11_abundantes)
summary (prueba_glm_nb)
# Ninguna significativa



# PROBAR PARA QUE SALGA GLM NORMAL. ME DEJA MÁXIMO 8
data_genes_prueba_abundantes_8 <- data.frame(
  row.names = c("GLDS-417_F13-ISST_359604_S18_1_contigs.fasta", "GLDS-417_F14-ISST_359599_S13_2_contigs.fasta",
                "GLDS-417_F15-ISST_359600_S14_3_contigs.fasta", "GLDS-417_F16-ISST_359597_S11_4_contigs.fasta",
                "GLDS-417_F17-ISST_359603_S17_5_contigs.fasta", "GLDS-417_F3-ISST_359605_S19_6_contigs.fasta",   
                "GLDS-417_F4-ISST_359602_S16_7_contigs.fasta", "GLDS-417_F5-ISST_359598_S12_8_contigs.fasta",
                "GLDS-417_F6-ISST_359601_S15_9_contigs.fasta", "GLDS-417_F7-ISST_359606_S20_10_contigs.fasta", 
                "GLDS-417_GC13-ISST_359592_S6_11_contigs.fasta", "GLDS-417_GC14-ISST_359588_S2_12_contigs.fasta", 
                "GLDS-417_GC15-ISST_359587_S1_13_contigs.fasta", "GLDS-417_GC16-ISST_359593_S7_14_contigs.fasta", 
                "GLDS-417_GC17-ISST_359596_S10_15_contigs.fasta", "GLDS-417_GC3-ISST_359595_S9_16_contigs.fasta",  
                "GLDS-417_GC4-ISST_359594_S8_17_contigs.fasta", "GLDS-417_GC5-ISST_359590_S4_18_contigs.fasta",  
                "GLDS-417_GC6-ISST_359591_S5_19_contigs.fasta", "GLDS-417_GC7-ISST_359589_S3_20_contigs.fasta"),
  Grupo = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  VanG_1 = c(2, 1, 0, 2, 0, 2, 0, 2, 0, 1, 0, 2, 3, 2, 2, 0, 3, 2, 0, 1),
  VanR_G_1 = c(3, 2, 1, 0, 4, 2, 0, 1, 1, 2, 1, 3, 4, 2, 1, 1, 6, 4, 1, 4),
  VanY_D_1 = c(2, 4, 1, 1, 0, 2, 1, 1, 3, 2, 2, 0, 3, 1, 2, 1, 5, 1, 2, 4),
  aadE_1 = c(1, 3, 2, 1, 2, 2, 3, 3, 1, 1, 2, 3, 2, 2, 0, 1, 1, 2, 1, 1),
  tet_32_2 = c(2, 3, 2, 1, 6, 1, 3, 3, 2, 1, 5, 3, 3, 2, 2, 2, 4, 2, 5, 2),
  tetO_1 = c(2, 2, 2, 3, 6, 5, 2, 2, 3, 2, 5, 2, 1, 4, 1, 1, 2, 0, 2, 3),
  tetW_1 = c(2, 1, 1, 2, 10, 4, 1, 2, 1, 3, 1, 2, 1, 3, 2, 2, 4, 2, 1, 3),
  VanY_G_1 = c(0, 1, 2, 2, 2, 3, 3, 4, 1, 8, 1, 4, 1, 0, 1, 1, 0, 3, 3, 3)
)

# Sumar todas las filas
sum_row <- colSums(data_genes)

# Ordenar por abundancia descendente
sum_row <- sort(sum_row, decreasing = TRUE)

# Obtener los 11 genes más abundantes
top_8_genes <- names(sum_row)[1:8]

# Mostrar los 11 genes más abundantes
print(top_8_genes)

modelo_glm_prueba <- glm(Grupo ~ ., family = "binomial", data = data_genes_prueba_abundantes_8)
summary (modelo_glm_prueba)
# Ninguna significativa






