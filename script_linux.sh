#!/bin/bash

# Paso 1. Activar ambiente de conda
conda activate prueba_bacteriano

# Paso 2. Calidad de las secuencias crudas
fastqc *.fastq.gz

# Paso 3. Procesamiento de los datos
fastp -i GLDS-417_metagenomics_F14-ISST_359599_S13_R1_raw.fastq.gz -I GLDS-417_metagenomics_F14-ISST_359599_S13_R2_raw.fastq.gz -o R1_clean.fq.gz -O R2_clean.fq.gz --cut_by_quality3 25 --cut_by_quality5 25 --cut_mean_quality 25 -l 100 --qualified_quality_phred 25 -h out_FastP.html

# Paso 4. Calidad de las secuencias procesadas
fastqc *clean.fq.gz

# Paso 5. Downsampling del 10% de las lecturas (BBTools)
reformat.sh in=R1_clean.fq.gz in2=R2_clean.fq.gz out=R1_clean_10.fq.gz out2=R2_clean_10.fq.gz samplerate=0.1

# Paso 6. Emsamblaje metagenómico
metaspades.py -k 21,33,55 -1 R1_clean_10.fq.gz -2 R2_clean_10.fq.gz -t 1 -m 15 -o out_metaspades_10_1_hilos

# Paso 7. Clasificación taxonómica de secuencias metagenómicas
kraken2 --db /home/asier.ocana/Documentos/TFM/inicio/kraken2_db/kraken2_Standard/ --paired R1_clean_10.fq.gz R2_clean_10.fq.gz --output kraken2_results_10_definitivo_1_hilo.txt --report kraken2_report.txt --threads 1 --gzip-compressed --memory-mapping

# Paso 8. Binning scaffolds
metabat2 -i scaffolds.fasta -o metabat2_scaffold_10_1_hilo

# Paso 9. Binning contigs
metabat2 -i contigs.fasta -o metabat2_contigs_10_1_hilo

# Paso 10. Agrupar contigs en un solo archivo 
cat /home/asier.ocana/Documentos/TFM/inicio/grupos_muestrales/14/7_metabat2/contigs/metabat2_contigs_10_1_hilo.*.fa > contigs_total.fa

# Paso 11. Agrupar scaffolds en un solo archivo
cat /home/asier.ocana/Documentos/TFM/inicio/grupos_muestrales/14/7_metabat2/scaffolds/metabat2_scaffold_10_1_hilo.*.fa > scaffolds_total.fa

# Paso 12. Agrupar contigs y scaffolds
cat scaffolds_total.fa contigs_total.fa > total.fa


# Paso 12. KRAKEN DE LOS CONTIGS
kraken2 --db /home/asier.ocana/Documentos/TFM/inicio/kraken2_db/kraken2_Standard/ contigs_total.fa --output kraken2_contigs_total.txt --threads 1 --report kraken2_report_contigs.txt --memory-mapping

# Paso 13. BRACKEN CONTIGS
bracken -d /home/asier.ocana/Documentos/TFM/inicio/kraken2_db/kraken2_Standard/ -i kraken2_report_contigs.txt -o salida_bracken_contigs.txt -w salida_bracken_contigs_report.txt

# Paso 14. Kraken de contigs y scaffolds. USAR ESTOS ARCHIVOS
kraken2 --db /home/asier.ocana/Documentos/TFM/inicio/kraken2_db/kraken2_Standard/ total.fa --output kraken_total.txt --threads 1 --report kraken2_report_total.txt --memory-mapping

# Paso 15. PRUEBA DE BRACKEN CONTIGS+SCAFFOLDS
bracken -d /home/asier.ocana/Documentos/TFM/inicio/kraken2_db/kraken2_Standard/ -i kraken2_report_total.txt -o salida_bracken_total.txt -w salida_bracken_total_report.txt

# Paso 16. Convertir el bracken de los contigs a formato .biom
awk -F'\t' '{print $1 "\t" $2 "\t" $3}' /home/asier.ocana/Documentos/TFM/inicio/grupos_muestrales/14/9_braken_contigs/salida_bracken_contigs.txt > salida_bracken_contigs_selected.txt

biom convert -i salida_bracken_contigs_selected.txt -o salida_bracken_contigs.biom --table-type="OTU table" --to-json

# Paso 17. Lo mismo pero para contigs y scaffolds para pasar a formato .biom
awk -F'\t' '{print $1 "\t" $2 "\t" $3}' /home/asier.ocana/Documentos/TFM/inicio/grupos_muestrales/14/10_bracken_contigs-scaffolds/salida_bracken_total.txt > salida_bracken_total_selected.txt

biom convert -i salida_bracken_total_selected.txt -o salida_bracken_total.biom --table-type="OTU table" --to-json


# Paso 18. Unir todas las muestras del report de kraken total en una en formato biom
cat /home/asier.ocana/Documentos/TFM/inicio/kraken-biom/kraken_report_biom/GLDS-417_*1.txt GLDS-417_*2.txt GLDS-417_*3.txt GLDS-417_*4.txt GLDS-417_*5.txt GLDS-417_*6.txt GLDS-417_*7.txt GLDS-417_*8.txt GLDS-417_*9.txt GLDS-417_*10.txt  > /home/asier.ocana/Documentos/TFM/inicio/kraken-biom/kraken_report_biom/grupo_F_muestras.txt
cat /home/asier.ocana/Documentos/TFM/inicio/kraken-biom/kraken_report_biom/GLDS-417_*11.txt GLDS-417_*12.txt GLDS-417_*13.txt GLDS-417_*14.txt GLDS-417_*15.txt GLDS-417_*16.txt GLDS-417_*17.txt GLDS-417_*18.txt GLDS-417_*19.txt GLDS-417_*20.txt  > /home/asier.ocana/Documentos/TFM/inicio/kraken-biom/kraken_report_biom/grupo_GC_muestras.txt

kraken-biom grupo_F_muestras.txt grupo_GC_muestras.txt --fmt json

# Paso 19. OJO. AGRUPAR POR MUESTRA
# Cambié el nombre de todos los archivos kraken2_report por su nombre de muestra. USAR ESTE ARCHIVO PARA R (ESTUDIO GENERAL)
kraken-biom GLDS* --fmt json


# Paso 20. Abricate contigs. Uso de base de datos RESFINDER
abricate -db resfinder contigs.fasta > abricate_contigs_10_1_hilo

# Paso 21. Abricate scaffolds.
abricate -db resfinder scaffolds.fasta > abricate_scaffold_10_1_hilo

# Paso 22. Unir resultados del abricate por su grupo correspondiente
cat abricate_contigs_10_1_grupo1-10 > grupo_F.csv 
cat abricate_contigs_10_1_grupo11-20 > grupo_GC.csv 
 
# Paso 23. Unir abricates por muestra. CAMBIAR NOMBRES A LOS CONTIGS POR NOMBRE MUESTRA
abricate -db resfinder GLDS-417_F*_1-10_contigs.fasta > GLDS-417_F*_1-10_abricate
abricate -db resfinder GLDS-417_GC*_11-20_contigs.fasta > GLDS-417_GC*_11-20_abricate

# Paso 24. Agrupar archivos por grupos muestrales con cat. USAR ESTOS ARCHIVOS PARA R (RESISTENCIAS ANTIMICROBIANAS)
cat GLDS-417_F13-ISST_359604_S18_1_abricate GLDS-417_F14-ISST_359599_S13_2_abricate GLDS-417_F15-ISST_359600_S14_3_abricate GLDS-417_F16-ISST_359597_S11_4_abricate GLDS-417_F17-ISST_359603_S17_5_abricate GLDS-417_F17-ISST_359603_S17_5_abricate GLDS-417_F3-ISST_359605_S19_6_abricate GLDS-417_F4-ISST_359602_S16_7_abricate GLDS-417_F5-ISST_359598_S12_8_abricate GLDS-417_F6-ISST_359601_S15_9_abricate GLDS-417_F7-ISST_359606_S20_10_abricate > grupo_F_abricate.csv
cat GLDS-417_GC13-ISST_359592_S6_11_abricate GLDS-417_GC14-ISST_359588_S2_12_abricate GLDS-417_GC15-ISST_359587_S1_13_abricate GLDS-417_GC16-ISST_359593_S7_14_abricate GLDS-417_GC17-ISST_359596_S10_15_abricate GLDS-417_GC3-ISST_359595_S9_16_abricate GLDS-417_GC4-ISST_359594_S8_17_abricate GLDS-417_GC5-ISST_359590_S4_18_abricate GLDS-417_GC6-ISST_359591_S5_19_abricate GLDS-417_GC7-ISST_359589_S3_20_abricate > grupo_GC_abricate.csv

cat grupo_F_abricate.csv grupo_GC_abricate.csv > abricate_total_muestras.csv









