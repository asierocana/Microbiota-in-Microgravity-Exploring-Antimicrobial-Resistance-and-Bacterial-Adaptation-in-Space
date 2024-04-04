# Microbiota-in-Microgravity-Exploring-Antimicrobial-Resistance-and-Bacterial-Adaptation-in-Space
Metagenomic analysis of microbial communities from space and Earth simulations. Employed bioinformatics tools in Linux environment for data processing, taxonomic classification, and antimicrobial resistance identification. it includes a study in RStudio.


# Master's Thesis Project - International University of Valencia
# Microbiota in Microgravity: Exploring Antimicrobial Resistance and Bacterial Adaptation in Space

This project aims to conduct a metagenomic analysis of fecal samples from female mice in space under microgravity conditions and samples on Earth with these extreme conditions simulated. A variety of bioinformatics tools will be employed in a Linux environment. The analysis encompasses several steps, from the quality assessment of raw sequences to the identification of antimicrobial resistances.

# Workflow Description
1. Conda Environment Activation: The Conda environment named "prueba_bacteriano" is activated to ensure the proper execution of the tools used in the analysis.

2. Quality Assessment of Raw Sequences: FastQC is used to assess the quality of raw DNA sequences obtained from the samples.

3. Data Processing: Raw sequences are processed using Fastp to remove adapters and filter out low-quality sequences, generating clean sequences ready for the next step.

4.Reassessment of Quality: FastQC is again employed to evaluate the quality of processed sequences and confirm the effectiveness of the cleaning process.

5. Downsampling of 10% of Reads: The amount of data is reduced by randomly selecting 10% of the processed reads.

6. Metagenomic Assembly: A metagenomic assembly of the reduced sequences is performed using MetaSPAdes to reconstruct possible genomes of microorganisms present in the sample.

7. Taxonomic Classification: Kraken2 is used to perform taxonomic classification of the metagenomic sequences, identifying the microorganisms present in the sample.

8. Binning of Scaffolds and Contigs: MetaBAT2 is employed to group the sequences into bins representing possible individual genomes of microorganisms present in the sample.

9. Grouping of Contigs and Scaffolds: Contigs and scaffolds are grouped into individual files and then merged to create a total file.

10. Analysis of Antimicrobial Resistances: An analysis of antimicrobial resistances is conducted using the RESFINDER database and Abricate.

# Project Execution
To execute the project, follow these steps:

Ensure Conda and the mentioned tools are installed in the working environment.

Execute each step of the workflow sequentially, ensuring to provide the paths of input and output files as necessary.

Verify the results of each step to ensure the quality and accuracy of the analysis.

Perform any additional analysis or data processing as necessary to meet the project objectives.

# References
FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Fastp: https://github.com/OpenGene/fastp
MetaSPAdes: https://github.com/ablab/spades
Kraken2: https://ccb.jhu.edu/software/kraken2/
MetaBAT2: https://bitbucket.org/berkeleylab/metabat
Abricate: https://github.com/tseemann/abricate

# Author
This project was carried out by Asier Ocaña Aguirrezabalaga as part of the Master's Thesis in Bioinformatics conducted at the International University of Valencia in April 2024.
