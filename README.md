# Myeloproliferative neoplasm patient cohort data visualization

This web application allows the user to explore interactively data from a cohort of 113 myeloproliferative neoplasm (MPN) blood donors and 15 healthy blood controls.  
To ensure that the dependencies are fulfilled, you can install the packages needed by running the install_dependencies.R script:

	Rscript install_dependencies.R


## Features

### Cohort exploration

Represent and compare graphically several variables describing the cohort.

### Variants exploration

Filter variants to select only the somatic mutations with sufficient evidence. Represent mutated gene per patient and co-occurrence of variants, as well as variants per disease and disease subtypes.

### Aberrations and gene fusions exploration

Represent co-occurences of variants and genomic aberrations, occurences of aberrations per disease.  
Display summary of genetic events per patient and per disease.  
A network, sketching the co-occurences identified for aberrations, mutations and diseases in the other tabs, is also available.


## Changelog

### v0.1

Implement basic cohort description plots

### v0.2

Implement basic variants description plots

### v0.3

Implement full CSS layout

### v0.4

Implement aberration description plots, and patient-wise Circos plots.

### v0.5

Implement disease-wise Circos plots.

### v0.6

Implement co-occurence networks.

## To-do list

### Cohort exploration

* Take into account NAs, non relevant-points and redundancy
* Sort legend on pie charts

### Variants exploration

* Correct color scheme of one by one matrices

### Gene expression exploration

### Differential splicing antigen summary

### Interface and texts

* Complete introduction text