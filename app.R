library(shiny)
library(shinyBS)
library(plotly)
library(ggplot2)

data = read.table("naseq_cohort_sample_annotation_fschischlik_2017_10.tsv", sep="\t", stringsAsFactors = FALSE)
head(data)
