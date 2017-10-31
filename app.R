library(shiny)
library(shinyBS)
library(plotly)
library(ggplot2)


data = read.table("naseq_cohort_sample_annotation_fschischlik_2017_10.tsv", sep="\t",
	stringsAsFactors = FALSE, header=T, comment.char="") # Import data file

data = data[data$flag == 1, -1] # Remove null-flagged patients

data$patient = factor(unlist(data.frame(strsplit(data$unique.sample.id, "#"))[1,])) # Extract unique patient ID from unique sample ID


# Define client UI
shinyUi <- navbarPage(title = "MPN cohort data vizualization",
	# Starting tab
	tabPanel(title = "What's this?"),

	# Cohort presentation tab
	tabPanel(title = "Explore the cohort"),
)

# Define R server
shinyServer <- function(input, output) {}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)

