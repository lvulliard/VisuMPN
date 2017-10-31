library(shiny)
library(shinyBS)
library(plotly)
library(ggplot2)


data = read.table("naseq_cohort_sample_annotation_fschischlik_2017_10.tsv", sep="\t",
	stringsAsFactors = FALSE, header=T, comment.char="") # Import data file

data = data[data$flag == 1, -1] # Remove null-flagged patients

# Convert some columns to dates
data$date.of.birth = as.Date(data$date.of.birth,format='%d.%m.%Y')
data$bleeding.date.dmy = as.Date(data$bleeding.date.dmy,format='%d.%m.%Y')
data$extraction.date = as.Date(data$extraction.date,format='%d.%m.%Y')
data$date.of.diagnosis = as.Date(data$date.of.diagnosis,format='%d.%m.%Y')
data$date.diagnosis.blast.phase = as.Date(data$date.diagnosis.blast.phase,format='%d.%m.%Y')
data$date.of.death = as.Date(data$date.of.death,format='%d.%m.%Y')

# Convert some columns to factors
data$batch = as.factor(data$batch)
data$diagnosis = as.factor(data$diagnosis)
data$genotype = as.factor(data$genotype)
data$sex = as.factor(data$sex)

# Create new columns by data transformation
data$patient = factor(unlist(data.frame(strsplit(data$unique.sample.id, "#"))[1,])) # Extract unique patient ID from unique sample ID
data$age.days = data$bleeding.date.dmy - data$date.of.birth
data$age.years = as.numeric(round(data$age.days / 365.25))

dataTypes = data.frame(dataColNames = names(data), description = c("Unique Sample ID", "Sample ID", "Follow-up", "Batch",
	"Diagnosis", "Genotype", "JAK2 mutant", "CALR mutant", "MPL mutant", "JAK2 burden", "CALR burden", "Date of birth", "Sex",
	"RIN", "Bleeding date", "Extraction date", "Extraction method", "Date of diagnosis", "Date of diagnosis of blast phase",
	"Date of death", "Hemoglobine level at diagnosis", "White blood cell level at diagnosis", "Platelet level at diagnosis",
	"Thrombotic events after diagnosis", "Date of first thrombotic event after diagnosis", "Library type", "Read length",
	"Percentage of mapping reads", "Mapped reads", "Set number", "PLEX number", "Index LSR1", "Flowcell number", "Patient ID",
	"Age when diagnosed (in days)", "Age when diagnosed (in years)"),
	dataColNum = seq(length(data)))

# Define client UI
shinyUi <- navbarPage(title = "MPN cohort data vizualization",
	# Starting tab
	tabPanel(title = "What's this?",
		p(paste("This web application allows the user to explore interactively data from a cohort of 113 myeloproliferative neoplasm",
			"(MPN) blood donors and 15 healthy blood controls."))
	),

	# Cohort presentation tabs
	navbarMenu(title = "Explore the cohort",
		# Violin plots tab
		tabPanel(title = "Clinical data",
			# Sidebar with plot options
			sidebarPanel(p("Lorem ipsum.")
			),
			# Main displaying panel
			mainPanel(
				tabsetPanel(
					tabPanel(title = "Plot",
						p("Lorem ipsum.")
					),
					tabPanel(title = "Data",
						p("Lorem ipsum.")
					)
				)
			)
		)
	)
)

# Define R server
shinyServer <- function(input, output) {}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)

