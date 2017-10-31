library(shiny)
library(shinyBS)
library(plotly)
library(ggplot2)


# Define constants
color.palette = c()
color.palette$main = "#FF6600"
color.palette$second = "#FFFFFF"
color.palette$bg = "#FFFFFF"

dataCohort = read.table("naseq_cohort_sample_annotation_fschischlik_2017_10.tsv", sep="\t",
	stringsAsFactors = FALSE, header=T, comment.char="") # Import data file

dataCohort = dataCohort[dataCohort$flag == 1, -1] # Remove null-flagged patients

# Convert some columns to dates
dataCohort$date.of.birth = as.Date(dataCohort$date.of.birth,format='%d.%m.%Y')
dataCohort$bleeding.date.dmy = as.Date(dataCohort$bleeding.date.dmy,format='%d.%m.%Y')
dataCohort$extraction.date = as.Date(dataCohort$extraction.date,format='%d.%m.%Y')
dataCohort$date.of.diagnosis = as.Date(dataCohort$date.of.diagnosis,format='%d.%m.%Y')
dataCohort$date.diagnosis.blast.phase = as.Date(dataCohort$date.diagnosis.blast.phase,format='%d.%m.%Y')
dataCohort$date.of.death = as.Date(dataCohort$date.of.death,format='%d.%m.%Y')

# Convert some columns to factors
dataCohort$batch = as.factor(dataCohort$batch)
dataCohort$diagnosis = as.factor(dataCohort$diagnosis)
dataCohort$genotype = as.factor(dataCohort$genotype)
dataCohort$sex = as.factor(dataCohort$sex)

# Create new columns by data transformation
# Extract unique patient ID from unique sample ID
dataCohort$patient = factor(unlist(data.frame(strsplit(dataCohort$unique.sample.id, "#"))[1,])) 
dataCohort$age.days = dataCohort$bleeding.date.dmy - dataCohort$date.of.birth
dataCohort$age.years = as.numeric(round(dataCohort$age.days / 365.25))

# Summarize and sort data by type
dataCohortTypes = data.frame(dataColNames = names(dataCohort), description = c("Unique Sample ID", "Sample ID", "Follow-up", "Batch",
	"Diagnosis", "Genotype", "JAK2 mutant", "CALR mutant", "MPL mutant", "JAK2 burden", "CALR burden", "Date of birth", "Sex",
	"RIN", "Bleeding date", "Extraction date", "Extraction method", "Date of diagnosis", "Date of diagnosis of blast phase",
	"Date of death", "Hemoglobine level at diagnosis", "White blood cell level at diagnosis", "Platelet level at diagnosis",
	"Thrombotic events after diagnosis", "Date of first thrombotic event after diagnosis", "Library type", "Read length",
	"Percentage of mapping reads", "Mapped reads", "Set number", "PLEX number", "Index LSR1", "Flowcell number", "Patient ID",
	"Age when diagnosed (in days)", "Age when diagnosed (in years)"), line = seq(length(dataCohort)), stringsAsFactors = FALSE)

quantitativeVar = c(10, 11, 21, 22, 23, 35, 36)
qualitativeVar = c(5,6,13, 24)

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
			sidebarPanel(
				bsCollapse(
					bsCollapsePanel("Variables",
						selectInput(inputId = "plotVariable",
			  				label = "Data to vizualize",
			  				choices = dataCohortTypes[quantitativeVar,2],
			  				selected = 1,
			  				multiple = FALSE
		  				),
		  				selectInput(inputId = "plotFactor",
			  				label = "Factor by",
			  				choices = c("Nothing", dataCohortTypes[qualitativeVar,2]),
			  				selected = 1,
			  				multiple = FALSE
		  				), 
		  				style = "primary"
					),
					bsCollapsePanel("Graphical parameters",
						sliderInput(inputId = "nbBins",
							label = "Number of bins:",
							min = 1,
							max = 35, 
							value = 10
						), 
						style = "info"
					),
					multiple = TRUE,
					open = "Variables"
				)
				
			),
			# Main displaying panel
			mainPanel(
				tabsetPanel(
					tabPanel(title = "Plot",
						plotlyOutput("histCohort"),  
						verbatimTextOutput("infoClick") # Display click coordinates
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
shinyServer <- function(input, output) {

	# From the variable chosen in the drop-down listbox, return info on corresponding field in the cohort data
	plotVariableInfo <- reactive({
		return(dataCohortTypes[dataCohortTypes[,2] == input$plotVariable,])
	})

	# From the co-factor chosen in the drop-down listbox, return info on corresponding field in the cohort data
	plotCofactorInfo <- reactive({
		return(dataCohortTypes[dataCohortTypes[,2] == input$plotFactor,])
	})

	# Histogram of chosen cohort descriptive variables
	output$histCohort <- renderPlotly({
		dt = plotVariableInfo()
		if(input$plotFactor == "Nothing"){
			gp1 = ggplot(dataCohort, aes_string(dt[[1]])) + geom_histogram(color = color.palette$second, 
				fill = color.palette$main, bins = input$nbBins) + 
				xlab(dt[[2]]) + ylab("Counts")
		}
		else{
			dt2 = plotCofactorInfo()
			print(dt2)
			gp1 = ggplot(dataCohort, aes_string(dt[[1]], fill = dt2[[1]])) + 
				geom_histogram(color = color.palette$second, bins = input$nbBins) + 
				xlab(dt[[2]]) + ylab("Counts")
		}
		margpy1 <- list(l=45, r=5, b=40, t=5) # margins on each side
		gpy1 = ggplotly(gp1) %>% layout(margin=margpy1)
		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))

	})

	# Return text with coordinates of the object clicked
	output$infoClick <- renderText({
		click_event = event_data("plotly_click")
		if(is.null(click_event)){
 			"Click on a bar to get its corresponding value"
		}
		else{
			paste0("x=", click_event$x, "\ny=", click_event$y)
		}
	})
}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)

