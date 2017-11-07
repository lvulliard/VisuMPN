library(shiny)
library(shinyBS)
library(plotly)
library(ggplot2)
library(RColorBrewer)

# Define functions

# Allow addition of uneval types (to add ggplot aes and aes_string outputs)
`+.uneval` <- function(a,b) {
    `class<-`(modifyList(a,b), "uneval")
}


# Define constants
color.palette = c()
color.palette$main = "#FF6600"
color.palette$second = "#FFFFFF"
color.palette$bg = "#FFFFFF"


# Load cohort data

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
dataCohort$thrombotic.events.after.diagn[dataCohort$thrombotic.events.after.diagn == ""] <- "unknown"
dataCohort$thrombotic.events.after.diagn = as.factor(dataCohort$thrombotic.events.after.diagn)

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

# Load variant data
# Import data file
dataVariants = read.table("rnaseq_varcall_patient_mutation_fixids_03_2017.tsv", 
	sep="\t", header=T, comment.char="", stringsAsFactors = FALSE) 
# Remove null-flagged patients
dataCohort$sample.id.variant.file.format = paste(dataCohort$sample.id, dataCohort$batch, sep = "_")
dataVariants = dataVariants[as.character(dataVariants$UNIQ_SAMPLE_ID) %in% dataCohort$sample.id.variant.file.format, ] 

# Join the diagnosis and JAK2/CALR mutation status to the variant table
dataVariants = merge(x = dataVariants, y = dataCohort[,names(dataCohort) %in% c("sample.id.variant.file.format", "diagnosis", "jak2", "calr")],
	by.x = "UNIQ_SAMPLE_ID", by.y = "sample.id.variant.file.format", all.x = T)

# Number of normal variants in the data before filtering
nbIniNormVar = sum(dataVariants$diagnosis == "Control")

# Create unique ID for all variants
dataVariants$VAR_ID <- paste0(dataVariants$CHROM, ":", 
	dataVariants$POS, ":", 
	dataVariants$REF, ":", 
	dataVariants$ALT)

# Convert Minor Allele Frequencies
dataVariants[dataVariants$X1000g2015feb_eur==".",]$X1000g2015feb_eur <- NA
dataVariants[dataVariants$X1000g2015feb_amr==".",]$X1000g2015feb_amr <- NA
dataVariants[dataVariants$X1000g2015feb_eas==".",]$X1000g2015feb_eas <- NA
dataVariants[dataVariants$X1000g2015feb_sas==".",]$X1000g2015feb_sas <- NA
dataVariants[dataVariants$X1000g2015feb_afr==".",]$X1000g2015feb_afr <- NA

dataVariants$X1000g2015feb_eur <- as.numeric(dataVariants$X1000g2015feb_eur)
dataVariants$X1000g2015feb_amr <- as.numeric(dataVariants$X1000g2015feb_amr)
dataVariants$X1000g2015feb_eas <- as.numeric(dataVariants$X1000g2015feb_eas)
dataVariants$X1000g2015feb_sas <- as.numeric(dataVariants$X1000g2015feb_sas)
dataVariants$X1000g2015feb_afr <- as.numeric(dataVariants$X1000g2015feb_afr)

# Define client UI
shinyUi <- navbarPage(title = "MPN cohort data visualization",
	# Starting tab
	tabPanel(title = "What's this?",
		p(paste("This web application allows the user to explore interactively data from a cohort of 113 myeloproliferative neoplasm",
			"(MPN) blood donors and 15 healthy blood controls."))
	),

	# Cohort presentation tabs
	navbarMenu(title = "Explore the cohort",
		# Histograms tab
		tabPanel(title = "Clinical data - Histograms",
			# Sidebar with plot options
			sidebarPanel(
				bsCollapse(
					bsCollapsePanel("Variables",
						selectInput(inputId = "plotVariable",
							label = "Data to visualize:",
							choices = dataCohortTypes[quantitativeVar,2],
							selected = 1,
							multiple = FALSE
						),
						selectInput(inputId = "plotFactor",
							label = "Factor by:",
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
						verbatimTextOutput("infoClick"), # Display click coordinates
						id = "histCohortPlotTab"
					),
					tabPanel(title = "Data",
						dataTableOutput("histCohortTable"),
						id = "histCohortDataTab"
					),
					id = "histCohortTabs"
				)
			)
		),

		# Scatter plots tab
		tabPanel(title = "Clinical data - Scatter plots",
			# Sidebar with plot options
			sidebarPanel(
				bsCollapse(
					bsCollapsePanel("Variables",
						selectInput(inputId = "plotVariableX",
							label = "Data on x-axis:",
							choices = dataCohortTypes[quantitativeVar,2],
							selected = "Hemoglobine level at diagnosis",
							multiple = FALSE
						),
						selectInput(inputId = "plotVariableY",
							label = "Data on y-axis:",
							choices = dataCohortTypes[quantitativeVar,2],
							selected = "White blood cell level at diagnosis",
							multiple = FALSE
						),
						selectInput(inputId = "plotVariableCol",
							label = "Color by:",
							choices = c("Nothing", dataCohortTypes[qualitativeVar,2]),
							selected = 1,
							multiple = FALSE
						),
						style = "primary"
					),
					multiple = TRUE,
					open = "Variables"
				)
				
			),
			# Main displaying panel
			mainPanel(
				tabsetPanel(
					tabPanel(title = "Plot",
						plotlyOutput("scatCohort"),
						id = "ScatCohortPlotTab"
					),
					tabPanel(title = "Data",
						dataTableOutput("scatCohortTable"),
						id = "ScatCohortDataTab"
					),
					id = "ScatCohortTabs"
				)
			)
		),

		# Violin plots tab
		tabPanel(title = "Clinical data - Violin plots",
			# Sidebar with plot options
			sidebarPanel(
				bsCollapse(
					bsCollapsePanel("Variables",
						selectInput(inputId = "plotVariableViolin",
							label = "Data to visualize:",
							choices = dataCohortTypes[quantitativeVar,2],
							selected = 1,
							multiple = FALSE
						),
						selectInput(inputId = "plotFactorViolin",
							label = "Factor by:",
							choices = dataCohortTypes[qualitativeVar,2],
							selected = 1,
							multiple = FALSE
						), 
						style = "primary"
					),
					bsCollapsePanel("Graphical parameters",
						sliderInput(inputId = "jitterViolin",
							label = "Jittering parameter:",
							min = 0,
							max = 0.4, 
							value = 0.05
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
						plotlyOutput("violCohort"),
						id = "ViolCohortPlotTab"
					),
					tabPanel(title = "Data",
						dataTableOutput("violCohortTable"),
						id = "ViolCohortDataTab"
					),
					id = "ViolCohortTabs"
				)
			)
		),
		tabPanel(title = "Clinical data - Matrices and pie charts",
			# Sidebar with plot options
			sidebarPanel(
				bsCollapse(
					bsCollapsePanel("Variables",
						checkboxGroupInput(inputId = "plotVariablesPie",
							label = "Data to visualize:",
                     		choices = dataCohortTypes[qualitativeVar,2],
                     		selected = dataCohortTypes[qualitativeVar,2]
						), 
						style = "primary"
					),
					bsCollapsePanel("Graphical parameters",
						radioButtons("pieOrMatrix", label = "Visualization type:",
							choices = c("Matrix", "Pie chart"),
							selected = "Matrix", inline = TRUE),
						conditionalPanel(
							condition = "input.pieOrMatrix == \"Pie chart\"",
								checkboxInput(inputId = "pieLogScale", label = "Radial log-transformation", value = TRUE)),
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
						conditionalPanel(
							condition = "input.pieOrMatrix == \"Matrix\"",
								plotlyOutput("matrixCohort")),
						conditionalPanel(
							condition = "input.pieOrMatrix == \"Pie chart\"",
								plotOutput("pieCohort")),
						id = "PieCohortPlotTab"
					),
					tabPanel(title = "Data",
						dataTableOutput("pieCohortTable"),
						id = "PieCohortDataTab"
					),
					id = "PieCohortTabs"
				)
			)
		)
	),

	# Variants presentation tabs
	navbarMenu(title = "Explore the variants",
		tabPanel(title = "Variants data - Filtering",
			mainPanel(
				HTML(paste0("<p>From the RNA sequencing of the samples, variant calling is performed using GATK.<br/>",
					"These variants should be filtered to integrate the biological information available, in order ",
					"to avoid false positives. Moreover in the absence of matched normal samples, the litterature ",
					"should be investigated to differentiate as precisely as possible the somatic variants from ",
					"the germline variants.<br/>The healthy control samples are used to test the stringency of the ",
					"filters, and the variants predicted in normal samples are removed from further analyses.<br/>",
					"The following filters can be parametrized:</p>")),
				bsCollapse(
					bsCollapsePanel("Minor Allele Frequency",
						bsModal("modalMAFFilter", "Minor Allele Frequency", "modMAFFilterLink", size = "large",
    						HTML(paste0("You can filter by Minor Allele Frequency (MAF) estimated either in the whole 1000 genomes ",
    							"project or by super population.<br/>Variants often found in the healthy cohort of the 1000 genomes ",
    							"project are likely to be germline variants and not somatic variants.<br/>NB.: variants with",
    							"unknown MAF values are unaffected by this filter.")) 
    					), # Information pop-up
    					div(bsButton("modMAFFilterLink", label = " More info", icon = icon("info")),
    						style="float:right"),
						sliderInput(inputId = "varMAFFilterALL",
							label = "Maximum MAF in 1K Genomes Project in the whole cohort:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterAFR",
							label = "Maximum MAF in 1K Genomes Project for African super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterAMR",
							label = "Maximum MAF in 1K Genomes Project for American super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterEAS",
							label = "Maximum MAF in 1K Genomes Project for East Asian super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterEUR",
							label = "Maximum MAF in 1K Genomes Project for European super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterSAS",
							label = "Maximum MAF in 1K Genomes Project for South Asian super population:",
							min = 0,
							max = 1, 
							value = 0.01
						), 
						style = "info"
					)
				),
				id = "varFilters"
			),
			sidebarPanel(
				uiOutput("varFilterStringencyUI"),
				id = "varFiltersSidebar"
			)
		),
		tabPanel(title = "Variants data - Co-occurence of mutations"),
		tabPanel(title = "Variants data - Occurence of mutations per disease"),
		tabPanel(title = "Variants data - Data",
			dataTableOutput("varTable"),
			id = "varDataTab")
	)
)

# Define R server
shinyServer <- function(input, output) {

	# For a variable chosen in an input of the UI, return info on corresponding field in the cohort data
	variableInfo = function(var){
		return(reactive({dataCohortTypes[dataCohortTypes[,2] == var,]}))
	}

	# Histogram of chosen cohort descriptive variables
	output$histCohort <- renderPlotly({
		dt = variableInfo(input$plotVariable)()
		if(input$plotFactor == "Nothing"){
			gp1 = ggplot(dataCohort, aes_string(dt[[1]])) + geom_histogram(color = color.palette$second, 
				fill = color.palette$main, bins = input$nbBins) + 
				xlab(dt[[2]]) + ylab("Counts")
		}
		else{
			dt2 = variableInfo(input$plotFactor)()
			gp1 = ggplot(dataCohort, aes_string(dt[[1]], fill = dt2[[1]])) + 
				geom_histogram(color = color.palette$second, bins = input$nbBins) + 
				xlab(dt[[2]]) + ylab("Counts")
		}
		margpy1 <- list(l=45, r=5, b=40, t=5) # margins on each side
		gpy1 = ggplotly(gp1 + theme_light()) %>% layout(margin=margpy1)
		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))
	})

	# Scatter plot of chosen cohort descriptive variables
	output$scatCohort <- renderPlotly({
		dtX = variableInfo(input$plotVariableX)()
		dtY = variableInfo(input$plotVariableY)()
		if(input$plotVariableCol == "Nothing"){
			gp1 = ggplot(dataCohort, aes(text = unique.sample.id) + aes_string(dtX[[1]], dtY[[1]])) + geom_point(color = color.palette$main) + 
				xlab(dtX[[2]]) + ylab(dtY[[2]])
		}
		else{
			dtC = variableInfo(input$plotVariableCol)()
			gp1 = ggplot(dataCohort, aes(text = unique.sample.id) + aes_string(dtX[[1]], dtY[[1]])) + geom_point(aes_string(color = dtC[[1]])) + 
				xlab(dtX[[2]]) + ylab(dtY[[2]])
		}
		gpy1 = ggplotly(gp1 + theme_light())
		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))

	})

	# Violin plots of chosen cohort descriptive variables
	output$violCohort <- renderPlotly({
		dt = variableInfo(input$plotVariableViolin)()
		dt2 = variableInfo(input$plotFactorViolin)()
		gp1 = ggplot(dataCohort, aes_string(x = dt2[[1]], y = dt[[1]])) + geom_violin(aes_string(fill = dt2[[1]])) +
			geom_jitter(height = 0, width = input$jitterViolin) +	xlab(dt2[[2]]) + ylab(dt[[2]]) + theme(legend.position="none")
		margpy1 <- list(l=60, r=5, b=40, t=5) # margins on each side
		gpy1 = ggplotly(gp1 + theme_light()) %>% layout(margin=margpy1)

		# Manually modify hoverinfo
		nbFactors = length(levels(dataCohort[,dt2[[3]]]))
		for(i in 1:nbFactors){
			# Remove duplicated annotation of qualitative variable on density lines
			gpy1$x$data[[i]]$text = gsub("^[^>]*>", "", gpy1$x$data[[i]]$text)
		}
		# Add sample ID to each point
		gpy1$x$data[[nbFactors+1]]$text = paste0(gpy1$x$data[[nbFactors+1]]$text, paste0("<br />", dataCohort$unique.sample.id))


		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))
	})

	# Pie chart or matrix of chosen cohort descriptive variables
	output$matrixCohort <- renderPlotly({
		# Check if variables are selected
		nbVarSelected = length(input$plotVariablesPie) 
		if(nbVarSelected == 0){
			return(ggplot())
		}

		# Sort the data by columns in input vector and count total number of factors in plot
		sortedDataCohort = dataCohort
		nbFactors = 0
		for(i in rev(input$plotVariablesPie)){
			dt = variableInfo(i)()
			sortedDataCohort = sortedDataCohort[order(sortedDataCohort[,dt[[3]]]),]
			nbFactors = nbFactors + length(levels(sortedDataCohort[,dt[[3]]]))
		}

		gp1 = ggplot(sortedDataCohort, aes(text = unique.sample.id))
		for(i in 1:nbVarSelected){
			dt = variableInfo(input$plotVariablesPie[i])()
			gp1 = gp1 + geom_rect(aes_string(fill=dt[[1]]) +
				aes_(ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=i+1, xmin=i))
		}			
		gp1 = gp1 + theme(aspect.ratio=1) + 
			scale_x_continuous(breaks = 1:nbVarSelected + 0.5, labels = input$plotVariablesPie) +
			scale_fill_manual(guide = guide_legend(title = NULL), values = colorRampPalette(brewer.pal(12, "Set3"))(nbFactors))
  		
  		gpy1 = ggplotly(gp1 + theme_light()) 

		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))
	})

	output$pieCohort <- renderPlot({
		# Check if variables are selected
		nbVarSelected = length(input$plotVariablesPie) 
		if(nbVarSelected == 0){
			return(ggplot())
		}

		# Sort the data by columns in input vector and count total number of factors in plot
		sortedDataCohort = dataCohort
		nbFactors = 0
		for(i in rev(input$plotVariablesPie)){
			dt = variableInfo(i)()
			sortedDataCohort = sortedDataCohort[order(sortedDataCohort[,dt[[3]]]),]
			nbFactors = nbFactors + length(levels(sortedDataCohort[,dt[[3]]]))
		}

		gp1 = ggplot(sortedDataCohort, aes(text = unique.sample.id))
		for(i in 1:nbVarSelected){
			dt = variableInfo(input$plotVariablesPie[i])()
			gp1 = gp1 + geom_rect(aes_string(fill=dt[[1]]) +
				aes_(ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax= ifelse(input$pieLogScale, log(i+1), i+1),
				xmin=ifelse(input$pieLogScale, log(i), i)))
		}			
		gp1 = gp1 + xlim(c(ifelse(input$pieLogScale, 0, 1), ifelse(input$pieLogScale, log(nbVarSelected+1), nbVarSelected+1))) + theme(aspect.ratio=1) + 
			scale_fill_manual(guide = guide_legend(ncol = nbVarSelected, title = NULL), values = colorRampPalette(brewer.pal(12, "Set3"))(nbFactors)) +
			scale_x_discrete() + # Remove radial legend
			coord_polar(theta="y")
  		
		gp1 + theme_light()
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

	# Return tables with selected data

	output$histCohortTable <- renderDataTable({
		dt = variableInfo(input$plotVariable)()
		colToExport = c(1,dt[[3]])
		if(input$plotFactor != "Nothing"){
			dt2 = variableInfo(input$plotFactor)()
			colToExport = append(colToExport, dt2[[3]])
		}
		tmpframe = data.frame(dataCohort[,colToExport])
		names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
		tmpframe
	})

	output$scatCohortTable <- renderDataTable({
		dt1 = variableInfo(input$plotVariableX)()
		dt2 = variableInfo(input$plotVariableY)()
		colToExport = c(1,dt1[[3]], dt2[[3]])
		if(input$plotVariableCol != "Nothing"){
			dt3 = variableInfo(input$plotVariableCol)()
			colToExport = append(colToExport, dt3[[3]])
		}
		tmpframe = data.frame(dataCohort[,colToExport])
		names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
		tmpframe
	})

	output$violCohortTable <- renderDataTable({
		dt = variableInfo(input$plotVariableViolin)()
		dt2 = variableInfo(input$plotFactorViolin)()
		colToExport = c(1,dt[[3]],dt2[[3]])
		tmpframe = data.frame(dataCohort[,colToExport])
		names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
		tmpframe
	})

	output$pieCohortTable <- renderDataTable({
		colToExport = c(1,unlist(dataCohortTypes[dataCohortTypes[,2] %in% input$plotVariablesPie,3]))
		tmpframe = data.frame(dataCohort[,colToExport])
		names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
		tmpframe
	})

	output$varTable <- renderDataTable({
		filteredDataVariants()
	})

	# Samples (including normal) filtered for the selected parameters, as variant IDs
	filteredSamples <- reactive({
		filtered = c()
		filtered = append(filtered, filterVariant_MAF(filtered)())
		return(filtered)
	})

	# Patient samples left after removing normal variants, as variant IDs
	filteredPatients <- reactive({
		filtered = filteredSamples()
		# Filter out normal sample variants
		filtered = append(filtered, unique(dataVariants[dataVariants$diagnosis == "Control",]$VAR_ID))
		return(filtered)
	})

	# Table of filtered patient variants
	filteredDataVariants <- reactive({
		return(dataVariants[! dataVariants$VAR_ID %in% filteredPatients(),])
	})

	# Output how well the filters remove healthy controls from the dataset
	output$varFilterStringencyUI <- renderUI({
		nbControlVariants = sum(dataVariants$diagnosis == "Control")
		nbControlUnfilteredVariants = sum(dataVariants[! dataVariants$VAR_ID %in% filteredSamples(),]$diagnosis == "Control")

		if(nbControlUnfilteredVariants == 0){
			stringencyPanelName = "Filters seem to have an appropriate stringency"
			stringencyPanelType = "success"
		} else if(nbControlUnfilteredVariants < nbControlVariants) {
			stringencyPanelName = "Filters are insufficiently stringent"
			stringencyPanelType = "warning"
		} else {
			stringencyPanelName = "Filters fail to remove any normal sample"
			stringencyPanelType = "danger"			
		}

		bsCollapsePanel(stringencyPanelName,
			p(paste0(nbControlUnfilteredVariants, " out of ", nbControlVariants,
				" normal variants are kept by applying these filters.")),
			style = stringencyPanelType
		)
	})

	# Filters on the variant data

	# Filter by MAF threshold
	filterVariant_MAF <- function(filtered){
		subset = dataVariants[! dataVariants$VAR_ID %in% filtered,]
		return(reactive({
			unique(subset[!((subset$X1000g2015feb_all < input$varMAFFilterALL | is.na(subset$X1000g2015feb_all)) &
				(subset$X1000g2015feb_afr < input$varMAFFilterAFR | is.na(subset$X1000g2015feb_afr)) &
				(subset$X1000g2015feb_amr < input$varMAFFilterAMR | is.na(subset$X1000g2015feb_amr)) &
				(subset$X1000g2015feb_eas < input$varMAFFilterEAS | is.na(subset$X1000g2015feb_eas)) &
				(subset$X1000g2015feb_eur < input$varMAFFilterEUR | is.na(subset$X1000g2015feb_eur)) &
				(subset$X1000g2015feb_sas < input$varMAFFilterSAS | is.na(subset$X1000g2015feb_sas))),]$VAR_ID)
		}))
	}

}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)

