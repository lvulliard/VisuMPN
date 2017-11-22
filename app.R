library(shiny)
library(shinyBS)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(heatmaply)

# Define functions

# Allow addition of uneval types (to add ggplot aes and aes_string outputs)
`+.uneval` <- function(a,b) {
	`class<-`(modifyList(a,b), "uneval")
}


# Define constants
color.palette = c()
color.palette$main = "#00ace6"
color.palette$second = "#FFFFFF"
color.palette$bg = "#FFFFFF"
color.palette$function_multi = colorRampPalette(brewer.pal(9, "Pastel1"))
# Bimodal color palette going from blue to white (rate x) then slowly to orange (rate x/10) 
color.palette$function_bimod = colorRampPalette(c("#00bfff","#b3ecff",colorRampPalette(c("#ffffff", "#ffcc80", "#ffad33", "#ff9900"))(21) ))


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
dataVariants$X1000g2015feb_eur <- as.numeric(dataVariants$X1000g2015feb_eur) # "." values are replaced by NAs
dataVariants$X1000g2015feb_amr <- as.numeric(dataVariants$X1000g2015feb_amr)
dataVariants$X1000g2015feb_eas <- as.numeric(dataVariants$X1000g2015feb_eas)
dataVariants$X1000g2015feb_sas <- as.numeric(dataVariants$X1000g2015feb_sas)
dataVariants$X1000g2015feb_afr <- as.numeric(dataVariants$X1000g2015feb_afr)

# Compute ALT counts
dataVariants$ALT.COUNT = as.numeric(gsub("^[^,]*,", "", dataVariants$AD))

# Convert SIFT and CADD scores
dataVariants$CADD_phred <- as.numeric(dataVariants$CADD_phred) # "." values are replaced by NAs
dataVariants$SIFT_score <- as.numeric(dataVariants$SIFT_score)

# Create subclasses of ET and PMF based on CALR and JAK2 mutations
dataVariants$subdiagnosis <- paste(dataVariants$diagnosis, ifelse(dataVariants$jak2, "JAK2-mutated", "JAK2-wt"), ifelse(dataVariants$calr, "CALR-mutated", "CALR-wt"), sep = "-")

# Explanatory text displayed in variant heatmaps
modalVarHeatmapText = HTML(paste0("The heatmaps display association between features through odds-ratios observed in the data.",
								"<br/>An odds-ratio of 0 corresponds to mutual exclusion whereas an infitine odds-ratio corresponds ",
								"to a mutual inclusion. For displaying purposes, the infinite values are replaced here by the ",
								"biggest values observed in the data plus one.<br/>Because of the finite sample size, the computed ",
								"association might not be significant, and therefore an exact Fisher test is performed, giving the ",
								"probability of observing such an association by chance if the underlying distribution does not ",
								"favor any inclusion or exclusion (p-value).<br/>To take into account the number of tests performed ",
								"a Benjamini-Hochberg FDR correction is performed. You can then choose to display the odds-ratios ",
								"only for a specified false discovery rate."))

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
						div(downloadButton('histCohortDL', 'Download'),style="float:right"),
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
						div(downloadButton('scatCohortDL', 'Download'),style="float:right"),
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
						div(downloadButton('violCohortDL', 'Download'),style="float:right"),
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
					tabPanel(title = "Filtered data",
						div(downloadButton('pieCohortDL', 'Download'),style="float:right"),
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
								"project are likely to be germline variants and not somatic variants.<br/>NB.: variants with ",
								"unknown MAF values are unaffected by this filter.")) 
						), # Information pop-up
						div(bsButton("modMAFFilterLink", label = " More info", icon = icon("info")),
							style="float:right"),
						sliderInput(inputId = "varMAFFilterALL",
							label = "Maximum MAF in 1000 Genomes Project in the whole cohort:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterAFR",
							label = "Maximum MAF in 1000 Genomes Project for African super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterAMR",
							label = "Maximum MAF in 1000 Genomes Project for American super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterEAS",
							label = "Maximum MAF in 1000 Genomes Project for East Asian super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterEUR",
							label = "Maximum MAF in 1000 Genomes Project for European super population:",
							min = 0,
							max = 1, 
							value = 0.01
						),
						sliderInput(inputId = "varMAFFilterSAS",
							label = "Maximum MAF in 1000 Genomes Project for South Asian super population:",
							min = 0,
							max = 1, 
							value = 0.01
						), 
						style = "info"
					),
					bsCollapsePanel("Annotation and evidence",
						bsModal("modalAnnFilter", "Annotation and evidence", "modAnnFilterLink", size = "large",
							HTML(paste0("You can filter by counts of reads supporting the presence of the variant (ALT counts).<br/>",
								"The more mutant sequencing reads mapping at the variant location, the stronger the evidence ",
								"of the variant.<br/>You may choose to keep variants with low ALT counts if they are annotated ",
								"in the COSMIC database, which means that the variant has already been identified as a cancer ",
								"germline variant in a previous study.<br/>Moreover, even for high ALT counts, you may want to ",
								"filter out variants having a rsID from dbSNP and no COSMIC ID, meaning that the variants have ",
								"already been previously observed but not yet linked to cancer and are therefore likely to be ",
								"germline variants.")) 
						), # Information pop-up
						div(bsButton("modAnnFilterLink", label = " More info", icon = icon("info")),
							style="float:right"),
						sliderInput(inputId = "varAnnFilterALT",
							label = "Minimum ALT count:",
							min = 1,
							max = 50, 
							value = 3
						),
						radioButtons("varAnnFilterDB", label = "Annotation filtering:",
							choices = c("None", "Keep all variants with COSMIC annotations", 
								"Remove variants with dbSNP annotations and keep all variants with COSMIC annotations"),
							selected = "Keep all variants with COSMIC annotations", inline = TRUE),
						style = "info"
					),
					bsCollapsePanel("Non-canonical transcripts",
						bsModal("modalCanFilter", "Non-canonical transcripts", "modCanFilterLink", size = "large",
							HTML(paste0("A small amount of the variants are predicted to have an effect only on a non-canonical ",
								"transcript, as annotated on the ENSEMBL database.<br/>These variants could be false positives ",
								"and you may choose to filter them out."))
						), # Information pop-up
						div(bsButton("modCanFilterLink", label = " More info", icon = icon("info")),
							style="float:right"),
						checkboxInput(inputId = "canonicalTranscriptFilter", label = "Filter out non-canonical transcripts",
							value = TRUE),
						style = "info"
					),
					bsCollapsePanel("Additional filters",
						bsModal("modalAddFilter", "Non-canonical transcripts", "modAddFilterLink", size = "large",
							HTML(paste0("Additionally you may want to apply extra filters to keep variants of interest only.<br/>",
								"SIFT score is a prediction tool measuring the impact of non-synonymous point mutations on ",
								"protein function. The output score corresponds to the probability of the variant to be tolerated ",
								"by the organism. A score belows 0.05 therefore means that the variant is deleterious with ",
								"probability 95%.<br/>",
								"CADD PHRED-like scaled C-score predicts whether a mutation will be deleterious or not, by ",
								"aggregating the output of different prediction algorithms such as SIFT. Higher scores mean ",
								"that the variants are more likely to be pathogenic: a score of 10 indicates a predicted effect ",
								"similar to the top 10% deleterious substitions observable on the human genome whereas a score of ",
								"20 corresponds to the top 1% deleterious mutations.<br/>",
								"Finally the variants can be filtered by mutant allele frequency, to ensure the variant is really ",
								"present on one or both chromosomes."))
						), # Information pop-up
						div(bsButton("modAddFilterLink", label = " More info", icon = icon("info")),
							style="float:right"),
						sliderInput(inputId = "varSIFTFilter",
							label = "Maximum SIFT score:",
							min = 0,
							max = 1, 
							value = 1
						),
						sliderInput(inputId = "varCADDFilter",
							label = "Minimum CADD PHRED-like scaled C-score:",
							min = 0,
							max = 40, 
							value = 0
						),
						sliderInput(inputId = "varFREQFilter",
							label = "Minimum variant allele frequency:",
							min = 0,
							max = 1, 
							value = 0
						)
					)
				),
				id = "varFilters"	
			),
			sidebarPanel(
				uiOutput("varFilterStringencyUI"),
				id = "varFiltersSidebar"
			)
		),
		tabPanel(title = "Variants data - Binary matrix",
			plotlyOutput("varBinMat", height = "600px")
		),
		tabPanel(title = "Variants data - Co-occurrence of mutations",
			tabsetPanel(
				tabPanel(title = "Plot",
					mainPanel(plotlyOutput("varCoOc", height = "550px", width =  "680px")),
					sidebarPanel(
						div(bsButton("modVarHeatmapsLink1", label = " More info", icon = icon("info")),
							style="float:right"),
						sliderInput(inputId = "alphaCoOc",
							label = "Alpha risk (with Benjamini-Hochberg FDR correction):",
							min = 0,
							max = 1, 
							value = 0.04
						),
						sliderInput(inputId = "nbRepCoOc",
							label = "Minimal amount of patients with variants per gene:",
							min=1,
							max=18,
							value = 4
						)
					)
				),
				tabPanel(title = "Data - Odds-ratios",
					div(downloadButton('varCoOcORDL', 'Download'),style="float:right"),
					dataTableOutput("varCoOcORTable")
				),
				tabPanel(title = "Data - Raw p-values",
					div(downloadButton('varCoOcPvalDL', 'Download'),style="float:right"),
					dataTableOutput("varCoOcPvalTable")
				),
				id = "varCoOcTabs"
			)
		),
		tabPanel(title = "Variants data - Occurrence of mutations per disease",
			tabsetPanel(
				tabPanel(title = "Plot",
					mainPanel(plotlyOutput("varDisOc", height = "550px", width =  "680px")),
					sidebarPanel(
						div(bsButton("modVarHeatmapsLink2", label = " More info", icon = icon("info")),
							style="float:right"),
						sliderInput(inputId = "alphaDisOc",
							label = "Alpha risk (with Benjamini-Hochberg FDR correction):",
							min = 0,
							max = 1, 
							value = 0.4
						),
						sliderInput(inputId = "nbRepDisOc",
							label = "Minimal amount of mutations per gene:",
							min=1,
							max=18,
							value = 3
						)
					)
				),
				tabPanel(title = "Data - Odds-ratios",
					div(downloadButton('varDisOcORDL', 'Download'),style="float:right"),
					dataTableOutput("varDisOcORTable")
				),
				tabPanel(title = "Data - Corrected p-values",
					div(downloadButton('varDisOcPvalDL', 'Download'),style="float:right"),
					dataTableOutput("varDisOcPvalTable")
				),
				id = "varDisOcTabs"
			)	
		),
		tabPanel(title = "Variants data - Occurrence of mutations in disease subtypes",
			tabsetPanel(
				tabPanel(title = "Plot",
					mainPanel(plotlyOutput("varSubDisOc", height = "550px", width =  "800px")),
					sidebarPanel(
						div(bsButton("modVarHeatmapsLink3", label = " More info", icon = icon("info")),
							style="float:right"),
						sliderInput(inputId = "alphaSubDisOc",
							label = "Alpha risk (with Benjamini-Hochberg FDR correction):",
							min = 0,
							max = 1, 
							value = 0.98
						),
						sliderInput(inputId = "nbRepSubDisOc",
							label = "Minimal amount of mutations per gene:",
							min=1,
							max=11,
							value = 7
						),
						radioButtons("subTypeAll", label = "Subtypes:",
							choices = c("PMF & ET", "All"),
							selected = "PMF & ET", inline = TRUE)
					)
				),
				tabPanel(title = "Data - Odds-ratios",
					div(downloadButton('varSubDisOcORDL', 'Download'),style="float:right"),
					dataTableOutput("varSubDisOcORTable")
				),
				tabPanel(title = "Data - Corrected p-values",
					div(downloadButton('varSubDisOcPvalDL', 'Download'),style="float:right"),
					dataTableOutput("varSubDisOcPvalTable")
				),
				id = "varSubDisOcTabs"
			)	
		),
		tabPanel(title = "Variants data - Data",
			div(downloadButton('varDL', 'Download'),style="float:right"),
			dataTableOutput("varTable"),
			id = "varDataTab"),
		bsModal("modalVarHeatmaps", "Variants heatmaps", "modVarHeatmapsLink1", size = "large",
							modalVarHeatmapText),
		bsModal("modalVarHeatmaps", "Variants heatmaps", "modVarHeatmapsLink2", size = "large",
							modalVarHeatmapText),
		bsModal("modalVarHeatmaps", "Variants heatmaps", "modVarHeatmapsLink3", size = "large",
							modalVarHeatmapText) # Information pop-up on fisher tests
	)
)

# Define R server
shinyServer <- function(input, output) {

	# For a variable chosen in an input of the UI, return info on corresponding field in the cohort data
	variableInfo <- function(var){
		return(reactive({dataCohortTypes[dataCohortTypes[,2] == var,]}))
	}

	histCohortObject <- reactive({
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
		gpy1 = style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))
	})

	# Histogram of chosen cohort descriptive variables
	output$histCohort <- renderPlotly({
		histCohortObject()
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
		# Columns with one point or less don't have a density component
		nbFactors = sum(table(dataCohort[!is.na(dataCohort[,dt[[3]]]),dt2[[3]]]) > 1) 
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
		factorGroups = c()
		for(i in rev(input$plotVariablesPie)){
			dt = variableInfo(i)()
			sortedDataCohort = sortedDataCohort[order(sortedDataCohort[,dt[[3]]]),]
			nbFactor = length(levels(sortedDataCohort[,dt[[3]]]))
			nbFactors = nbFactors + nbFactor
			factorGroups = append(factorGroups, rep(dt[[2]], nbFactor))
		}

		gp1 = ggplot(sortedDataCohort, aes(text = unique.sample.id))
		for(i in 1:nbVarSelected){
			dt = variableInfo(input$plotVariablesPie[i])()
			gp1 = gp1 + geom_rect(aes_string(fill=dt[[1]]) +
				aes_(ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=i+1, xmin=i))
		}			
		gp1 = gp1 + theme(aspect.ratio=1) + 
			scale_x_continuous(breaks = 1:nbVarSelected + 0.5, labels = input$plotVariablesPie) +
			scale_fill_manual(guide = guide_legend(title = NULL), values = color.palette$function_multi(nbFactors))
 		
		gpy1 = ggplotly(gp1 + theme_light()) 

		# Group legends
		factorGroups = rev(factorGroups)
		elementsLegend = sapply(gpy1$x$data, function(x) x$showlegend)
		for(i in 1:sum(elementsLegend)){
			gpy1$x$data[[which(elementsLegend)[i]]]$legendgroup = factorGroups[i]
		}

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
			scale_fill_manual(guide = guide_legend(ncol = nbVarSelected, title = NULL), values = color.palette$function_multi(nbFactors)) +
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
			binsCenters = histCohortObject()$x$data[[1]]$x
			complement = ""
			if(length(binsCenters) > 1){
				binsize = binsCenters[2] - binsCenters[1]
				values = variableInfo(input$plotVariable)()[[1]]
				patientInBin = dataCohort[ (dataCohort[,values] <= click_event$x + binsize/2)&
					(dataCohort[,values]  >= click_event$x - binsize/2),]$unique.sample.id
				patientInBin = patientInBin[!is.na(patientInBin)]
				complement = paste(patientInBin, collapse="\n")
				complement = paste0("\nSamples: ", complement)
			}
			paste0("x=", click_event$x, "\ny=", click_event$y, complement)
		}
	})

	# Return tables with selected data and generate corresponding files

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

	output$histCohortDL <- downloadHandler(
		filename = "cohort_histogram.csv",
		content = function(file) {
			dt = variableInfo(input$plotVariable)()
			colToExport = c(1,dt[[3]])
			if(input$plotFactor != "Nothing"){
				dt2 = variableInfo(input$plotFactor)()
				colToExport = append(colToExport, dt2[[3]])
			}
			tmpframe = data.frame(dataCohort[,colToExport])
			names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
			write.csv(tmpframe, file)
		}
	)

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

	output$scatCohortDL <- downloadHandler(
		filename = "cohort_scatter_plot.csv",
		content = function(file) {
			dt1 = variableInfo(input$plotVariableX)()
			dt2 = variableInfo(input$plotVariableY)()
			colToExport = c(1,dt1[[3]], dt2[[3]])
			if(input$plotVariableCol != "Nothing"){
				dt3 = variableInfo(input$plotVariableCol)()
				colToExport = append(colToExport, dt3[[3]])
			}
			tmpframe = data.frame(dataCohort[,colToExport])
			names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
			write.csv(tmpframe, file)
		}
	)

	output$violCohortTable <- renderDataTable({
		dt = variableInfo(input$plotVariableViolin)()
		dt2 = variableInfo(input$plotFactorViolin)()
		colToExport = c(1,dt[[3]],dt2[[3]])
		tmpframe = data.frame(dataCohort[,colToExport])
		names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
		tmpframe
	})

	output$violCohortDL <- downloadHandler(
		filename = "cohort_violin_plot.csv",
		content = function(file) {
			dt = variableInfo(input$plotVariableViolin)()
			dt2 = variableInfo(input$plotFactorViolin)()
			colToExport = c(1,dt[[3]],dt2[[3]])
			tmpframe = data.frame(dataCohort[,colToExport])
			names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
			write.csv(tmpframe, file)
		}
	)

	output$pieCohortTable <- renderDataTable({
		colToExport = c(1,unlist(dataCohortTypes[dataCohortTypes[,2] %in% input$plotVariablesPie,3]))
		tmpframe = data.frame(dataCohort[,colToExport])
		names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
		tmpframe
	})

	output$pieCohortDL <- downloadHandler(
		filename = "cohort_pie_chart.csv",
		content = function(file) {
			colToExport = c(1,unlist(dataCohortTypes[dataCohortTypes[,2] %in% input$plotVariablesPie,3]))
			tmpframe = data.frame(dataCohort[,colToExport])
			names(tmpframe) = unlist(dataCohortTypes[colToExport,2])
			write.csv(tmpframe, file)
		}
	)


	# Variant exploration

	output$varTable <- renderDataTable({
		filteredDataVariants()
	})

	output$varDL <- downloadHandler(
		filename = "variants_filtered.csv",
		content = function(file) {
			write.csv(filteredDataVariants(), file)
		}
	)

	# Samples (including normal) filtered for the selected parameters, as variant IDs
	filteredSamples <- reactive({
		filtered = c()
		filtered = append(filtered, filterVariant_MAF(filtered)())
		filtered = append(filtered, filterVariant_Ann(filtered)())
		filtered = append(filtered, filterVariant_Can(filtered)())
		filtered = append(filtered, filterVariant_SIFT(filtered)())
		filtered = append(filtered, filterVariant_CADD(filtered)())
		filtered = append(filtered, filterVariant_Freq(filtered)())
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

	# Filter by ALT score and annotation
	filterVariant_Ann <- function(filtered){
		subset = dataVariants[! dataVariants$VAR_ID %in% filtered,]
		return(reactive({
			if(input$varAnnFilterDB == "None"){
				subset[subset$ALT.COUNT < input$varAnnFilterALT,]$VAR_ID		
			} else if(input$varAnnFilterDB == "Keep all variants with COSMIC annotations") {
				subset[(subset$ALT.COUNT < input$varAnnFilterALT) & (subset$cosmic70 == "."),]$VAR_ID
			}
			else{
				append(subset[(subset$ALT.COUNT < input$varAnnFilterALT) & (subset$cosmic70 == "."),]$VAR_ID,
					subset[(subset$snp129 != ".") & (subset$cosmic70 == "."),]$VAR_ID)
			}
		}))
	}

	# Filter out non-canonical transcripts
	filterVariant_Can <- function(filtered){
		subset = dataVariants[! dataVariants$VAR_ID %in% filtered,]
		return(reactive({
			if(input$canonicalTranscriptFilter){
				unique(subset[subset$IS_CANONICAL_TRANSCRIPT=="N",]$VAR_ID)
			} 
		}))
	}

	# Filter on SIFT score
	filterVariant_SIFT <- function(filtered){
		subset = dataVariants[! dataVariants$VAR_ID %in% filtered,]
		return(reactive({
			unique(subset[(subset$SIFT_score > input$varSIFTFilter)&!(is.na(subset$SIFT_score)),]$VAR_ID)
		}))
	}

	# Filter on CADD score
	filterVariant_CADD <- function(filtered){
		subset = dataVariants[! dataVariants$VAR_ID %in% filtered,]
		return(reactive({
			unique(subset[(subset$CADD_phred < input$varCADDFilter)&!(is.na(subset$CADD_phred)),]$VAR_ID)
		}))
	}

	# Filter on mutant frequency
	filterVariant_Freq <- function(filtered){
		subset = dataVariants[! dataVariants$VAR_ID %in% filtered,]
		return(reactive({
			unique(subset[(subset$VARIANT_FREQUENCY < input$varFREQFilter)&!(is.na(subset$VARIANT_FREQUENCY)),]$VAR_ID)
		}))
	}

	# Gene mutations co-occurrence
	coOcOR <- reactive({
		dataset = filteredDataVariants()
		variantsPerPatient = table(dataset$GENESYMBOL, dataset$UNIQ_SAMPLE_ID) != 0 # Discard the number of mutation per patient

		n = dim(variantsPerPatient)[2] # Number of patients with mutations

		# Remove genes with insufficient number of variants observed in the dataset
		variantsPerPatient = variantsPerPatient[rowSums(variantsPerPatient) >= input$nbRepCoOc, ]
	
		# Test if at least 2 genes remain
		if(length(variantsPerPatient) <= n) {print("No data to display.");return(list("Nothing to display.","Nothing to display."))}

		nbGenes = dim(variantsPerPatient)[1] # Number of genes that can be mutated several times
		n = dim(variantsPerPatient)[2] # Number of patients with mutations

		fisherTests = sapply(1:nbGenes, function(y) sapply(y:nbGenes, function(x) fisher.test(variantsPerPatient[y,], 
			variantsPerPatient[x,]))) # Test is co-occurrence observed can be expected by chance for random association rates
		
		pvalMat <- ORMat <- matrix(nrow =  nbGenes, ncol = nbGenes)
		rownames(pvalMat) <- rownames(ORMat) <- colnames(pvalMat) <- colnames(ORMat) <- rownames(variantsPerPatient)

		nbTests = (nbGenes-1)*nbGenes/2 # Non-self cooccurrence tests
		all_pval = rep(NA, nbTests) # Store all non-diagonal p-values for BH correction
		for(x in 1:nbGenes){
			for(y in 1:(nbGenes-x+1)){
				fTest = fisherTests[[x]][,y]
				pvalMat[x,x+y-1] <- pvalMat[x+y-1,x] <- fTest$p.value # Matrix of non-adjusted p-values
				ORMat[x,x+y-1] <- ORMat[x+y-1,x] <- fTest$estimate # Matrix of odds-ratios
				
				if(y!=1){all_pval[y-1+(x-1)*(2*nbGenes - x)/2] = fTest$p.value}
			}
		}

		ORMat[ORMat == Inf] <- -1 # Replace infinite OR (mutual association) by 2nd biggest value + 1
		ORMat[ORMat == -1] <- max(ORMat) + 1

		# Benjamini-Hochberg
		all_pval_bh = p.adjust(all_pval, method="BH")

		for(x in 1:nbGenes){
			for(y in 1:(nbGenes-x+1)){
				# If not a diagonal case and corrected p-value too high
				if((y!=1)&&(all_pval_bh[y-1+(x-1)*(2*nbGenes - x)/2] > input$alphaCoOc)){
					# Remove value from matrix
					ORMat[x,x+y-1] <- ORMat[x+y-1,x] <- NA
				}
			}
		}

		# Remove rows without any information
		genesToKeep = colSums(ORMat, na.rm=T) != max(ORMat, na.rm=T)
		
		# Test if at least 2 genes remain
		if(sum(genesToKeep) < 2) {print("No data to display.");return(list("Nothing to display.","Nothing to display."))}

		ORMat = ORMat[genesToKeep, genesToKeep]

		return(list(ORMat, pvalMat))		
	})

	output$varCoOc <- renderPlotly({
		ORMat = coOcOR()[[1]]

		if(class(ORMat) != "matrix"){return()} # No data to display

		heatmaply(ORMat, symm = T, na.rm = T, colors = color.palette$function_bimod, show_grid = T, dendrogram = "none",
			na.value=color.palette$function_bimod(256)[256], margins = c(75,75,NA,0), limits = c(0,20), grid_gap=2, 
			label_names = c("Row", "Column", "OR"), colorbar_len = 1, key.title = "Odds-ratios")
	})

	output$varCoOcORTable <- renderDataTable({
		tab = coOcOR()[[1]]
		return(cbind(gene = rownames(tab), tab))
	})

	output$varCoOcPvalTable <- renderDataTable({
		tab = coOcOR()[[2]]
		return(cbind(gene = rownames(tab), tab))
	})

	output$varCoOcORDL <- downloadHandler(
		filename = "variants_gene_gene_OR.csv",
		content = function(file) {
			write.csv(coOcOR()[[1]], file)
		}
	)	

	output$varCoOcPvalDL <- downloadHandler(
		filename = "variants_gene_gene_pval.csv",
		content = function(file) {
			write.csv(coOcOR()[[2]], file)
		}
	)

	# Common fisher test function for mutation per subtype and per disease
	fisherMutationMatrices <- function(variantsPerDisease, alphaRisk){
		n = dim(variantsPerDisease)
		nbMut = sum(variantsPerDisease)

		ORMat <- pvalMat <- as.data.frame.matrix(variantsPerDisease)

		# Fisher's exact tests on contigency tables for each pair of gene and disease
		for(x in 1:n[1]){
			for(y in 1:n[2]){
				cntgTab = data.frame(c(variantsPerDisease[x,y], sum(variantsPerDisease[,y]) - variantsPerDisease[x,y]),
					c(sum(variantsPerDisease[x,]) - variantsPerDisease[x,y], 
						nbMut + variantsPerDisease[x,y] - sum(variantsPerDisease[,y]) - sum(variantsPerDisease[x,]) ))
				fTest = fisher.test(cntgTab)
				ORMat[x,y] = fTest$estimate
				pvalMat[x,y] = fTest$p.value
			}
		}

		ORMat[ORMat == Inf] <- -1 # Replace infinite OR (mutual association) by 2nd biggest value + 1
		ORMat[ORMat == -1] <- max(ORMat) + 1

		pvalMat = matrix(p.adjust(unlist(pvalMat), method="BH"), ncol=ncol(pvalMat), byrow = F) # Benjamini-Hochberg FDR
		colnames(pvalMat) = colnames(ORMat)
		rownames(pvalMat) = rownames(ORMat)
		is.na(ORMat) <- pvalMat > alphaRisk # Convert non-significant odds-ratios to NA

		# Remove rows without any information
		ORMat[!is.na(ORMat)] <- ORMat[!is.na(ORMat)]+0.1 # Ensure that the non-NA values are positive
		genesToKeep = colSums(ORMat, na.rm=T) > 0 # Check the sums by columns and rows to detect empty features
		diagnosesToKeep = rowSums(ORMat, na.rm=T) > 0
		ORMat[!is.na(ORMat)] <- ORMat[!is.na(ORMat)]-0.1 # Reverse to original OR values

		return(list(genesToKeep, diagnosesToKeep, ORMat, pvalMat))
	}



	# Gene mutations per subtype

	subDisOcOR <- reactive({
		dataset = filteredDataVariants()
		variantsPerDisease = table(dataset$subdiagnosis, dataset$GENESYMBOL)
		variantsPerDisease = variantsPerDisease[rowSums(variantsPerDisease) > 0, colSums(variantsPerDisease) > 0]
		
		# Remove genes with insufficient number of variants observed in the dataset
		variantsPerDisease = variantsPerDisease[,colSums(variantsPerDisease) >= input$nbRepSubDisOc]
	
		# Remove JAK2 and CALR
		variantsPerDisease = variantsPerDisease[,!(colnames(variantsPerDisease) %in% c("JAK2", "CALR"))]

		if(input$subTypeAll == "PMF & ET"){
			# Keep only PMF and ET subtypes
			variantsPerDisease = variantsPerDisease[rownames(variantsPerDisease) %in% c("ET-JAK2-wt-CALR-mutated", "ET-JAK2-mutated-CALR-wt",
				"PMF-JAK2-mutated-CALR-wt", "PMF-JAK2-wt-CALR-mutated"),]
		}

		# Test if at least 1 gene remains
		n = dim(variantsPerDisease)
		if(!n[2]) {print("No data to display.");return(list("Nothing to display.","Nothing to display."))}

		fisherList = fisherMutationMatrices(variantsPerDisease, input$alphaSubDisOc)
		genesToKeep = fisherList[[1]]
		diagnosesToKeep = fisherList[[2]]
		ORMat = fisherList[[3]]
		pvalMat = fisherList[[4]]
		
		# Test if at least 1 OR is left
		if(sum(genesToKeep) < 1) {print("No data to display.");return(list("Nothing to display.","Nothing to display."))}

		ORMat = ORMat[diagnosesToKeep, genesToKeep, drop=F] # Keep data.frame format even for single elements
		pvalMat = pvalMat[diagnosesToKeep, genesToKeep, drop=F]

		return(list(ORMat, pvalMat))
	})

	output$varSubDisOc <- renderPlotly({
		ORMat = subDisOcOR()[[1]]
	
		if(class(ORMat) != "data.frame"){return()} # No data to display

		heatmaply(ORMat, na.rm = T, colors = color.palette$function_bimod, show_grid = T,
			na.value=color.palette$function_bimod(256)[256], dendrogram = "none", key.title = "Odds-ratios",
			margins = c(75,200,NA,0), limits = c(0,20), grid_gap=2, label_names = c("Subtype", "Gene", "OR"), colorbar_len = 1)		
	})

	output$varSubDisOcORTable <- renderDataTable({
		tab = subDisOcOR()[[1]]
		return(cbind(diagnosis = rownames(tab), tab))
	})

	output$varSubDisOcPvalTable <- renderDataTable({
		tab = subDisOcOR()[[2]]
		return(cbind(diagnosis = rownames(tab), tab))
	})

	output$varSubDisOcPvalDL <- downloadHandler(
		filename = "variants_gene_subtype_pval.csv",
		content = function(file) {
			write.csv(subDisOcOR()[[2]], file)
		}
	)

	output$varSubDisOcORDL <- downloadHandler(
		filename = "variants_gene_subtype_OR.csv",
		content = function(file) {
			write.csv(subDisOcOR()[[1]], file)
		}
	)

	# Gene mutations per disease

	disOcOR <- reactive({
		dataset = filteredDataVariants()
		variantsPerDisease = table(dataset$diagnosis, dataset$GENESYMBOL)
		variantsPerDisease = variantsPerDisease[rowSums(variantsPerDisease) > 0, colSums(variantsPerDisease) > 0]
		
		# Remove genes with insufficient number of variants observed in the dataset
		variantsPerDisease = variantsPerDisease[,colSums(variantsPerDisease) >= input$nbRepDisOc]
	
		# Test if at least 1 gene remains
		n = dim(variantsPerDisease)
		if(!n[2]) {print("No data to display.");return(list("Nothing to display.","Nothing to display."))}

		fisherList = fisherMutationMatrices(variantsPerDisease, input$alphaDisOc)
		genesToKeep = fisherList[[1]]
		diagnosesToKeep = fisherList[[2]]
		ORMat = fisherList[[3]]
		pvalMat = fisherList[[4]]
		
		# Test if at least 1 OR is left
		if(sum(genesToKeep) < 1) {print("No data to display.");return(list("Nothing to display.","Nothing to display."))}

		ORMat = ORMat[diagnosesToKeep, genesToKeep, drop=F] # Keep data.frame format even for single elements
		pvalMat = pvalMat[diagnosesToKeep, genesToKeep, drop=F]

		return(list(ORMat, pvalMat))
	})

	output$varDisOc <- renderPlotly({
		ORMat = disOcOR()[[1]]

		if(class(ORMat) != "data.frame"){return()} # No data to display

		heatmaply(ORMat, na.rm = T, colors = color.palette$function_bimod, show_grid = T, key.title = "Odds-ratios",
			na.value=color.palette$function_bimod(256)[256], dendrogram = "none", margins = c(75,75,NA,0), limits = c(0,20), 
			grid_gap=2, label_names = c("Disease", "Gene", "OR"), colorbar_len = 1)		
	})

	output$varDisOcORTable <- renderDataTable({
		tab = disOcOR()[[1]]
		return(cbind(diagnosis = rownames(tab), tab))
	})

	output$varDisOcPvalTable <- renderDataTable({
		tab = disOcOR()[[2]]
		return(cbind(diagnosis = rownames(tab), tab))
	})

	output$varDisOcPvalDL <- downloadHandler(
		filename = "variants_gene_disease_pval.csv",
		content = function(file) {
			write.csv(disOcOR()[[2]], file)
		}
	)

	output$varDisOcORDL <- downloadHandler(
		filename = "variants_gene_disease_OR.csv",
		content = function(file) {
			write.csv(disOcOR()[[1]], file)
		}
	)


	# Variant per sample matrix
	output$varBinMat <- renderPlotly({
		dataset = filteredDataVariants()
		# Sort samples by diagnosis then by mutational load
		mutLoads = table(dataset$UNIQ_SAMPLE_ID)
		dataset$mutLoad = sapply(dataset$UNIQ_SAMPLE_ID, function(x) -mutLoads[names(mutLoads) == x][[1]])
		dataset$UNIQ_SAMPLE_ID = factor(dataset$UNIQ_SAMPLE_ID, levels = unique((dataset$UNIQ_SAMPLE_ID)[order(dataset$diagnosis, dataset$mutLoad)]))
		# Sort gene symbols by frequency
		dataset$GENESYMBOL = factor(dataset$GENESYMBOL, levels = names(sort(table(dataset$GENESYMBOL))))
		dataset = dataset[order(dataset$diagnosis),]

		gp1 = ggplot(dataset, aes(x = UNIQ_SAMPLE_ID, y = GENESYMBOL, fill = CADD_phred, text = diagnosis)) + 
			geom_raster() +
			labs(x = "Sample ID", y = "Gene Symbol", caption = "Hover to get values") +
			theme(axis.text.x = element_blank(), axis.text.y = element_blank())
		gpy1 = ggplotly(gp1) 
		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))
		
		gp2 = ggplot(dataset, aes(x = UNIQ_SAMPLE_ID, y = 1, fill = diagnosis)) + geom_raster() +
			labs(x = "Sample ID", y = "Diagnosis", caption = "Hover to get values") +
			scale_fill_manual(guide = guide_legend(title = NULL), values = color.palette$function_multi(length(levels(dataset$diagnosis)))) +
			theme(axis.text.x = element_blank(), axis.text.y = element_blank())
		gpy2 = ggplotly(gp2)

		subplot(gpy1, gpy2, nrows = 2, shareX = T, shareY = T, heights = c(0.9,0.1))
	})
}

# Start Shiny app
shinyApp(ui = shinyUi, server = shinyServer)
