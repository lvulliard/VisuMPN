table(substr(levels(as.factor(substr(data$unique.sample.id, 1, 4))), 1, 1))

gp1 = ggplot(dataCohort, aes(factor(sex), jak2.burden)) + geom_violin(aes(fill = factor(sex))) + geom_jitter(height = 0, width = 0.05)

gp2 = ggplot(dataCohort, aes_string(x = "sex", y = "jak2") + aes(text = unique.sample.id)) + geom_point()
ggplotly(gp2)

gp3 = ggplot(with(dataCohort, dataCohort[order(diagnosis, sex, genotype, thrombotic.events.after.diagn),])) + 
  geom_rect(aes(fill=diagnosis, ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=10, xmin=9, text = unique.sample.id)) +
  geom_rect(aes(fill=thrombotic.events.after.diagn, ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=9, xmin=7, text = unique.sample.id)) +
  geom_rect(aes(fill=sex, ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=7, xmin=4, text = unique.sample.id)) +
  geom_rect(aes(fill=genotype, ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=4, xmin=0, text = unique.sample.id)) +
  xlim(c(0, 10)) + 
  theme(aspect.ratio=1) + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(17)) +
  coord_polar(theta="y")  
ggplotly(gp3) 

		for(i in rev(orderVect)){
			d2 = d2[order(d2[,i]),]
		}

gp = ggplot(dataCohort)
for(i in 1:4){
	gp = append(gp, tail(gp, 1) + geom_rect(aes_string(fill=names(dataCohort)[qualitativeVar[i]]) +
		aes(ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=(i+1), xmin=(i)))) 
}			
gp = append(gp, tail(gp, 1) + xlim(c(1, 5)) + theme(aspect.ratio=1) + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(17)))
tail(gp, 1)

plotAllLayers<-function(){
  p<-ggplot(data=dataCohort)
  for(i in 1:4){ 
    p<-p+ geom_rect(aes_string(fill=names(dataCohort)[qualitativeVar[i]]) +
		aes_(ymax=1:dim(dataCohort)[1], ymin=1:dim(dataCohort)[1]-1, xmax=i+1, xmin=i)) 
  }
  return(p)
}

plotAllLayers()




		gp1 = ggplot(dataCohort, aes_string(x = "sex", y = "hb.at.diagn" )) + geom_violin(aes_string(fill = "sex")) + geom_jitter(height = 0, width = 0.1) + xlab("Sex") + ylab("Hb level") + theme(legend.position="none")
		margpy1 <- list(l=60, r=5, b=40, t=5) # margins on each side
		gpy1 = ggplotly(gp1 + theme_light()) %>% layout(margin=margpy1)
		gpy1$x$data[[i]]$text = gsub("^[^>]*>", "", gpy1$x$data[[i]]$text)
		gpy1$x$data[[j]]$text = paste0(gpy1$x$data[[j]]$text, paste0("<br />", dataCohort$unique.sample.id))
		
		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))



regmatches(gpy1$x$data[[i]]$text,gregexpr("(?<=:).*",foo,perl=TRUE))



# Fisher test

variantsPerPatient = table(dataVariants$GENESYMBOL, dataVariants$UNIQ_SAMPLE_ID) != 0
n = length(colnames(variantsPerPatient))

ctngVariants = variantsPerPatient %*% t(variantsPerPatient)
ctngVariants[c(2,35),c(2,35)]

n_JAK.AEBP = ctngVariants[2,35]
n_JAK = ctngVariants[35,35]
n_AEBP = ctngVariants[2,2]

table_JAK.AEBP = data.frame(AEBPp = c(n_JAK.AEBP, n_AEBP - n_JAK.AEBP), AEBPn = c(n_JAK - n_JAK.AEBP, n + n_JAK.AEBP - n_AEBP - n_JAK))
rownames(table_JAK.AEBP) = c("JAKp", "JAKn")

fisher.test(variantsPerPatient[2,], variantsPerPatient[35,])
fisher.test(table_JAK.AEBP)
fisher.test(t(table_JAK.AEBP))

OR = n_JAK.AEBP * (n + n_JAK.AEBP - n_JAK - n_AEBP) / ((n_JAK - n_JAK.AEBP)*(n_AEBP - n_JAK.AEBP))

p_val = 0
for (x in n_JAK.AEBP:min(n_JAK, n_AEBP)){
	table_JAK.AEBP = data.frame(AEBPp = c(x, n_AEBP - x), AEBPn = c(n_JAK - x, n + x - n_AEBP - n_JAK))
	p_val = p_val + choose(n_JAK, x)*choose(sum(table_JAK.AEBP[2,]),table_JAK.AEBP[2,1])/choose(n,n_AEBP)
} # Not just tables with higher interaction but tables with lower p-values

p_val = c()
for (x in 0:min(n_JAK, n_AEBP)){
	table_JAK.AEBP = data.frame(AEBPp = c(x, n_AEBP - x), AEBPn = c(n_JAK - x, n + x - n_AEBP - n_JAK))
	p_val = append(p_val, choose(n_JAK, x)*choose(sum(table_JAK.AEBP[2,]),table_JAK.AEBP[2,1])/choose(n,n_AEBP))
} 
p_val = sort(p_val)

n_JAK.AEBP = ctngVariants[2,35]
table_JAK.AEBP = data.frame(AEBPp = c(n_JAK.AEBP, n_AEBP - n_JAK.AEBP), AEBPn = c(n_JAK - n_JAK.AEBP, n + n_JAK.AEBP - n_AEBP - n_JAK))
rownames(table_JAK.AEBP) = c("JAKp", "JAKn")
p_val_obs = choose(n_JAK, n_JAK.AEBP)*choose(sum(table_JAK.AEBP[2,]),table_JAK.AEBP[2,1])/choose(n,n_AEBP)

sum(p_val[p_val <= p_val_obs]) # Real p-value

makeCtngTab = function(x, y){
	n_JAK.AEBP = ctngVariants[x,y]
	n_JAK = ctngVariants[y,y]
	n_AEBP = ctngVariants[x,x]

	return(data.frame(AEBPp = c(n_JAK.AEBP, n_AEBP - n_JAK.AEBP), AEBPn = c(n_JAK - n_JAK.AEBP, n + n_JAK.AEBP - n_AEBP - n_JAK)))
}


tVar = c()
tTab = c()
for(nbRep in 1:5000){
	t1 = Sys.time()
	f1 = apply(variantsPerPatient, 1, function(x) fisher.test(variantsPerPatient[2,], x))
	t2 = Sys.time()


	t3 = Sys.time()
	f2 = sapply(1:dim(ctngVariants)[1], function(x) fisher.test(makeCtngTab(x,2)) )
	t4 = Sys.time()

	tVar = append(tVar, as.numeric(t2-t1))
	tTab = append(tTab, as.numeric(t4-t3))
}

wilcox.test(tVar, tTab) # Using variantsPerPatient is faster

t5 = Sys.time()
L = apply(variantsPerPatient, 1, function(y) apply(variantsPerPatient, 1, function(x) fisher.test(y, x)))
t6 = Sys.time()


t7 = Sys.time()
nbGenes = dim(variantsPerPatient)[1]
L2 = sapply(1:nbGenes, function(y) sapply(y:nbGenes, function(x) fisher.test(variantsPerPatient[y,], variantsPerPatient[x,])))
t8 = Sys.time() # Should be less than 1s


# > all_pval[pvalOrder][2346]
# [1] 1
# > all_pval[pvalOrder][2346] == 1
# [1] FALSE
# > all_pval[pvalOrder][2346] >= 1
# [1] TRUE

 distfun: function used to compute the distance (dissimilarity) between
          both rows and columns. Defaults to dist.

hclustfun: function used to compute the hierarchical clustering when
          Rowv or Colv are not dendrograms. Defaults to hclust.


	variantsPerDisease = table(dataVariants$diagnosis, dataVariants$GENESYMBOL)

subsetDrivers = variantsPerPatient[rownames(variantsPerPatient) %in% c("JAK2", "CALR"),]
subsetDrivers[,colSums(subsetDrivers) == 2] #  H_0040A_CC2 H_0059A_TK are double mutants

heatmaply(B, limits = c(0,20), grid_gap=2, na.value="#ff0000", dendrogram=F, na.rm=T, 
	grid_color="#ffffff", label_names = c("Gene_r", "Gene_c", "OR"), fontsize_row = 15, fontsize_col = 15, colorbar_len = 0.6)

,
limits = c(0,20), grid_gap=2, label_names = c("Gene_r", "Gene_c", "OR"), colorbar_len = 0.7



colorRampPalette(c("#990033","#b3ecff",colorRampPalette(c("#ffffff", "#ffcc80", "#ffad33", "#ff9900"))(21) ))

colorbar_ypos
guide_colorbar(title = "Odds-ratio", 
	  UNIQ_SAMPLE_ID GENESYMBOL CHROM       POS ID
1        112c_JM       JAK2     9   5073770  .
2        112c_JM     NOTCH1     9 139418328  .
3        112c_JM        CBL    11 119148891  .
4        112c_JM       TP53    17   7577499  .
5        112c_JM      U2AF1    21  44514777  .

	
		mutLoads = table(dataset$UNIQ_SAMPLE_ID)
		dataset$mutLoad = sapply(dataset$UNIQ_SAMPLE_ID, function(x) -mutLoads[names(mutLoads) == x][[1]])
		dataset$UNIQ_SAMPLE_ID = factor(dataset$UNIQ_SAMPLE_ID, levels = unique((dataset$UNIQ_SAMPLE_ID)[order(dataset$diagnosis, dataset$mutLoad)]))
		# Sort gene symbols by frequency
		dataset$GENESYMBOL = factor(dataset$GENESYMBOL, levels = names(sort(table(dataset$GENESYMBOL))))
		dataset = dataset[order(dataset$diagnosis),]
		gp1 = ggplot(dataset, aes(x = UNIQ_SAMPLE_ID, y = GENESYMBOL, fill = CADD_phred)) + 
			geom_raster() +
			labs(x = "Sample ID", y = "Gene Symbol", caption = "Hover to get values") +
			theme(axis.text.x = element_blank(), axis.text.y = element_blank())
		gpy1 = ggplotly(gp1) 

		gp2 = ggplot(dataset, aes(x = UNIQ_SAMPLE_ID, y = 1, fill = diagnosis)) + geom_raster() +
			labs(x = "Sample ID", y = "Diagnosis", caption = "Hover to get values") +
			scale_fill_manual(guide = guide_legend(title = NULL), values = color.palette$function_multi(length(levels(dataset$diagnosis)))) +
			theme(axis.text.x = element_blank(), axis.text.y = element_blank())
		gpy2 = ggplotly(gp2)

		subplot(gpy1, gpy2, nrows = 2, shareX = T, shareY = T, heights = c(0.9,0.1))

		scale_x_continuous(breaks = 1:nbVarSelected + 0.5, labels = input$plotVariablesPie) +
		
 		style(gpy1, hoverinfo = "text", hoverlabel = list(bgcolor = color.palette$bg))

table(dataCohort[!is.na(dataCohort[,dt[[3]]]),dt2[[3]]])


dt = dataCohortTypes[13,] # Simulate dt and i increase to get gpy1
elementsLegend = sapply(gpy1$x$data, function(x) x$showlegend)
gpy1$x$data[elementsLegend]
sapply(gpy1$x$data[elementsLegend], function(x) x$name)
sapply(gpy1$x$data[elementsLegend], function(x) x$legendgroup)
sapply(gpy1$x$data[elementsLegend], function(x) {x = })

# Patient ID
pat = str_extract(dataAberrations$unique.sample.id, "^.\\d*")
# Visit
str_extract(dataAberrations$unique.sample.id, ".(?=#)")
# Batch
str_extract(dataAberrations$unique.sample.id, ".$")

A = levels(as.factor(dataAberrations$sample.id))
B = levels(as.factor(gsub("_[^_]*$", "", dataVariants$UNIQ_SAMPLE_ID)))
# All patients having variants have been tested for aberrations

table(dataAberrations$unique.sample.id %in% dataCohort$unique.sample.id)
table(dataAberrations$sample.id %in% dataCohort$sample.id)

dataAberrations$unique.sample.id.no.batch[which(!(dataAberrations$unique.sample.id.no.batch %in% dataCohort$unique.sample.id.no.batch))]

shinyUi <- fluidPage(mainPanel(p("test")))
shinyServer <- function(input, output) {}
shinyApp(shinyUi, shinyServer)

shinyUi <- fluidPage(mainPanel(
	plotOutput("plot1")
))
shinyServer <- function(input, output) {
	output$plot1 <- renderPlot({chorddiag(titanic.mat, type = "bipartite", 
		groupColors = groupColors,
		tickInterval = 50)
	})
}
shinyApp(shinyUi, shinyServer)


library(circlize)
library(migest)
library(dplyr)
 
### Make data
m <- data.frame(order = 1:6,
            country = c("Ausralia", "India", "China", "Japan", "Thailand", "Malaysia"),
            V3 = c(1, 150000, 90000, 180000, 15000, 10000),
            V4 = c(35000, 1, 10000, 12000, 25000, 8000),
            V5 = c(10000, 7000, 1, 40000, 5000, 4000),
            V6 = c(7000, 8000, 175000, 1, 11000, 18000),
            V7 = c(70000, 30000, 22000, 120000, 1, 40000),
            V8 = c(60000, 90000, 110000, 14000, 30000, 1),
            r = c(255,255,255,153,51,51),
            g = c(51, 153, 255, 255, 255, 255),
            b = c(51, 51, 51, 51, 51, 153),
            stringsAsFactors = FALSE)
df1 <- m[, c(1,2, 9:11)]
m <- m[,-(1:2)]/1e04
m <- as.matrix(m[,c(1:6)])
dimnames(m) <- list(orig = df1$country, dest = df1$country)
#Sort order of data.frame and matrix for plotting in circos
df1 <- arrange(df1, order)
df1$country <- factor(df1$country, levels = df1$country)
m <- m[levels(df1$country),levels(df1$country)]
 
 
### Define ranges of circos sectors and their colors (both of the sectors and the links)
df1$xmin <- 0
df1$xmax <- rowSums(m) + colSums(m)
n <- nrow(df1)
df1$rcol<-rgb(df1$r, df1$g, df1$b, max = 255)
df1$lcol<-rgb(df1$r, df1$g, df1$b, alpha=200, max = 255)
 
### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()
 
### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)
 
### Sector details
circos.initialize(factors = df1$country, xlim = cbind(df1$xmin, df1$xmax))
 
### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$country, track.height=0.1,
                      #panel.fun for each sector
                      panel.fun = function(x, y) {
                      #select details of current sector
                      name = get.cell.meta.data("sector.index")
                      i = get.cell.meta.data("sector.numeric.index")
                      xlim = get.cell.meta.data("xlim")
                      ylim = get.cell.meta.data("ylim")
 
                      #text direction (dd) and adjusmtents (aa)
                      theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                      dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                      aa = c(1, 0.5)
                      if(theta < 90 || theta > 270)  aa = c(0, 0.5)
 
                      #plot country labels
                      circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
 
                      #plot main sector
                      circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                  col = df1$rcol[i], border=df1$rcol[i])
 
                      #blank in part of main sector
                      circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, 
                                  col = "white", border = "white")
 
                      #white line all the way around
                      circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
 
                      #plot axis
                      circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5), 
                                  minor.ticks=1, labels.away.percentage = 0.15)
                    })
 
### Plot links (inner part)
### Add sum values to df1, marking the x-position of the first links
### out (sum1) and in (sum2). Updated for further links in loop below.
df1$sum1 <- colSums(m)
df1$sum2 <- numeric(n)
 
### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
           timevar="dest", time=rownames(m),  v.names = "m")
df2 <- arrange(df2,desc(m))
 
### Keep only the largest flows to avoid clutter
df2 <- subset(df2, m > quantile(m,0.6))
 
### Plot links
for(k in 1:nrow(df2)){
    #i,j reference of flow matrix
    i<-match(df2$orig[k],df1$country)
    j<-match(df2$dest[k],df1$country)
 
#plot link
circos.link(sector.index1=df1$country[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
            sector.index2=df1$country[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
            col = df1$lcol[i])
 
#update sum1 and sum2 for use when plotting the next link
df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
}



tracks = list()

addBioCircosSNPTrack <- function(tracklist, trackname, chromosomes, positions, values)

BioCircosTracklist <- function(){
	x = list()
	class(x) <- c("BioCircosTracklist", class(x))
  	return(x)
}