library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt")
head(seqdata)
dim(seqdata)

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)
head(countdata)

##filtering lowly counts
# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
table(rowSums(thresh))# There are 11433 genes that have TRUEs in all 12 samples.
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)

#Challenge
plot(myCPM[,2],countdata[,2])
plot(myCPM[,2],countdata[,2],ylim=c(0,50),xlim=c(0,3))
abline(v=0.5, col='blue')
abline(h=10, col='blue')

##Convert counts to DGEList object
y <- DGEList(counts.keep)
# have a look at y
y
# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples

##Quality control
#library sizes and distribution plots
y$samples$lib.size
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=2)
# Add a title to the plot
title("Barplot of library sizes")
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#Multidimensional scaling plots
plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$CellType)
# Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)
# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(sampleinfo$Status)
col.status <- c("blue","red","dark green")[sampleinfo$Status]
col.status
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

# There is a sample info corrected file in your data directory
# Old sampleinfo
sampleinfo
# I'm going to write over the sampleinfo object with the corrected sample info
sampleinfo <- read.delim("data/SampleInfo_Corrected.txt")
sampleinfo
# Redo the MDSplot with corrected information
par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","dark green")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$CellType))
title("Cell type")
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","dark green"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")
# Dimension 3 appears to separate pregnant samples from the rest. Dim4?
plotMDS(y,dim=c(3,4),col=col.status,pch=16,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col=col.status,pch=16)
legend("bottomright",legend=levels(sampleinfo$CellType),pch=c(1,4))
#Glimma
labels <- paste(sampleinfo$SampleName, sampleinfo$CellType, sampleinfo$Status)
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group <- factor(group)
glMDSPlot(y, labels=labels, groups=group, folder="mds")


##Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[sampleinfo$CellType]
# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()


##Normalisation for composition bias
# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")

# save(group,y,logcounts,sampleinfo,file="day1objects.Rdata")
# load("day1objects.Rdata")
# objects()
# Look at group variable again
group
# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v
# What is contained in this object?
names(v)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

# Fit the linear model
#fit a linear model for each gene using the lmFit function in limma
fit <- lmFit(v)
names(fit)
# make the comparison of interest
cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate,levels=design)
cont.matrix
#apply the contrasts matrix to the fit object to get the statistics and estimated parameters of our comparison that we are interested in. 
fit.cont <- contrasts.fit(fit, cont.matrix)
#call the eBayes function, which performs empirical Bayes shrinkage on the variances, and estimates moderated t-statistics and the associated p-values.
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
# limma decideTests function to generate a quick summary of DE genes for the contrasts.
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
#output the top 10 genes 
topTable(fit.cont,coef="B.PregVsLac",sort.by="p")
# This will give the same output
topTable(fit.cont,coef=1,sort.by="p")

#Adding annotation and saving the results
source("https://bioconductor.org/biocLite.R")
biocLite("org.Mm.eg.db")
columns(org.Mm.eg.db)
ann <- select(org.Mm.eg.db,keys=rownames(fit.cont),columns=c("ENTREZID","SYMBOL","GENENAME"))
# Have a look at the annotation
head(ann)
# double check that the ENTREZID column matches exactly to our fit.cont rownames.
table(ann$ENTREZID==rownames(fit.cont))

fit.cont$genes <- ann
topTable(fit.cont,coef="B.PregVsLac",sort.by="p")
limma.res <- topTable(fit.cont,coef="B.PregVsLac",sort.by="p",n="Inf")
write.csv(limma.res,file="B.PregVsLacResults.csv",row.names=FALSE)

# Plots after testing for DE
# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.PregVsLac"])

# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL)

par(mfrow=c(1,3))
# Let's look at the first gene in the topTable, Wif1, which has a rowname 24117
stripchart(v$E["24117",]~group)
# This plot is ugly, let's make it better
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's use nicer colours
nice.col <- brewer.pal(6,name="Dark2")
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="Wif1")

# Testing relative to a threshold (TREAT)
# Let's decide that we are only interested in genes that have a absolute logFC of 1.
# This corresponds to a fold change of 2, or 0.5 (i.e. double or half).
# We can perform a treat analysis which ranks our genes according to p-value AND logFC.
# This is easy to do after our analysis, we just give the treat function the fit.cont object and specify our cut-off.
fit.treat <- treat(fit.cont,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)
# Notice that much fewer genes are highlighted in the MAplot
par(mfrow=c(1,2))
plotMD(fit.treat,coef=1,status=res.treat[,"B.PregVsLac"])
abline(h=0,col="grey")
plotMD(fit.treat,coef=2,status=res.treat[,"L.PregVsLac"])
abline(h=0,col="grey")



## Gene Set Testing
#Gene ontology testing with goana
go <- goana(fit.cont, coef="B.PregVsLac",species = "Mm")
topGO(go, n=10)
colnames(seqdata)
m <- match(rownames(fit.cont),seqdata$EntrezGeneID)
gene_length <- seqdata$Length[m]
head(gene_length)
# Rerun goana with gene length information
go_length <- goana(fit.cont,coef="B.PregVsLac",species="Mm",covariate=gene_length)
topGO(go_length, n=10)

#CAMERA gene set testing using the Broad’s curated gene sets
# Load in the mouse c2 gene sets
# The R object is called Mm.c2
load("data/mouse_c2_v5.rdata")
# Have a look at the first few gene sets
names(Mm.c2)[1:5]
# Number of gene sets in C2
length(Mm.c2)
#map the Entrez gene ids between the list of gene sets and our voom object. 
c2.ind <- ids2indices(Mm.c2, rownames(v))
gst.camera <- camera(v,index=c2.ind,design=design,contrast = cont.matrix[,1],inter.gene.cor=0.05)
gst.camera[1:5,]
table(gst.camera$FDR < 0.05)
write.csv(gst.camera,file="gst_BPregVsLac.csv")

#ROAST gene set testing
# ROAST doesn’t care about what the other genes in the experiment are doing, which is different to camera and goana.
H.camera[1:10,]
grep("MYC_",names(c2.ind))
# Let's save these so that we can subset c2.ind to test all gene sets with MYC in the name
myc <- grep("MYC_",names(c2.ind))
# What are these pathways called?
names(c2.ind)[myc]
myc.rst <- roast(v,index=c2.ind[myc],design=design,contrast=cont.matrix[,1],nrot=999)
myc.rst[1:15,]

#Visualising gene set tests: Barcode and enrichment plots
# Have a look at the logFCs and t-statistics in fit.cont
names(fit.cont)
head(fit.cont$coefficients)
head(fit.cont$t)
par(mfrow=c(1,1))
# barcode plot with logFCs
barcodeplot(fit.cont$coeff[,1], index=c2.ind[["MENSSEN_MYC_TARGETS"]], 
            main="LogFC: MENSSEN_MYC_TARGETS")
# barcode plot using t-statistics
barcodeplot(fit.cont$t[,1], index=c2.ind[["MENSSEN_MYC_TARGETS"]], 
            main="T-statistic: MENSSEN_MYC_TARGETS")





