source("http://bioconductor.org/biocLite.R")
biocLite("limma")

library(edgeR)
library(limma)
library(EDASeq)
library(Glimma)


# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
targets <- read.delim("data/targets.txt")

head(seqdata)
dim(seqdata)
targets

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)

# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)

colnames(countdata)
table(colnames(countdata)==targets$SampleName)

#There are a few ways to filter out lowly expressed genes. When there are biological replicates in each group, in this case we have a sample size of 2 in each group, we favour filtering on a minimum counts per million threshold present in at least 2 samples. Two represents the smallest sample size for each group in our experiment. In this dataset, we choose to retain genes if they are expressed at a counts-per-million (CPM) above 0.5 in at least two samples.

# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)

#A CPM of 0.5 is used as it corresponds to a count of 10-15 for the library sizes in this data set. If the count is any smaller, it is considered to be very low, indicating that the associated gene is not expressed in that sample. A requirement for expression in two or more libraries is used as each group contains two replicates. This ensures that a gene will be retained if it is only expressed in one group. Smaller CPM thresholds are usually appropriate for larger libraries. As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5. You should filter with CPMs rather than filtering on the counts directly, as the latter does not account for differences in library sizes between samples.

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])

# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)


y <- DGEList(counts.keep)
# have a look at y
y

# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples

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

#From the boxplots we see that overall the density distributions of raw log-intensities are not identical but still not very different. If a sample is really far above or below the blue horizontal line we may need to investigate that sample further.

# You need to give plotRLE a matrix of counts
plotRLE(y$counts,las=2)
title("RLE plot of unnormalised counts")

#Ideally, the distributions should be centred around zero and be as tight as possible. It looks a bit different to the boxplot of the logcounts, but overall the patterns look the same.

#An MDSplot is a visualisation of a principle components analysis, which determines the greatest sources of variation in the data. A principle components analysis is an example of an unsupervised analysis, where we don't need to specify the groups. If your experiment is well controlled and has worked well, what we hope to see is that the greatest sources of variation in the data are the treatments/groups we are interested in.
plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(targets$CellType)

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[targets$CellType]
col.cell

# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("bottomleft",fill=c("purple","orange"),legend=levels(targets$CellType))
# Add a title
title("Cell type")

# Similarly for status
levels(targets$Status)
col.status <- c("blue","red","dark green")[targets$Status]
col.status

plotMDS(y,col=col.status)
legend("bottomleft",fill=c("blue","red","dark green"),legend=levels(targets$Status),cex=0.8)
title("Status")

# Let's combine status and cell type info on one plot using plotting characters and colours
# Make a new group variable that joins together cell type and status with the paste command
par(mfrow=c(1,1))
group <- factor(paste(targets$CellType,targets$Status,sep="."))
group
# the table function will tell us how many samples we have in each group
table(group)

points <- c(0,1,2,15,16,17)
colors <- rep(c("blue", "darkgreen", "red"), 2)
plotMDS(y, col=colors[group], pch=points[group])
legend("bottomleft", legend=levels(group), pch=points, col=colors, ncol=2,cex=0.7)

#plotMDS plots the first two dimensions as a default, however you can plot higher dimensions using the dim argument.
plotMDS(y,dim=c(3,4),col=colors[group], pch=points[group])
legend("topright", legend=levels(group), pch=points, col=colors, ncol=2,cex=0.8)

pca <- prcomp(logcounts)
plot(pca$sdev)

labels <- paste(targets$SampleName, targets$CellType, targets$Status)
glMDSPlot(y, labels=labels, groups=group, folder="mds")

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

library(gplots)
library(RColorBrewer)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()

# Apply normalisation to DGEList object
y <- calcNormFactors(y)

library(limma)
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")

# Look at group variable again
group

# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design

library(stringi)
