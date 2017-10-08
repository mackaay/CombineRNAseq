##create new repository 
##in this folder, create a folder named as 'data'
## unzip all files into data folder

# Read the data into R
small_counts <- read.table("data/small_counts.txt", header = TRUE)
print(small_counts)
dim(small_counts)


#subseting
small_counts$Sample_1
small_counts[, 1]
small_counts[, c("Sample_1")]
small_counts[, c("Sample_1", "Sample_3")]
small_counts[1:3, c("Sample_1", "Sample_3")]
#vectiorisation
small_counts$Sample_1 * 2
log(small_counts)
sum(small_counts$Sample_1)
sum(small_counts$Sample_2)
# Sum the counts for each sample
#Note the MARGIN argument tells apply which direction you want to go in. MARGIN = 1 will apply sum to each row, while MARGIN = 2 will apply sum to each column.
sample_sums = apply(small_counts, MARGIN = 2, sum)
print(sample_sums)

# data type
ResultsTable_small <- read.table("data/ResultsTable_small.txt", header=TRUE)
print(ResultsTable_small)

str(ResultsTable_small)
str(ResultsTable_small$SYMBOL)
typeof(ResultsTable_small$SYMBOL)

mydata <- c("case", "control", "control", "case")
factor_ordering_example <- factor(mydata, levels = c("control", "case"))
str(factor_ordering_example)
ResultsTable_small <- read.table("data/ResultsTable_small.txt", stringsAsFactors = FALSE, header=TRUE)

#sorting
ResultsTable_small[order(ResultsTable_small$logFC), ]
#subseting with logical statement
ResultsTable_small$logFC > 3
ResultsTable_small$logFC[ResultsTable_small$logFC > 3]
ResultsTable_small[ResultsTable_small$logFC > 3, ]
ResultsTable_small[ResultsTable_small$logFC > 3 | ResultsTable_small$logFC < -3, ]
ResultsTable_small[abs(ResultsTable_small$logFC) > 3, ]

my_genes <- c("Smad7", "Wif1", "Fam102b", "Tppp3")
ResultsTable_small$SYMBOL %in% my_genes
ResultsTable_small[ResultsTable_small$SYMBOL %in% my_genes, ]

match(my_genes, ResultsTable_small$SYMBOL)
ResultsTable_small[match(my_genes, ResultsTable_small$SYMBOL), ]


