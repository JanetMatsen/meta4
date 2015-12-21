library(DESeq2)

# Get all arguments after --args
args <- commandArgs(trailingOnly=TRUE)
print(args)


# Load in the un-normalized table with all samples and all organisms.  
tsvFile <- args[1]  
# we will restrict to normalizing all rows of the 88 samples for 1 genome at a time. 
genome <- args[2] 

print(tsvFile)
print(genome)

masterD <- read.table(tsvFile, sep="\t", header=T, quote="", row.names=2)

# drop some columns in masterD and d
# second locus_tag column came from merging in SQL twice: SQL has a 61 table limit! 
masterD$locus_tag.2 <- NULL
# Select just the data for the genome of interest. 
masterD <- masterD[masterD$genome == genome, ]
countData <- masterD
# Delete the genome and product columns; they aren't read counts and DESeq doesn't want them. 
countData$genome <- NULL
countData$product <- NULL
head(countData)
# remove all 0 rows
countData <- countData[rowSums(countData[, -1]) > 0, ]
head(countData)
head(colData)


# Read in the experiment info; imperative for DESeq's normalization scheme. 
# ---------FIX: added sample_info to current dir for development on local computer-----------------
colData <- read.table("../sample_info.xls", sep="\t", header=T, quote="", row.names=1)
#colData <- read.table("./sample_info.xls", sep="\t", header=T, quote="", row.names=1)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ week + O2)
# remove rows with all zeros before applying DESeq. 
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds <- DESeq(dds)
#res <- results(dds)
rld <- rlog(dds)
head(assay(rld), 3)
vsd <- varianceStabilizingTransformation(dds)
head(assay(vsd), 3)
normCounts <- assay(vsd)

geneInfo <- masterD[c("genome", "product")]
# join geneinfo with normCounts
