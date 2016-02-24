###################
#   WGNCA
###################

#### Data input, cleaning and pre-processing 

# they suggest doing varianceStabilizing in Deseq before inputting data
# for now lets log transform the RPKMs

data = read.table("GSE77166_RPKM.txt", sep ='\t', stringsAsFactors = F, header = T)
head(data)
dimen = dim(data)[2]

loggy = function (data){
  x = log2(data + 1)
  return(x)
}


logdat = loggy(data[,3:dimen])
head(logdat)

logdatfull = cbind(data[,1:2], logdat)
head(logdatfull)

### the dataframe has zeros where the genes arent annotated 
View(table(logdatfull$genename.chicken))

## replace 0s with x1,x2,...
sum(logdatfull$genename.chicken == 0)
# 2170 are equal to 0


for (y in 1:nrow(logdatfull))
{
  if (logdatfull[y,2] == 0)
  {
    logdatfull[y,2] = paste("x", y, sep = "")
  }
}

freqs = as.data.frame(table(logdatfull$genename.chicken))

## take the mean of genes repeats 
sub1 = subset(freqs, freqs[,2] ==1)
uniques = subset(logdatfull, logdatfull$genename.chicken %in% sub1[,1])
rownames(uniques) = uniques[,2]
uniques = uniques[,-c(1,2)]
head(uniques)
sub2 = subset(freqs, freqs[,2] > 1)

for (y in 1:nrow(sub2))
     {
       name = as.character(sub2[y,1])
       sub3 = subset(logdatfull, logdatfull$genename.chicken == name)
       meanie = as.data.frame(t(colMeans(sub3[,3:ncol(sub3)])))
       rownames(meanie) = name
       if (y == 1) {
         nonuniques = (meanie)
       }
       if (y > 1) {
         nonuniques = as.data.frame(rbind(nonuniques, meanie))
       }
       
}

dim(nonuniques)

curated = rbind(uniques, nonuniques)

### now lets try WGCNA

source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")
library(WGCNA)


# check if file is okay
gsg = goodSamplesGenes(curated, verbose = 3)
gsg$allOK
# theyre all okay

# now we are clustering the samples using hierachical clustering
# doing this to detect outliers
tree = hclust(dist(curated), method = "average")

# plotting time 

plot(tree, main = "Sample clustering to detect outliers")
# plot cutoff line
abline(h = 15, col = "red");
# everything is under the cutoff 

### create data expression table 

expressed = as.data.frame(t(curated))

##### creating trait table 
chickennames = rows