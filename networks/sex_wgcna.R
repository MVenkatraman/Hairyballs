###################
#   WGNCA
###################

#### Data input, cleaning and pre-processing ####


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

### WGCNA

library(WGCNA)

# check if file is okay
gsg = goodSamplesGenes(curated, verbose = 3)
gsg$allOK
# theyre all okay

### create data expression table 

expressed = as.data.frame(t(curated))

# now we are clustering the samples using hierachical clustering
# doing this to detect outliers
tree = hclust(dist(expressed), method = "average")

# plotting time 

plot(tree, main = "Sample clustering to detect outliers")
# plot cutoff line
abline(h = 15, col = "red");
# everything is under the cutoff 


##### creating trait table for sexes
traits = as.data.frame(list(1:8))
rownames(traits) = rownames(expressed)
colnames(traits) = "Population"
traits[,1] = c(1, 1, 2, 1, 2, 2, 1, 2)
#F F	M	F	M	M	F	M


### now re cluster 
tree2 = hclust(dist(expressed), method = "average")
traitColors = numbers2colors(traits[,1], signed = FALSE);

plotDendroAndColors(tree2, traitColors, groupLabels = names(traits), main = "Sample dendrogram and trait heatmap")

### save data
save(expression, traits, file = "chicken_input.RData")

##### Module construction #####
###############################


###### choose a soft threshold power ######
# the threshold for co-expression similarity

# first pick a set of soft t powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(expressed, powerVector = powers, verbose = 5)

# now we want to plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
# we're adding a cutoff line
abline(h=0.90,col="blue")


# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## so i choose the value 16 b/c it is the lowest value for which the scalefree top
# reaches 0.9

## now calculate adjacencies using the soft t power of 16

softPower =16
adjacency = adjacency(expressed, power = softPower)

# to minimize noise we convert adjacency to topographical overlap network. 
TOM = TOMsimilarity(adjacency);
# now we calculate dissimilarity
dissTOM = 1-TOM

### clustering the tom
# WGCNA masks hclust with a faster hclust
geneTree = hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
# hang = The fraction of the plot height by which labels should hang below the rest of the plot

## so we have a hanging tree that has densely connected areas that are correlated genes
# we need to cut it up to only get the correlated genes

#### We use dynamic tree cut

# lets choose an big sized module = 30
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

table(dynamicMods)
# it gave me 72 modules, the above function tells you the size of the modules
# label 0 is unassinged genes 

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# lets plot dendogram with colors
par(mfrow=c(1,1))
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# literally the coolest thing ever

### now lets merge modules whose expression profile is similar
# calculating eigengenes, which is kinda just pc1
MEList = moduleEigengenes(expressed, colors = dynamicColors)
MEs = MEList$eigengenes

# calculate dissimilarity of modules
MEDiss = 1-cor(MEs)

# now we cluster the eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# we need to pick a correlation threshold
# we pick .75

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")

## merging correlated modules
merge = mergeCloseModules(expressed, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# give the merged modules colors
mergedColors = merge$colors

# Eigengenes of the new merg
mergedMEs = merge$newMEs

### so lets see what happened with the merging, we plot geneTree again
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# renames mergecolors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "networkconstruct.RData")

##### Trait Correlation ######
###############################

nGenes = ncol(expressed)
nSamples = nrow(expressed)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(expressed, moduleColors)$eigengenes
# correlate trait and eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### Gene Significance and Module Membership
# GS = absolute correlation between gene and trait
# MM = the correlation of the module eigengene and the gene expression profile

# Define variable weight containing the weight column of datTrait
pop = as.data.frame(traits$Population)
names(pop) = "pop"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

# calculate MM
geneModuleMembership = as.data.frame(cor(expressed, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# calculate GS
geneTraitSignificance = as.data.frame(cor(expressed, pop, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(pop), sep="");
names(GSPvalue) = paste("p.GS.", names(pop), sep="")

### identify genes with high GS and MM
## we will look at the magenta module

module = "plum1"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for population",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(expressed)[moduleColors=="magenta"]


######### export data
genenames = names(expressed)

geneInfo0 = data.frame(moduleColor = moduleColors, 
                       geneTraitSignificance, GSPvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, pop, use = "p")))

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.pop));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")

#######################
# Subset the modules 
#######################

# get the modules that have genes that have GS p values < 0.05

head(geneInfo)
dimdim = dim(geneInfo)

sub_module = subset(geneInfo$moduleColor, geneInfo$p.GS.pop < 0.05)
mods = unique(sub_module)
length(mods)

important = subset(geneInfo, geneInfo$moduleColor %in% mods)

### now subset this set by modules that show strong module-trait relationships

moddy = as.data.frame(moduleTraitCor)
modtraithigh = subset(rownames(moddy) , moddy[,1] >= 0.5 | moddy[,1] <= -0.5)

# get rid of the ME part
modtraithigh = gsub("ME","",modtraithigh)

finalsub = subset(important, important$moduleColor %in% modtraithigh)
finalmods = unique(finalsub$moduleColor)
finalmods


# look at the size of the modules
length(subset(finalsub$moduleColor, finalsub$moduleColor == "plum1"))

##################################
# Extract data for analysis
##################################

# finalmods is the list of modules we want

# extract subset by color (i) and top hub genes (j). and get rid of self correlation
# and drop low weights (k)

subsetTOM = function(i, j, k=0){
  damod = as.character(i)
  inModule = (moduleColors==damod)
  danames = names(expressed)
  modjames = danames[inModule]
  modTOM = TOM[inModule, inModule]
  head(modTOM)
  dimnames(modTOM) = list(modjames, modjames)
  nTop = j
  IMConn = softConnectivity(expressed[, modjames])
  top = (rank(-IMConn) <= nTop)
  modTOMTOM = modTOM[top,top]
  tomtom <- as.data.frame(as.table(modTOMTOM))
  colnames(tomtom) = c("Gene1", "Gene2", "Weights")
  tomtom = tomtom[!(tomtom$Gene1 == tomtom$Gene2), ]
  tomtom = tomtom[!(tomtom$Weight < k), ]
  return(tomtom)
}


##################################
# Finding the correct threshold 
##################################
# k in subsetTOM

darkseagreen4 <- subsetTOM("darkseagreen4", 100)

hist(darkseagreen4$Weight, breaks = 30 , col = "salmon", axes = F)
axis(1, at = seq(from = 0, to = 1, by = 0.01), tick = T)


# extract 
darkseagreen4 <- subsetTOM("darkseagreen4", 30, 0.04)
View(darkseagreen4)

write.table(darkseagreen4, file = "darkseagreen4full.txt", sep = "\t", row.names = FALSE, quote = FALSE)
