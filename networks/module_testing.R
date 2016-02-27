#################################
#     Module Testing 
#################################
# Madhvi Venkatraman
# 2/27/2015

sexmods= read.table("sex_modules.txt", sep ='\t', stringsAsFactors = F, header = T)
row.names(sexmods) = sexmods$X
d = dim(sexmods)
sexmods = sexmods[2:d[2]]
View(sexmods)

par(mfrow=c(2,2))
pca1 <- prcomp(sexmods[, 1:4])
summary(pca1)
plot(pca1$x[,1], pca1$x[,2], col = 3:7, pch = 19, xlab = "PC1", ylab = "PC2", main = "dark red", cex = 2)

pca2 <- prcomp(sexmods[, 5:8])
summary(pca2)
plot(pca2$x[,1], pca2$x[,2], col = 3:7, pch = 19, xlab = "PC1", ylab = "PC2", main = "sea green", cex = 2)

pca3 <- prcomp(sexmods[, 17:20])
summary(pca3)
plot(pca3$x[,1], pca3$x[,2], col = 3:7, pch = 19, xlab = "PC1", ylab = "PC2", main = "magenta", cex = 2)

pca4 <- prcomp(sexmods[, 21:24])
summary(pca4)
plot(pca4$x[,1], pca4$x[,2], col = 3:7, pch = 19, xlab = "PC1", ylab = "PC2", main = "maroon", , cex = 2)
