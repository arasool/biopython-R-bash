#The original data consisted of microarrays from 128 different individuals with acute
#lymphoblastic leukemia (ALL). We will look at subtype BCR/ABL versus control samples
#("NEG").
#Do the following exploratory data analyses in R. Submit your R script with comments and
#outputs.
#1) do a hierarchical clustering using Euclidean distance and complete linkage and plot your
#result. (7 marks)
#2) show a heatmap of the data (6 marks)
#3) do a kmeans clustering with k=5. (7 marks)
#4) perform a PCA. Use the pairs() function to compare all principal components against each
#other and show the results. (10 marks)
#5) perform a GSEA on the expression. Show the top 5 significantly enriched pathways. (20
#marks)



load("Data_mining_ass1-2.RData")
library("ISLR")
library("gplots")

# Apply ALL.labs names to ALL.data
colnames(ALL.data)<-ALL.labs
head(ALL.data) 

# START OF PART 1: Hierarchical Clustering
#going to use a subset of data because of the warning I receieved by using all of the data.
# First 50 rows(genes), all columns(individuals):
ALL.dataSUBSET <- as.matrix(ALL.data[1:50,])

# do hierarchical clustering using complete linkage
#can be tested on the ALL.data as well. will see individuals clustered.
hc.complete=hclust(dist(ALL.dataSUBSET, method="euclidean"), method = "complete")

par(mfrow=c(1,1)) 
# plot dendrogram  (note that R overloads the plot function)
plot(hc.complete,main="Complete Linkage", xlab="", sub="", cex=.9)

# END OF PART 1

# START OF PART 2: Heat Map 
#individual dendrograms on the side, genes dendrograms on the top.
hv <-heatmap(t(ALL.dataSUBSET), trace="none",dendrogram="column",keysize = 1,scale = c("none"),symkey=F,symm=F,cexRow=1, cexCol=1, sepcolor="black")

# END OF PART 2

# START OF PART 3: K-Means

hc.out=hclust(dist(ALL.data))
hc.clusters=cutree(hc.out,5)  # cut dendrogram at height to get 5 clusters

set.seed(2)
km.out=kmeans(scale(ALL.data), 5)
km.clusters=km.out$cluster
table(km.clusters,hc.clusters)

# END OF PART 3

# START OF PART 4: PCA
#prcomp -> pass in the matrix of the data, scale the data where variance is set to 1.
pr.out=prcomp(ALL.data, scale=TRUE)

#assign random colors to various kinds of data we will visualize.
Cols=function(vec){
  cols=rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

par(mfrow=c(1,1)) # sub-windows for plots.
# plot the principal components score vectors 
#pr.out$x variable to save PCA results, they are indexed.
#plotting all the principal components, varied in different colors.Not sure why they are not in diagonal formation.
plot(pr.out$x, col=Cols(ALL.data), pch=19,xlab="NEG",ylab="BCR/ABL")
summary(pr.out)   # summary of the PCA output gives variance explained
plot(pr.out)
pve=100*pr.out$sdev^2/sum(pr.out$sdev^2) # percentage of variance.
par(mfrow=c(1,1))
plot(pve,  type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="brown3")

#Use the pairs() function to compare first 5 principal components against each other.
#pr.out is the variable where they are stored, so I am trying to extract and compare the first 5.
#the pdf with the plots will appear in your working directory.
pdf("pairsAll.pdf", width=20, height=20)
pairs(pr.out$x[,c(1,2,3,4,5)])
dev.off()



# END OF PART 4

# START OF PART 5

library(gage)
install.packages(gageData)
library(gageData)
df=as.matrix(ALL.data)
data(kegg.gs)
gage(df,gsets=kegg.gs)

#this gage set produces a chart with statistical values such as mean of BCR/ABL and NEG

#calculating p-value, and obtaining hsa values with location/function of the gene
ALL.data.kegg.p <- gage(ALL.data, gsets = kegg.gs)

head(ALL.data.kegg.p$greater[, 1:5], 5)
head(ALL.data.kegg.p$less[, 1:5], 5)

#significantly up-regulated and down-regulated gene-sets
ALL.data.kegg.sig<-sigGeneSet(ALL.data.kegg.p, outname="ALL.data.kegg")


# END OF PART 5