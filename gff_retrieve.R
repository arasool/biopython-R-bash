#Ayesha Rasool
#this is part 4 of the midterm.
#the resultant file is called "final.gff", which is the output of this script/function inside.
setwd()
'''
Using the Tair10 database on orion.bio.nyu.edu, create a function that queries the
Tair10 database to produce the original GFF file that was used to load Chromosome 1
features. For the last column in the GFF file, you only need to retrieve the Name or
ID of the features.

'''

library("RMySQL")

# Assuming the host is not busy, which it was when this script was created:
# attribute = data.frame(dbGetQuery(con, "SELECT * from attribute"))
# feature = data.frame(dbGetQuery(con, "SELECT * from feature"))
# typelist = data.frame(dbGetQuery(con, "SELECT * from typelist"))
# locationlist = data.frame(dbGetQuery(con, "SELECT * from locationlist"))

# If host is busy and cannot connect to database, use local tables:
#attribute=read.csv(file="attribute.csv",sep=",", header=FALSE)
#feature=read.csv(file ="feature.csv", sep=",", header=FALSE)
#typelist=read.csv(file="typelist.csv", sep=",", header=FALSE)
#locationlist=read.csv(file="locationlist.csv", sep=",", header=FALSE)
sqlquery <- read.csv(file="sqlquery.csv", sep=";", header=TRUE)
sourcefeature <- read.csv(file="sourcefeature.csv", sep=";", header=TRUE)

makeGFFFile <- function() {
  
  require(reshape)
  library(stringr)
  con=dbConnect(MySQL(), user='ga1009', password='mkatari@nyu', dbname='tair10', host='orion.bio.nyu.edu')
  sqlquery <- data.frame(dbGetQuery(con, "SELECT DISTINCT locationlist.seqname, feature.start, feature.end, CASE feature.strand WHEN -1 THEN '-' WHEN 0 THEN '.' WHEN 1 THEN '+' END AS strand, attribute.id FROM feature, attribute, locationlist WHERE locationlist.seqname = 'Chr1'AND feature.id = attribute.id AND feature.seqid = locationlist.id"))
  sourcefeature <- data.frame(dbGetQuery(con, "SELECT DISTINCT attribute.id, typelist.tag FROM attribute, locationlist, typelist, feature WHERE locationlist.seqname = 'Chr1'AND locationlist.id = feature.seqid AND typelist.id = feature.typeid AND feature.id = attribute.id;"))
  
  # Split the sourcefeature on the : 
  #df = data.frame(transform(sourcefeature, tag = colsplit(tag, split = ":", names = c('feature', 'source'))))
  df <- data.frame(str_split_fixed(sourcefeature$tag, ":", 2))
  
  # Create data.frame of dots:
  dots <- data.frame(rep(".", nrow(sqlquery)))
  
  # Put everything together:
  finaldata <- data.frame(cbind(sqlquery[1], df[2], df[1], sqlquery[2], sqlquery[3], dots[1], sqlquery[4], dots[1], sqlquery[5]))
  colnames(finaldata) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "annotation")
  
  # Output file:
  write.table(finaldata, "final.gff", row.names=FALSE)
  
}

makeGFFFile()

# Clean:
rm(sqlquery)
rm(sourcefeature)
rm(df)
rm(dots)
rm(finaldata)
