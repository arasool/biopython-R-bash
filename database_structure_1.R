#Ayesha Rasool HW # 3
#two files for this HW:
#hw3_AR.R, expvalue.sqlite
# this work is done on locally installed R
# Set working directory:
setwd("/Users/AR/Desktop/BBDHW")
getwd()

'''Part 1: Loading Expression data into the Database. (15)

The file expvalues.txt is an example of how you may find normalized expression data. Each row represents a probe placed on the microarray chip (hopefully representing a gene uniquely) and each column is a different experiement.

The goal of this homework is to create a database containing the information represented in the file and then use the database to make queries. 

First load the file expvalues.txt into your R or python workspace. Now write a function in R or python that will take a dataframe(pandas dataframe in python) of experimental values and load them into the following tables:

Experiment
	expid primary key
	expname

Probes
	probeid primary key
	probename

Data
	dataid primary key
	expid foreign key
	probeid foreign key
	expvalue

It is a good idea to separate the problem into parts and tackle each part separately. Creating functions for the different tasks will help you in the future homeworks, the midterm, and the final project so don’t be lazy ... do it right.

Part 2: Querying the database to calculate average expression (5)
Write a function in R or Python that calculates average value for each probeset across all experiments. For example the average value of 244901_at in all experiments is 224.8. 

There are many ways to do this but you MUST query the database for the experimental values.
'''

#make sure to have expvalues.tab3.txt in the wd.

# Read data from file into a data frame and assign it to a variable
df1=read.table(file="expvalues.tab3.txt", sep="\t", header=TRUE)


# Probe:
#this table will refer to the first column of expvalues
probe=unique(df1["X"])
#naming the column that contains probe names
colnames(probe) <- "probename"
#generating probeid
probe$probeid=seq.int(nrow(probe))
#switch the order of columns
probe <- probe[c(2,1)]

# Experiment:
# this table refers to columns 2 - 7 of expvalues
experiment <- matrix(colnames(df1[2:7]),ncol=1,byrow=TRUE)
#converting the list into dataframe
experiment <- data.frame(matrix(unlist(experiment),ncol=1,byrow=TRUE))
#naming the column that contains experiment names
colnames(experiment) <- "expname"
#generating expid
experiment$expid=seq.int(nrow(experiment))
#switch the order of columns
experiment <- experiment[c(2,1)]

#installing "reshape" package, and use it for flexibly restructuring and aggregating data.(online) 
install.packages("reshape")
library(reshape)

#Data
#going to use melt function on the expvalues.txt file.
#The melt function takes data in wide format and stacks a set of columns into a single column of data (online).
meltdf=melt(df1,id="X")
#melt provides the columns with default names, I am renaming them based on their reference I made in my tables.
colnames(meltdf)=c("probename","expname","expvalue")
#we are going to "join" or merge tables based on probename and expname
data=merge(merge(meltdf,probe, by ="probename"),experiment, by ="expname")[3:5]
#generating the dataid
data$dataid=seq.int(nrow(data))
#rearranging the table into the desired order
data <- data[c(4, 3, 2, 1)]

#now sort the data table by probeid
sort.data.frame <- function(x, decreasing=FALSE, by=1, ... ){
  f <- function(...) order(...,decreasing=decreasing)
  i <- do.call(f,x[by])
  x[i,,drop=FALSE]
}
data <- sort(data, by="probeid")
#row.names column appeared, to refer to actual row number before sort, going to remove it.
row.names(data)<-NULL


#installing RSQLite package
install.packages("RSQLite")
#connecting to SQLITE and creating a database expvalue.sqlite
drv<-dbDriver("SQLite")
con<-dbConnect(drv, "expvalue.sqlite")
#listing the tables in the database, there should be nothing in there.
dbListTables(con)
#start adding or writing the tables that we created to the database.
dbWriteTable(con, "probe", probe)
dbWriteTable(con, "experiment", experiment)
dbWriteTable(con, "data", data)
#checking to make sure all tables are present
dbListTables(con)
dbListFields(con, "experiment")
#storing my query of calculating the average in a variable called AVERAGE
AVERAGE<-dbGetQuery(con, "SELECT probename, AVG(expvalue) as AVERAGE
FROM probe, data
WHERE probe.probeid = data.probeid
GROUP BY data.probeid")


#the following two database reference are not needed anymore, ignore these commands.
remove(meltdf)
remove(df1)