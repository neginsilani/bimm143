#' ---
#' title: "Class 5: R graphics"
#' author: "Negin Silani"
#' date: "April 16, 2019"
#' ---

#class 5 R graphics

# 2A. Line plot
weight<- read.table("bimm143_05_rstats/weight_chart.txt",
                    header = TRUE)

plot( weight$Age, weight$Weight, xlab="Age(month)", 
      ylab="Weight (kg)", pch=18, typ="b",
      main="Some title")


# 2B. Barplot    
feat<- read.table("bimm143_05_rstats/feature_counts.txt", 
           sep = "\t",
           header= TRUE)

 
# I need to argue with this plot to make it a bit nicer, 
par(mar=c(5,11,1,1))
barplot(feat$Count, names.arg = feat$Feature, horiz=TRUE,
        las=1, main="Some title")
      
# section 3 color
counts<- read.table("bimm143_05_rstats/male_female_counts.txt", 
                    sep = "\t",
                    header= TRUE)
counts<- read.delim("bimm143_05_rstats/male_female_counts.txt")

barplot(counts$Count, names.arg = counts$Sample, las=2,
        col=rainbow(10) )
barplot(counts$Count, names.arg = counts$Sample, las=2,
        col= c(1,2))

barplot(counts$Count, names.arg = counts$Sample, las=2,
        col=rainbow(10))
