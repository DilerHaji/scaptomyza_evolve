tab <- read.delim("EliSorayaPop3&4 Data - updated min2 20230721ADG.csv", header = T, sep = ",", na.string = "NA")

# Factoring variables 
tab$PhenotypeCageSpecies <- as.factor(tab$PhenotypeCageSpecies)
tab$Group <- as.factor(tab$Group)
tab$PopSpecies <- as.factor(tab$PopSpecies)
tab$PhenotypeCageSpecies <- as.factor(tab$PhenotypeCageSpecies)
tab$PopNumber <- as.factor(tab$PopNumber)

# Defining a alternative replicate label (evolved condition + replicate population)
tab$rep2 <- paste(tab$PopSpecies, tab$PopNumber, sep = "")
tab$rep <- as.factor(substr(tab$Rep, nchar(tab$Rep), nchar(tab$Rep)))

# Counting the number of trays associated with each replicate (should be 6) 
tab$trays <- as.numeric(table(tab$Rep)[match(tab$Rep, names(table(tab$Rep)))])

# Re-leving factors for the evolved condition and the test condition so that the reference is M
tab$PopSpecies <- relevel(tab$PopSpecies, ref = "M")
tab$PhenotypeCageSpecies <- relevel(tab$PhenotypeCageSpecies, ref = "M")

# Total, max, and mean number of flies from 11Feb to 20Feb
tab$numFlies <- as.numeric(apply(tab[, 7:16], 1, function(x){ sum(x[!is.na(x)]) }))
tab$maxFlies <- as.numeric(apply(tab[, 7:16], 1, function(x){ max(x[!is.na(x)]) }))
tab$meanFlies <- as.numeric(apply(tab[, 7:16], 1, function(x){ mean(x[!is.na(x)]) }))

# normalized total fly counts by the number of trays 
tab$normFlies <- tab$numFlies/tab$trays

# Excluding outlier value (this does not have an overall impact on the results, but improves plot visualization) 
tab <- tab[!tab$normFlies == max(tab$normFlies), ]

# metadata 
full <- expand.grid(
  Group = levels(tab$Group),
  rep = levels(tab$rep),
  PopNumber = levels(tab$PopNumber),
  PopSpecies = levels(tab$PopSpecies),
  PhenotypeCageSpecies = levels(tab$PhenotypeCageSpecies),
  trays = 1
)