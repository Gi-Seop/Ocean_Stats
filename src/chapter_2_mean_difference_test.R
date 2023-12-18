

## Read input data
DO <- read.table("../data/BK1417_DO.txt", header=FALSE)

DOS <- DO$V1
DOB <- DO$V2

## randomness test {randtests}
library(randtests)
library(trend)

runs.test(DOS)
runs.test(DOB)

bartels.rank.test(DOS)
bartels.rank.test(DOB)

# Normality test

shapiro.test(DOS)
shapiro.test(DOB)

library(nortest)
ad.test(DOS)
ad.test(DOB)

## Variance test + Mean difference test
var.test(DOS, DOB)
t.test(DOS, DOB, var.equal=TRUE)

## Power test of the mean-differerce test

library(pwr)

ES <- (mean(DOS) - mean(DOB))/sd(c(DOS, DOB))

pwr.t.test(n = length(DOS), d = ES, sig.level = 0.05, power = NULL, 
    type = "paired", alternative = "two.sided")
