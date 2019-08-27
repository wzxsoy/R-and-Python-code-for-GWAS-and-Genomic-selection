qualdat = read.csv("pi.csv", header=T)

## Check to ensure data imported correctly
str(qualdat)
head(qualdat)
tail(qualdat)

## Attach dataset
attach(qualdat)

## Examine distribution of brix data
hist(DSI, col="gold")
boxplot(yield~Loc, xlab="Location", ylab="Yield", main="Yield by Location", col="pink")

# Rename variables for ease of use
Yield= as.numeric(yield)
LINE = as.factor(Line)
LOC = as.factor(Loc)
YEAR = as.factor(Year)
REP = as.factor(Rep)

## Calculate variance components
# requires lme4 package
library(lme4)

# Linear Model with random effects for variance components
brixvarcomp = lmer(Yield~ (1|LINE) + (1|LOC) +(1|REP%in%LOC)+ (1|LINE:LOC))
# Extract variance components
summary(brixvarcomp)

## BLUPS
# fit the model
brixmodel= lmer(Yield~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))

# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="2011yield-BLUPS-new.csv")




## Examine distribution of brix data
hist(mat, col="gold")
boxplot(mat~Loc, xlab="Location", ylab="Maturity", main="Maturity by Location", col="pink")

# Rename variables for ease of use
MAT= as.numeric(mat)
LINE = as.factor(Line)
LOC = as.factor(Loc)
YEAR = as.factor(Year)
REP = as.factor(Rep)

## Calculate variance components
# requires lme4 package
library(lme4)

# Linear Model with random effects for variance components
brixvarcomp = lmer(MAT~ (1|LINE) + (1|LOC) + (1|LINE:LOC))
# Extract variance components
summary(brixvarcomp)

## BLUPS
# fit the model
brixmodel = lmer(MAT~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))

# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="2007mat-BLUPS.csv")




## Examine distribution of brix data
hist(height, col="gold")
boxplot(height~Loc, xlab="Location", ylab="Height", main="Height by Location", col="pink")

# Rename variables for ease of use
Height= as.numeric(height)
LINE = as.factor(Line)
LOC = as.factor(Loc)
YEAR = as.factor(Year)
REP = as.factor(Rep)

## Calculate variance components
# requires lme4 package
library(lme4)

# Linear Model with random effects for variance components
brixvarcomp = lmer(Height~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))
# Extract variance components
summary(brixvarcomp)

## BLUPS
# fit the model
brixmodel= lmer(Height~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))

# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="2007height-BLUPS.csv")



## Examine distribution of brix data
hist(lodging, col="gold")
boxplot(lodging~Loc, xlab="Location", ylab="Lodging", main="Lodging by Location", col="pink")

# Rename variables for ease of use
Lodging= as.numeric(lodging)
LINE = as.factor(Line)
LOC = as.factor(Loc)
YEAR = as.factor(Year)
REP = as.factor(Rep)

## Calculate variance components
# requires lme4 package
library(lme4)

# Linear Model with random effects for variance components
brixvarcomp = lmer(Lodging~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))
# Extract variance components
summary(brixvarcomp)

## BLUPS
# fit the model
brixmodel = lmer(Lodging~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))

# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="2007Lodging-BLUPS.csv")





## Examine distribution of brix data
hist(pro, col="gold")
boxplot(pro~Loc, xlab="Location", ylab="Protein content(%)", main=" Protein content by Location", col="pink")

# Rename variables for ease of use
Pro= as.numeric(pro)
LINE = as.factor(Line)
LOC = as.factor(Loc)
YEAR = as.factor(Year)
REP = as.factor(Rep)

## Calculate variance components
# requires lme4 package
library(lme4)

# Linear Model with random effects for variance components
brixvarcomp = lmer(Pro~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))
# Extract variance components
summary(brixvarcomp)

## BLUPS
# fit the model
brixmodel = lmer(Pro~ (1|LINE) + (1|LOC)+ (1|LINE:LOC))

# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="2007pro-old-BLUPS.csv")



## Examine distribution of brix data
hist(oil, col="gold")
boxplot(oil~Loc, xlab="Location", ylab="Oil content(%)", main=" Oil content by Location", col="pink")

# Rename variables for ease of use
Oil= as.numeric(oil)
LINE = as.factor(Line)
LOC = as.factor(Loc)
YEAR = as.factor(Year)
REP = as.factor(Rep)

## Calculate variance components
# requires lme4 package
library(lme4)

# Linear Model with random effects for variance components
brixvarcomp = lmer(Oil~ (1|LINE) + (1|LOC)+(1|REP%in%LOC) + (1|LINE:LOC))
# Extract variance components
summary(brixvarcomp)

## BLUPS
# fit the model
brixmodel= lmer(Oil~ (1|LINE) + (1|LOC) +(1|REP%in%LOC)+ (1|LINE:LOC))

# estimate BLUPS
brixblup = ranef(brixmodel)
# look at output structure
str(brixblup)
# extract blup for line
brixlineblup = brixblup$LINE
# see the structure of the blup for each line
str(brixlineblup)
# save the brixlineblup output to a separate .csv file
write.csv(brixlineblup, file="2007oil-BLUPS.csv")






