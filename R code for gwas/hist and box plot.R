qualdat = read.csv("phist.csv", header=T)

## Check to ensure data imported correctly
str(qualdat)
head(qualdat)
tail(qualdat)

## Attach datase
attach(qualdat)
windows.options(width=5, height=5)
## Examine distribution of brix data
hist(DSI, col=" light blue",breaks=12)



boxplot(yield~Loc, xlab="Location", ylab="Degrees Brix", main="Degrees Brix by Location", col="pink")


