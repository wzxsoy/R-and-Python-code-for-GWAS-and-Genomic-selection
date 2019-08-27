

install.packages("rrBLUP")
library(rrBLUP)

options(max.print=5.5E3)
Markers <- as.matrix (read.table(file="Panel1_cnn.txt"), header=F)
head(Markers)
Pheno <-as.matrix (read.table(file ="indu all.txt", header=TRUE))
head (Pheno)
# Check the dimensions of the matrix
dim(Markers)
dim(Pheno)

impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed
head(Markers_impute)
dim(Markers)
dim(Markers_impute)


traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 50))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_50t.csv")
mean(accuracy)





traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 50))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[, 6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[, 1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[, 6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_50t6.csv")
mean(accuracy)



traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 100))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_100t6.csv")
mean(accuracy)



traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 200))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_200t6.csv")
mean(accuracy)



traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 300))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_300t6.csv")
mean(accuracy)



traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 400))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_400t6.csv")
mean(accuracy)




traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 500))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_500t6.csv")
mean(accuracy)



traits=1
cycles=50
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 600))
  test<-sample (setdiff(1:703,train),90)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}
accuracy
write.csv (accuracy, file= "validation_600t6.csv")
mean(accuracy)



traits=1
cycles=10
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="Panel1_cnn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,50)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 120))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50.csv")
mean(accuracy)



traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="Panel1_cnn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,50)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50.csv")
mean(accuracy)


################################
#code for different snp number for ayt number of marker from 50 to 2000

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,50)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50.csv")
mean(accuracy)


#100




traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,100)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_100.csv")
mean(accuracy)

#200


traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,200)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_200.csv")
mean(accuracy)

#400

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,400)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_400.csv")
mean(accuracy)

#600

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,600)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_600.csv")
mean(accuracy)

#800

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,800)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_800.csv")
mean(accuracy)

#1000

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,1000)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1000.csv")
mean(accuracy)



#1200


traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,1200)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1200.csv")
mean(accuracy)

#1400



traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,1400)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1400.csv")
mean(accuracy)


#1600

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,1600)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1600.csv")
mean(accuracy)


#1800

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,1800)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1800.csv")
mean(accuracy)


#2000

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,2000)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2000.csv")
mean(accuracy)


#3000

traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,3000)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_3000.csv")
mean(accuracy)

#4000
traits=1
cycles=20
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  
  Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,4000)]
  head(Markers)
  Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
  head (Pheno)
  # Check the dimensions of the matrix
  dim(Markers)
  dim(Pheno)
  
  impute=A.mat(Markers,max.missing=0.5,impute.method="mean",return.imputed=T)
  Markers_impute=impute$imputed
  head(Markers_impute)
  dim(Markers)
  dim(Markers_impute)
  
  train= as.matrix(sample(1:692, 520))
  test<-sample (setdiff(1:692,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  yield_answer<-mixed.solve(yield, Z=m_train, K=NULL, SE = FALSE, return.Hinv=FALSE)
  YLD = yield_answer$u
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])+yield_answer$beta
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4000.csv")
mean(accuracy)


