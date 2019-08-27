library(BGLR)
library(rrBLUP)

Markers <- as.matrix (read.table(file="Panel1_cn.txt"), header=F)
Marker50 <- Markers
head(Markers)
Pheno <-as.matrix (read.table(file ="indu all.txt", header=TRUE))
head (Pheno)
# Check the dimensions of the matrix
dim(Markers)
dim(Pheno)

impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
Markers_impute=impute$imputed
head(Markers_impute)
dim(Markers)
dim(Markers_impute)

traits=1
cycles=10
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 50))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_50t.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 100))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_100t1.csv")
mean(accuracy)
accuracy

#200t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 200))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_200t1.csv")
mean(accuracy)
accuracy

#300t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 300))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_300t1.csv")
mean(accuracy)
accuracy

#400t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 400))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_400t1.csv")
mean(accuracy)
accuracy


#500

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 500))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_500t1.csv")
mean(accuracy)
accuracy

#600

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:703, 600))
  test<-sample (setdiff(1:703,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_600t1.csv")
mean(accuracy)
accuracy