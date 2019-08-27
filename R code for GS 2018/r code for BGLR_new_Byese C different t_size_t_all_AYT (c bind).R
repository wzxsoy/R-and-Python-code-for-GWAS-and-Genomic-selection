library(BGLR)
library(rrBLUP)

Markers <- as.matrix (read.table(file="AYT_hapn.txt"), header=F)
Marker50 <- Markers
head(Markers)
Pheno <-as.matrix (read.table(file ="ayt_t.txt", header=TRUE))
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
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 50))
  test<-sample (setdiff(1:692,train),100)
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
write.csv (accuracy, file= "validation_50t1.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 100))
  test<-sample (setdiff(1:692,train),100)
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
  train= as.matrix(sample(1:692, 200))
  test<-sample (setdiff(1:692,train),100)
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
  train= as.matrix(sample(1:692, 300))
  test<-sample (setdiff(1:692,train),100)
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
  train= as.matrix(sample(1:692, 400))
  test<-sample (setdiff(1:692,train),100)
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
  train= as.matrix(sample(1:692, 500))
  test<-sample (setdiff(1:692,train),100)
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
  train= as.matrix(sample(1:692, 600))
  test<-sample (setdiff(1:692,train),90)
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

##Trait2



traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 50))
  test<-sample (setdiff(1:692,train),100)
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
write.csv (accuracy, file= "validation_50t2.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 100))
  test<-sample (setdiff(1:692,train),100)
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
write.csv (accuracy, file= "validation_100t2.csv")
mean(accuracy)
accuracy

#200t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 200))
  test<-sample (setdiff(1:692,train),100)
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
write.csv (accuracy, file= "validation_200t2.csv")
mean(accuracy)
accuracy

#300t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 300))
  test<-sample (setdiff(1:692,train),100)
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
write.csv (accuracy, file= "validation_300t2.csv")
mean(accuracy)
accuracy

#400t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 400))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_400t2.csv")
mean(accuracy)
accuracy


#500

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 500))
  test<-sample (setdiff(1:692,train),100)
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
write.csv (accuracy, file= "validation_500t2.csv")
mean(accuracy)
accuracy

#600

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 600))
  test<-sample (setdiff(1:692,train),90)
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
write.csv (accuracy, file= "validation_600t2.csv")
mean(accuracy)
accuracy


#Trait3 of indust.

traits=1
cycles=10
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 50))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_50t3.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 100))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_100t3.csv")
mean(accuracy)
accuracy

#200t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 200))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_200t3.csv")
mean(accuracy)
accuracy

#300t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 300))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_300t3.csv")
mean(accuracy)
accuracy

#400t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 400))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_400t3.csv")
mean(accuracy)
accuracy


#500

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 500))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_500t3.csv")
mean(accuracy)
accuracy

#600

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 600))
  test<-sample (setdiff(1:692,train),90)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_600t3.csv")
mean(accuracy)
accuracy


#trait 4 of indust.

traits=1
cycles=10
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 50))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_50t4.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 100))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_100t4.csv")
mean(accuracy)
accuracy

#200t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 200))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_200t4.csv")
mean(accuracy)
accuracy

#300t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 300))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_300t4.csv")
mean(accuracy)
accuracy

#400t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 400))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_400t4.csv")
mean(accuracy)
accuracy


#500

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 500))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_500t4.csv")
mean(accuracy)
accuracy

#600

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 600))
  test<-sample (setdiff(1:692,train),90)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_600t4.csv")
mean(accuracy)
accuracy


###Trait 5 of indus.


traits=1
cycles=10
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 50))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_50t5.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 100))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_100t5.csv")
mean(accuracy)
accuracy

#200t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 200))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_200t5.csv")
mean(accuracy)
accuracy

#300t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 300))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_300t5.csv")
mean(accuracy)
accuracy

#400t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 400))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_400t5.csv")
mean(accuracy)
accuracy


#500

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 500))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_500t5.csv")
mean(accuracy)
accuracy

#600

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 600))
  test<-sample (setdiff(1:692,train),90)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_600t5.csv")
mean(accuracy)
accuracy


###trait 6 of indus.



traits=1
cycles=10
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 50))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_50t6.csv")
mean(accuracy)
accuracy


#t_100 

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 100))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_100t6.csv")
mean(accuracy)
accuracy

#200t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 200))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_200t6.csv")
mean(accuracy)
accuracy

#300t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 300))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_300t6.csv")
mean(accuracy)
accuracy

#400t

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 400))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_400t6.csv")
mean(accuracy)
accuracy


#500

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 500))
  test<-sample (setdiff(1:692,train),100)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_500t6.csv")
mean(accuracy)
accuracy

#600

traits=1
cycles=30
accuracy = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles)
{
  train= as.matrix(sample(1:692, 600))
  test<-sample (setdiff(1:692,train),90)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,6])
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_600t6.csv")
mean(accuracy)
accuracy


a <- read.csv("t_size_50t1.csv")
b <- read.csv("t_size_100t1.csv")
c <- read.csv("t_size_200t1.csv")
d <- read.csv("t_size_300t1.csv")
e <- read.csv("t_size_400t1.csv")
f <- read.csv("t_size_500t1.csv")
g <- read.csv("t_size_600t1.csv")
h <- read.csv("t_size_670t1.csv")


accuracy_BA <- cbind(a[,2],b[,2],c[,2],d[,2],e[,2],f[,2],g[,2],h[,2])
write.csv (accuracy_BA, file= "accuracy__RKHS__t_all.csv")
write.csv (fm$fit, file= "Fit_RKHS_all.csv")
