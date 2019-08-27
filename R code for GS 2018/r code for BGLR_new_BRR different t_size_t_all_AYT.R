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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1]<-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=100)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
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
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy [r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
}

accuracy
write.csv (accuracy, file= "validation_600t6.csv")
mean(accuracy)
accuracy
library(BGLR)
library(rrBLUP)

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50t1.csv")
mean(accuracy)


#100




traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1001t.csv")
mean(accuracy)

#200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2001t.csv")
mean(accuracy)

#400

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4001t.csv")
mean(accuracy)

#600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_6001t.csv")
mean(accuracy)

#800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_8001t.csv")
mean(accuracy)

#1000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_10001t.csv")
mean(accuracy)



#1200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_12001t.csv")
mean(accuracy)

#1400



traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_14001t.csv")
mean(accuracy)


#1600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_16001t.csv")
mean(accuracy)


#1800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_18001t.csv")
mean(accuracy)


#2000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_20001t.csv")
mean(accuracy)


#3000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_30001t.csv")
mean(accuracy)

#4000
traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,1])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_40001t.csv")
mean(accuracy)

###TRAIT 2


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50t2.csv")
mean(accuracy)


#100




traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1002t.csv")
mean(accuracy)

#200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2002t.csv")
mean(accuracy)

#400

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4002t.csv")
mean(accuracy)

#600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_6002t.csv")
mean(accuracy)

#800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_8002t.csv")
mean(accuracy)

#1000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_10002t.csv")
mean(accuracy)



#1200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_12002t.csv")
mean(accuracy)

#1400



traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_14002t.csv")
mean(accuracy)


#1600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_16002t.csv")
mean(accuracy)


#1800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_18002t.csv")
mean(accuracy)


#2000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_20002t.csv")
mean(accuracy)


#3000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_30002t.csv")
mean(accuracy)

#4000
traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,2])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,2]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_40002t.csv")
mean(accuracy)

###Trait 3


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50t3.csv")
mean(accuracy)


#100




traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1003t.csv")
mean(accuracy)

#200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2003t.csv")
mean(accuracy)

#400

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4003t.csv")
mean(accuracy)

#600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_6003t.csv")
mean(accuracy)

#800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_8003t.csv")
mean(accuracy)

#1000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_10003t.csv")
mean(accuracy)



#1200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_12003t.csv")
mean(accuracy)

#1400



traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_14003t.csv")
mean(accuracy)


#1600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_16003t.csv")
mean(accuracy)


#1800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_18003t.csv")
mean(accuracy)


#2000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_20003t.csv")
mean(accuracy)


#3000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_30003t.csv")
mean(accuracy)

#4000
traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,3])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,3]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_40003t.csv")
mean(accuracy)

###Trait 4

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50t4.csv")
mean(accuracy)


#100




traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_100t4.csv")
mean(accuracy)

#200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_200t4.csv")
mean(accuracy)

#400

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_400t4.csv")
mean(accuracy)

#600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_600t4.csv")
mean(accuracy)

#800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_800t4.csv")
mean(accuracy)

#1000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1000t4.csv")
mean(accuracy)



#1200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1200t4.csv")
mean(accuracy)

#1400



traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1400t4.csv")
mean(accuracy)


#1600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1600t4.csv")
mean(accuracy)


#1800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1800t4.csv")
mean(accuracy)


#2000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2000t4.csv")
mean(accuracy)


#3000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_3000t4.csv")
mean(accuracy)

#4000
traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,4])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,4]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4000t4.csv")
mean(accuracy)

#####Trait5

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50t5.csv")
mean(accuracy)


#100




traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_100t5.csv")
mean(accuracy)

#200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_200t5.csv")
mean(accuracy)

#400

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_400t5.csv")
mean(accuracy)

#600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_600t5.csv")
mean(accuracy)

#800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_800t5.csv")
mean(accuracy)

#1000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1000t5.csv")
mean(accuracy)



#1200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1200t5.csv")
mean(accuracy)

#1400



traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1400t5.csv")
mean(accuracy)


#1600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1600t5.csv")
mean(accuracy)


#1800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1800t5.csv")
mean(accuracy)


#2000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2000t5.csv")
mean(accuracy)


#3000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_3000t5.csv")
mean(accuracy)

#4000
traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  yield=(Pheno_train[,5])
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,5]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4000t5.csv")
mean(accuracy)

####Trait 6

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_50t6.csv")
mean(accuracy)


#100




traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_100t6.csv")
mean(accuracy)

#200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_200t6.csv")
mean(accuracy)

#400

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_400t6.csv")
mean(accuracy)

#600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_600t6.csv")
mean(accuracy)

#800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_800t6.csv")
mean(accuracy)

#1000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1000t6.csv")
mean(accuracy)



#1200


traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1200t6.csv")
mean(accuracy)

#1400



traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1400t6.csv")
mean(accuracy)


#1600

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1600t6.csv")
mean(accuracy)


#1800

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_1800t6.csv")
mean(accuracy)


#2000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_2000t6.csv")
mean(accuracy)


#3000

traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_3000t6.csv")
mean(accuracy)

#4000
traits=1
cycles=30
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
  
  impute=A.mat(Marker50,max.missing=0.5,impute.method="mean",return.imputed=T)
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
  
  ETA<-list(
    list(X=m_train, model='BRR')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=1000, burnIn=500)
  YLD = fm$ETA[[1]]$b
  
  e = as.matrix(YLD)
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,6]
  accuracy[r,1] <-cor(pred_yield_valid, yield_valid, use="complete" )
  
}
accuracy
write.csv (accuracy, file= "SNP_4000t6.csv")
mean(accuracy)


###End.

