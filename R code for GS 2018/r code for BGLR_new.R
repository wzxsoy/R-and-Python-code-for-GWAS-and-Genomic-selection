
library(BGLR)
library(rrBLUP)

  Markers <- as.matrix (read.table(file="Panel1_cnn.txt"), header=F)
  Marker50 <- Markers[,sample (1:5403,5000)]
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
  
  train= as.matrix(sample(1:703, 520))
  test<-sample (setdiff(1:703,train),120)
  Pheno_train=Pheno[train,]
  m_train=Markers_impute[train,]
  Pheno_valid=Pheno[test,]
  m_valid=Markers_impute[test,]
  yield=(Pheno_train[,1])
  
  
  ETA<-list(
    list(X=m_train, model='BayesA')
  )
  
  fm<-BGLR(y=yield,ETA=ETA, nIter=12000, burnIn=2000)
  YLD = fm$ETA[[1]]$b
  e = as.matrix(YLD)
  
  pred_yield_valid =  m_valid %*% e
  pred_yield=(pred_yield_valid[,1])
  pred_yield
  yield_valid = Pheno_valid[,1]
  accuracy <-cor(pred_yield_valid, yield_valid, use="complete" )
 accuracy
 
 plot(pred_yield_valid, yield_valid)
 
 