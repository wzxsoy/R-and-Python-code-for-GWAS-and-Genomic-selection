### First step, run the following two functions:GAPIT.Numericalization and 

GAPIT.Numericalization <-
  function(x,bit=2,effect="Add",impute="None", Create.indicator = FALSE, Major.allele.zero = FALSE, byRow=TRUE){
    #Object: To convert character SNP genotpe to numerical
    #Output: Coresponding numerical value
    #Authors: Feng Tian and Zhiwu Zhang
    # Last update: May 30, 2011 
    ##############################################################################################
    if(bit==1)  {
      x[x=="X"]="N"
      x[x=="-"]="N"
      x[x=="+"]="N"
      x[x=="/"]="N"
      x[x=="K"]="Z" #K (for GT genotype)is replaced by Z to ensure heterozygose has the largest value
    }
    
    if(bit==2)  {
      x[x=="XX"]="N"
      x[x=="--"]="N"
      x[x=="++"]="N"
      x[x=="//"]="N"
      x[x=="NN"]="N"
    }
    
    n=length(x)
    lev=levels(as.factor(x))
    lev=setdiff(lev,"N")
    #print(lev)
    len=length(lev)
    #print(lev)
    
    
    
    #Genotype counts
    count=1:len
    for(i in 1:len){
      count[i]=length(x[(x==lev[i])])
    }
    
    
    
    if(Major.allele.zero){
      if(len>1 & len<=3){
        #One bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the second position
        if(bit==1){ 
          count.temp = cbind(count, seq(1:len))
          if(len==3) count.temp = count.temp[-3,]
          count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
          if(len==3)order =  c(count.temp[,2],3)else order = count.temp[,2]
        }
        
        #Two bit: Make sure that the SNP with the major allele is on the top, and the SNP with the minor allele is on the third position
        if(bit==2){ 
          count.temp = cbind(count, seq(1:len))
          if(len==3) count.temp = count.temp[-2,]
          count.temp <- count.temp[order(count.temp[,1], decreasing = TRUE),]
          if(len==3) order =  c(count.temp[1,2],2,count.temp[2,2])else order = count.temp[,2]
        }
        
        count = count[order]
        lev = lev[order]
        
      }   #End  if(len<=1 | len> 3)
    } #End  if(Major.allele.zero)
    
    
    
    #make two  bit order genotype as AA,AT and TT, one bit as A(AA),T(TT) and X(AT)
    if(bit==1 & len==3){
      temp=count[2]
      count[2]=count[3]
      count[3]=temp
    }
    position=order(count)
    
    
    #1status other than 2 or 3
    if(len<=1 | len> 3)x=-1
    
    #2 status
    if(len==2)x=ifelse(x=="N",NA,ifelse(x==lev[1],-1,1))
    
    #3 status
    if(bit==1){
      if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],-1,ifelse(x==lev[3],0,1)))
    }else{
      if(len==3)x=ifelse(x=="N",NA,ifelse(x==lev[1],-1,ifelse(x==lev[3],1,0)))
    }
    
    #print(paste(lev,len,sep=" "))
    #print(position)
    
    #missing data imputation
    if(impute=="Middle") {x[is.na(x)]=1 }
    
    if(len==3){
      if(impute=="Minor")  {x[is.na(x)]=position[1]  -1}
      if(impute=="Major")  {x[is.na(x)]=position[len]-1}
      
    }else{
      if(impute=="Minor")  {x[is.na(x)]=2*(position[1]  -1)}
      if(impute=="Major")  {x[is.na(x)]=2*(position[len]-1)}
    }
    
    #alternative genetic models
    if(effect=="Dom") x=ifelse(x==1,1,0)
    if(effect=="Left") x[x==1]=0
    if(effect=="Right") x[x==1]=2
    
    if(byRow) {
      result=matrix(x,n,1)
    }else{
      result=matrix(x,1,n)  
    }
    
    return(result)
  }
#end of GAPIT.Numericalization function



# Beginning of GAPIT.HapMap function


GAPIT.HapMap <-
  function(G,SNP.effect="Add",SNP.impute="Middle",heading=TRUE, Create.indicator = FALSE, Major.allele.zero = FALSE){
    #Object: To convert character SNP genotpe to numerical
    #Output: Coresponding numerical value
    #Authors: Feng Tian and Zhiwu Zhang
    # Last update: May 30, 2011 
    ##############################################################################################
    
    print(paste("Converting HapMap format to numerical under model of ", SNP.impute,sep=""))
    #gc()
    #GAPIT.Memory.Object(name.of.trait="HapMap.Start")
    
    #GT=data.frame(G[1,-(1:11)])
    if(heading){
      GT= t(G[1,-(1:11)])
      GI= G[-1,c(1,3,4)]
    }else{
      GT=NULL
      GI= G[,c(1,3,4)]
    }
    
    
    #Set column names
    if(heading)colnames(GT)="taxa"
    colnames(GI)=c("SNP","Chromosome","Position")
    
    #Initial GD
    GD=NULL
    bit=nchar(as.character(G[2,12])) #to determine number of bits of genotype
    #print(paste("Number of bits for genotype: ", bit))
    
    print("Perform numericalization")
    
    if(heading){
      if(!Create.indicator) GD= apply(G[-1,-(1:11)],1,function(one) GAPIT.Numericalization(one,bit=bit,effect=SNP.effect,impute=SNP.impute, Major.allele.zero=Major.allele.zero))
      if(Create.indicator) GD= t(G[-1,-(1:11)])
    }else{
      if(!Create.indicator) GD= apply(G[  ,-(1:11)],1,function(one) GAPIT.Numericalization(one,bit=bit,effect=SNP.effect,impute=SNP.impute, Major.allele.zero=Major.allele.zero))
      if(Create.indicator) GD= t(G[ ,-(1:11)])
    }
    
    #set GT and GI to NULL in case of null GD
    if(is.null(GD)){
      GT=NULL
      GI=NULL
    }
    
    print("The dimension of GD is:")
    print(dim(GD))
    write.csv(GD, file ="Numeric data.csv")
    
    
    if(!Create.indicator) {print(paste("Succesfuly finished converting HapMap which has bits of ", bit,sep="")) }
    return(list(GT=GT,GD=GD,GI=GI))
  }#end of GAPIT.HapMap function

###Second step, run the following code to read the hapmap data and convert it to numeric data (-1,0,1) format.

myG <- read.table("The name of your hapmap.txt", head = FALSE)

GAPIT.HapMap (myG)





