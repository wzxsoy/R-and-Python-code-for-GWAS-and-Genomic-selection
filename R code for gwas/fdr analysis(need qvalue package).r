 p <- scan("pvalue.txt")
 qobj <- qvalue(p, lambda=0, fdr.level=0.05)
max(qobj$pvalues[qobj$qvalues <= 0.05])


qobj <- qvalue(p)

qplot(qobj)
qwrite(qobj, filename="myresults.txt")
qobj <- qvalue(p, lambda=seq(0.2,0.8,0.2), robust=TRUE)
qobj <- qvalue(p, pi0.meth="bootstrap")