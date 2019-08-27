    qq = function(pvector, ...) {  
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")	
pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]  
  o = -log10(sort(pvector,decreasing=F))	
#e = -log10( 1:length(o)/length(o) )	
e = -log10( ppoints(length(pvector) ))
	plot(e,o,pch=20,cex=0.8, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), col="Blue",cex.axis=1.2)
	abline(0,1,col="red",lwd=2)}


 mah= function(dataframe, colors=c("black","#666666","#CC6600"), pch=20, ymax="max", cex.x.axis=1, limitchromosomes=1:20, suggestiveline=-log10(5.46E-05), genomewideline=-log10(5.46E-05), annotate=NULL, ...) {

    d=dataframe
    if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
    
    if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
    d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
    if (ymax=="max") ymax<-ceiling(max(d$logp))
    if (ymax<8) ymax<-8
    
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
        d$pos=d$BP
        ticks=floor(length(d$pos))/2+1
    } else {
        for (i in unique(d$CHR)) {
        	if (i==1) {
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
    		} else {
    			lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
    			d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
    		}
    		ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    	}
    }
    
    if (numchroms==1) {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), pch=20))
    }	else {
        with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", pch=20))
        axis(1, at=ticks, lab=unique(d$CHR), cex.axis=cex.x.axis)
        icol=1
        for (i in unique(d$CHR)) {
            with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], pch=20))
            icol=icol+1
    	}
    }
    
    if (!is.null(annotate)) {
        d.annotate=d[which(d$SNP %in% annotate), ]
        with(d.annotate, points(pos, logp, col="green3", pch=20)) 
    }
    
    if (suggestiveline) abline(h=suggestiveline, col="blue")
    if (genomewideline) abline(h=genomewideline, col="red")
}



dig<-read.table("blupmlmayt.txt")
  digp<-data.frame(dig)
  windows.options(width=5, height=5)
  qq(digp$P)

windows.options(width=100, height=30)
mah(digp, colors=c("#66CD00","#CC6600"), pch=20,cex=1.3,cex.axis=1.5)