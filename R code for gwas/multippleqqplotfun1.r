 # umesh rosyara 8/23/2012
 qqplotfun <- function (x1, x2) {
   if (!is.numeric(x1)) stop("D'oh! X1 P value vector is not numeric.")
   x1 <- x1[!is.na(x1) & x1<1 & x1>0]
   if (!is.numeric(x2)) stop("D'oh! X2 P value vector is not numeric.")
   x2 <- x2[!is.na(x2) & x2<1 & x2>0]
 x <- c(-log10(x1),-log10(x2))
  qx <- qqnorm(x)
  qx$gr <- c(rep(1, length(x1)), rep(2, length (x2)))
 df1 <- data.frame ( x = qx$x, y= qx$y, gr = qx$gr)
  qx1 <- df1[df1$gr==1,]
 qx2 <- df1[df1$gr==2,] 
 plot(qx1$x,qx1$y, xlab = "Theoritical Quantiles", ylab = "log10 (p-value)", ylim = c(0, max(x)), xlim = c(-3, 3))
 points (qx2$x,qx2$y, col = "red", pch = 18)
 qqline(x, col = 3)

 }
 # example 
x1 <- rnorm (100, 0.5, 1)
x2 <- rnorm (100, 0.3, 1)
qqplotfun (x1, x2)


 
 

 