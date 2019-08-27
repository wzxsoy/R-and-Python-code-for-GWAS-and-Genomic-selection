qualdat = read.csv("maf.csv", header=T)
attach(qualdat)
plot(density(LD),col="blue")
lines(density(im),col="red")

library("ggplot2")
my.df <-read.csv("maf.csv", header=T)
ggplot(my.df, aes(x=dat, color="red",fill="white")) + geom_bar(position="dodge")

ggplot(my.df, aes(x=rating, fill=cond,width=0.5)) + geom_histogram(binwidth=0.1, position = position_dodge(width=0.09))+scale_fill_manual(values=c("blue", "red"))+
theme(panel.background = element_blank()) +scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5))+
theme(axis.text.x = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=.5,face="plain"))

   

ggplot(my.df, aes(x=rating, fill=cond,width=.85)) + geom_histogram(binwidth=0.05, position="dodge")+ geom_bar(width = 1.8, position = position_dodge(width = 0.9))+
theme(panel.background = element_blank())
