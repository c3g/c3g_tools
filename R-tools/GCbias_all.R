args=commandArgs(TRUE)
#print(args)
for (j in 1:length(args)) {
    input<-args[j]
    t<-read.table(input)
    c<-cor.test(t$V4,t$V5,method="spearman",exact=FALSE)
    #c<-cor.test(t$V4,t$V5,method="pearson")
    print(input)
    print(c$estimate)
}
