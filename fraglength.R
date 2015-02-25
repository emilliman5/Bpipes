#!/biotools/R/bin/Rscript

argv<-commandArgs(TRUE)

value=scan(argv[1])
values=subset(value, value >=0)
png(argv[2], height=800, width=1600, units="px")
hist(values, breaks=800, xlim = c(25, 1000), main=argv[1], xlab="IntraMate Fragment Length (bp)")
legend("topright", cex=1.5, paste(names(summary(values)), format(summary(values), digits=2), collapse="\n"), bty="n")
dev.off()