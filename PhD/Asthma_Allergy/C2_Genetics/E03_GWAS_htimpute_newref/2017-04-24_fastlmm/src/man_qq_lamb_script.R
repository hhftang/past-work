library(ggplot2)
library(qqman)

args = commandArgs(T)
print(args)

task = args[1]
print(paste0("Performing ",task))

master = read.table(paste0("output/",task,"/master_",task,"_sort.txt"),header=F)
colnames(master) = c("SNP","CHR","BP","P")

print(paste0("Completed loading table into memory for ",task))

png(filename=paste0("output/",task,"/manplot_",task,".png"), width=1440, height=720)
manhattan(master, ylim=c(0,max(max(-log10(master$P)+1),10)),
	main = paste0("Manhattan Plot - ",task), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()
print(paste0("Completed manplot for ",task))

png(filename=paste0("output/",task,"/qqplot_",task,".png"), width=720,height=720)
qq(master$P, ylim=c(0,max(max(-log10(master$P)),10)),
	main = paste0("QQ plot -",task), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()
print(paste0("Completed qqplot for ",task))

library(GenABEL)

sink(paste0("output/",task,"/lambda_",task,".txt"))
estlambda(master$P,plot=F)
sink()
print(paste0("Completed lambda for ",task))

q()