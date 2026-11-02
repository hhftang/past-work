# In Bash
# plink --bfile cas_imp_3 --assoc --pheno pheno/SCORESUM_PC1_std.phe --out PC1
# awk 'BEGIN{OFS="\t";} {print $2, $1, $3, $9;}' PC1.qassoc > PC1.txt
# plink --bfile cas_imp_3 --linear --pheno pheno/SCORESUM_PC1_std.phe --out PC1
# awk 'BEGIN{OFS="\t";} {print $2, $1, $3, $9;}' PC1.assoc.linear > PC1.txt

# In R
library(ggplot2)
library(qqman)
# devtools::install_github("cran/GenABEL")
# library(GenABEL)

args = commandArgs(T)
print(args)

task = args[1]
print(paste0("Performing ",task))

print(paste0("Completed loading table into memory for ",task))

master.df <- read.table(paste0(task, ".txt"), header=T)

png(filename=paste0("output/manplot_", task, ".png"), width=1440, height=720)
manhattan(master.df, ylim=c(0,max(max(-log10(master.df$P)+1),10)),
	main = paste0("Manhattan Plot - ",task), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()
print(paste0("Completed manplot for ",task))

png(filename=paste0("output/qqplot_",task,".png"), width=720,height=720)
qq(master.df$P, ylim=c(0,max(max(-log10(master.df$P)),10)),
	main = paste0("QQ plot -",task), cex=1.2, cex.axis=1.1, cex.lab=1.1)
dev.off()

print(paste0("Completed qqplot for ",task))

sink(paste0("output/lambda_",task,".txt"))
pvalue <- master.df$P
chisq <- qchisq(1 - pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq) ) ) )
lci <- lambda - (SE.median / qchisq(0.5, 1))
uci <- lambda + (SE.median / qchisq(0.5, 1))
print(paste0("lambda = ", lambda))
print(paste0("lci = ", lci))
print(paste0("uci = ", uci))
sink()
print(paste0("Completed lambda for ",task))

# In bash, for locuszoom
cat "$TARGET" | awk -v chr="$CHR" '{if ($2 == chr) {print $2,$3,$4,$1}}' |\
	cut -f1 -d":" | sort -k 1,2 -g  > locuszoom_output.txt
echo 'Chr BP P.value MarkerName'| cat - locuszoom_output.txt \
	> temp && mv temp locuszoom_output.txt