load("output/2018-01-31_process_phyloseq_1.Rdata")

rm(list=setdiff(ls(), c("coast.biom.tr.sum.coast.qc2", 
												"coast.biom.tr.sum.coast.qc2.css", 
												"coast.biom.method1.css")))

coast.biom.qc <- coast.biom.tr.sum.coast.qc2
coast.biom.qc.css <- coast.biom.tr.sum.coast.qc2.css
rm(coast.biom.tr.sum.coast.qc2)
rm(coast.biom.tr.sum.coast.qc2.css)

save.image("data/coast_phyloseq_files_qc.Rdata")

