load("output/2018-02-06_process_phyloseq_1.Rdata")

rm(list=setdiff(ls(), c("cas.biom.tr.cas.qc", "cas.biom.tr.cas.qc.css", 
												"cas.biom.method1.css")))

cas.biom.qc <- cas.biom.tr.cas.qc
cas.biom.qc.css <- cas.biom.tr.cas.qc.css
rm(cas.biom.tr.cas.qc)
rm(cas.biom.tr.cas.qc.css)
save.image("data/cas_phyloseq_files_qc.Rdata")
