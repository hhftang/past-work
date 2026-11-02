# Simple plots

df <- with(master_GRS.df, 
					 data.frame(SUBJ_ID, SCORESUM_allergic_tert, 
					 					 cras_3, cras_4, cras_5, cras_10)) %>%
	melt(id.vars = c("SUBJ_ID", "SCORESUM_allergic_tert"),
			 variable.name = "timepoint",
			 value.name = "cras") %>%
	mutate(tertile = SCORESUM_allergic_tert,
				 timepoint = as.numeric(as.character(
		factor(timepoint, levels = c("cras_3", "cras_4", "cras_5", "cras_10"),
					 labels = c(3, 4, 5, 10))))) %>%
	mutate(tertile = factor(tertile, levels=c("Low", "Med", "High")))
df2 <- data.frame(with(df, table(timepoint, cras, tertile))) %>%
	group_by(timepoint, tertile) %>%
	mutate(sum = sum(Freq),
				 prop = Freq/sum, 
				 se = 1.96*(prop*(1- prop)/sum)^0.5,
				 lower = prop - se, 
				 upper = prop + se) %>%
	filter(cras == 1)

pdf("output/2019-01-09_GRS_simple_asthma.pdf", width=4, height=2.5)
df2 %>% 
	ungroup %>%
	ggplot(aes(x = timepoint, y = prop)) + 
	geom_line(aes(group = tertile, color = tertile)) +
	geom_point(aes(group = tertile, color = tertile, shape = tertile)) + 
	scale_color_manual(values = divPalette.alt2(3)) + 
	scale_fill_manual(values = divPalette.alt2(3)) + 
	labs(x = "Timepoint or age", y = "Proportion with asthma", 
			 color = "GRS tertile",  shape = "GRS tertile")
dev.off()

df <- with(master_GRS.df, 
					 data.frame(SUBJ_ID, SCORESUM_allergic_std, 
					 					 group))

pdf("output/2019-01-09_GRS_simple_npEM.pdf", width=4, height=2.5)
df %>% 
	filter(!is.na(group)) %>%
	ggplot(aes(x = factor(group, labels=c("CAS1", "CAS2", "CAS3")),
												y = SCORESUM_allergic_std, color = group, 
						 fill = group)) + 
	geom_violin(alpha=0.1) + 
	geom_jitter(width = 0.1, height = 0) + 
	scale_color_manual(values = c("#8da0cb","#a6d854","#fc8d62")) + 
	scale_fill_manual(values = c("#8da0cb","#a6d854","#fc8d62")) +
	labs(y = "GRS", x = "npEM cluster") +
	theme(legend.position="none")
dev.off()


df <- with(master_GRS.df, 
					 data.frame(SUBJ_ID, SCORESUM_allergic_std, 
					 					 illness.mpg.healthy.age2.bin))

pdf("output/2019-01-09_GRS_simple_bad_mpg.pdf", width=4, height=2.5)
df %>% 
	filter(!is.na(illness.mpg.healthy.age2.bin)) %>%
	ggplot(aes(x = factor(illness.mpg.healthy.age2.bin, 
												levels=c("Low", "High")),
						 y = SCORESUM_allergic_std, color = illness.mpg.healthy.age2.bin, 
						 fill = illness.mpg.healthy.age2.bin)) + 
	geom_violin(alpha=0.1) + 
	geom_jitter(width = 0.1, height = 0) + 
	scale_color_manual(values = c("#D73027","#4575B4")) + 
	scale_fill_manual(values = c("#D73027","#4575B4")) +
	labs(y = "GRS", x = "High/Low prop. of illness-associated MPG") + 
	theme(legend.position="none")
dev.off()