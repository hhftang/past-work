# 2017-02-10_virus_baseline_timeframe_check

# ----

# Use virus_y1.df and virus_y23.df (derived from virus NPA dataset 
# 2015-12-18_NPA_Virus_y1.csv (originally ICHR_Post_nasel_aspirates_Yr1_final.xlsx)
# and 2015-12-18_NPA_Virus_y23.csv (originally CAS_PNA_Virus_parsed_Cryobank.xlsx)

virus_y1_control.df <- with(virus_y1.df[which(virus_y1.df$Category == "control"),],
														data.frame(SUBJ_ID = subjid, age_control = age, LogNo_control = LogNoS))
virus_y1_infect.df <- with(virus_y1.df[which(virus_y1.df$Category == "Infectious"),],
													 data.frame(SUBJ_ID = subjid, age_infect = age, LogNo_infect = LogNoS))

virus_y1_control_infect.df <- merge(virus_y1_control.df, virus_y1_infect.df, 
																		by = "SUBJ_ID")
virus_y1_control_infect.df$age_diff <- virus_y1_control_infect.df$age_control - 
	virus_y1_control_infect.df$age_infect

virus_y1_control_infect.df$bad <- (virus_y1_control_infect.df$age_diff>0 &
																	 	virus_y1_control_infect.df$age_diff<=28)
virus_y1_LogNo.df <- data.frame(table(virus_y1_control_infect.df[
	which(virus_y1_control_infect.df$bad == T),]$LogNo_control))
sum(table(virus_y1_LogNo.df$Freq))

# 455 baselines
# 978 infection samples
# 1652 combinations of infections vs. baseline, matched by SUBJ_ID
# Of these 1652, 52 combinations involve baseline within 4 weeks post-infection. 
# 52 combinations map to 52 baselines.

virus_y2_control.df <- with(virus_y2.df[which(virus_y2.df$Category == "control"),],
														data.frame(SUBJ_ID = subjid, age_control = age, LogNo_control = LogNoG))
virus_y2_infect.df <- with(virus_y2.df[which(virus_y2.df$Category == "Infectious"),],
													 data.frame(SUBJ_ID = subjid, age_infect = age, LogNo_infect = LogNoG))

virus_y2_control_infect.df <- merge(virus_y2_control.df, virus_y2_infect.df, 
																		by = "SUBJ_ID")
virus_y2_control_infect.df$age_diff <- virus_y2_control_infect.df$age_control - 
	virus_y2_control_infect.df$age_infect

virus_y2_control_infect.df$bad <- (virus_y2_control_infect.df$age_diff>0 &
																	 	virus_y2_control_infect.df$age_diff<=28)
virus_y2_LogNo.df <- data.frame(table(virus_y2_control_infect.df[
	which(virus_y2_control_infect.df$bad == T),]$LogNo_control))
sum(table(virus_y2_LogNo.df$Freq))

# 312 baselines
# 913 infection samples
# 1093 combinations of infections vs. baseline, matched by SUBJ_ID
# Of these 1093, 35 combinations involve baseline within 4 weeks post-infection. 
# 35 combinations map to 34 baselines.

virus_y3_control.df <- with(virus_y3.df[which(virus_y3.df$Category == "control"),],
														data.frame(SUBJ_ID = subjid, age_control = age, LogNo_control = LogNoG))
virus_y3_infect.df <- with(virus_y3.df[which(virus_y3.df$Category == "Infectious"),],
													 data.frame(SUBJ_ID = subjid, age_infect = age, LogNo_infect = LogNoG))

virus_y3_control_infect.df <- merge(virus_y3_control.df, virus_y3_infect.df, 
																		by = "SUBJ_ID")
virus_y3_control_infect.df$age_diff <- virus_y3_control_infect.df$age_control - 
	virus_y3_control_infect.df$age_infect

virus_y3_control_infect.df$bad <- (virus_y3_control_infect.df$age_diff>0 &
																	 	virus_y3_control_infect.df$age_diff<=28)
virus_y3_LogNo.df <- data.frame(table(virus_y3_control_infect.df[
	which(virus_y3_control_infect.df$bad == T),]$LogNo_control))
sum(table(virus_y3_LogNo.df$Freq))

# 410 baselines
# 781 infection samples
# 1393 combinations of infections vs. baseline, matched by SUBJ_ID
# Of these 1393, 76 combinations involve baseline within 4 weeks post-infection. 
# 76 combinations map to 76 baselines.

# ----

# Use ari.df (derived from clinical dataset CAS_PNA_cyrobanked_171013.txt)

ari_control.df <- with(ari.df[which(ari.df$Category == "control"),],
														data.frame(SUBJ_ID = subjid, age_control = age, LogNo_control = LogNo))
ari_infect.df <- with(ari.df[which(ari.df$Category == "Infectious"),],
													 data.frame(SUBJ_ID = subjid, age_infect = age, LogNo_infect = LogNo))

ari_control_infect.df <- merge(ari_control.df, ari_infect.df, 
																		by = "SUBJ_ID")
ari_control_infect.df$age_diff <- ari_control_infect.df$age_control - 
	ari_control_infect.df$age_infect

ari_control_infect.df$bad <- (ari_control_infect.df$age_diff>0 &
																	 	ari_control_infect.df$age_diff<=28)
ari_LogNo.df <- data.frame(table(ari_control_infect.df[
	which(ari_control_infect.df$bad == T),]$LogNo_control))
sum(table(ari_LogNo.df$Freq))

# 1943 baselines
# 3657 infection timepoints (not necessarily samples)
# 29654 combinations of infections vs. baseline, matched by SUBJ_ID
# Of these 29654, 309 combinations involve baseline within 4 weeks post-infection. 
# 309 combinations map to 306 baselines.
