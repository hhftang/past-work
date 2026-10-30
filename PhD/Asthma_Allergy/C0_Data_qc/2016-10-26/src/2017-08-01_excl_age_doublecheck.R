for (k in 1:3){
	subset_df <- virus_resp.df %>%
		# Filter samples by desired timeperiod, and by exclusion age greater than  
		# timeperiod (if less, then value for timeperiod cannot be calculated)
		filter(age>365/12*(monthly_val[k]) &
					 	age<=365/12*(monthly_val[k+1]))
	df <- subset_df %>%
		# select only desirable variables
		select(SUBJ_ID, one_of(clin_var_to_count, path_var_to_count)) %>%
		# prepare to sum by subject
		group_by(SUBJ_ID) %>%
		# summate all desirable variables
		summarise_each(funs(sum(., na.rm=T))) %>%
		# Rename with timepoint suffix
		rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_", monthly_var[k])))
	for (j in 1:length(path_var_to_count)){
		df2 <- subset_df %>%
			# Filter samples by positive finding of desired pathogen
			filter_(interp(~ var == 1, var = as.name(path_var_to_count[j]))) %>%
			select(SUBJ_ID, one_of(clin_var_to_count)) %>%
			group_by(SUBJ_ID) %>%
			summarise_each(funs(sum(., na.rm=T))) %>%
			# Rename with pathogen prefix, timepoint suffix
			rename_(.dots = setNames(names(.)[-1], 
															 paste0(path_var_to_count[j],"_",
															 			 names(.)[-1], "_", monthly_var[k])))
		# Incrementally merge with previous dataframe
		df <- merge(df, df2, by="SUBJ_ID", all=T)
		# All NAs replaced with zeroes at this stage
		df[is.na(df)] <- 0
	}
	# Give complete dataframe a useful name
	assign(paste0("virus_resp_sum_noexcl_",monthly_var[k],".df"), df)
}

# Year-based summaries (1,2,3,4,5)
# ----
for (k in 1:5){
	subset_df <- virus_resp.df %>%
		filter(age>365*(k-1) &
					 	age<=365*k)
	if (k <= 3){ 
		df <- subset_df %>%
			select(SUBJ_ID, one_of(clin_var_to_count, path_var_to_count)) %>%
			group_by(SUBJ_ID) %>%
			summarise_each(funs(sum(., na.rm=T))) %>%
			rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_", k)))
		for (j in 1:length(path_var_to_count)){
			df2 <- subset_df %>%
				filter_(interp(~ var == 1, var = as.name(path_var_to_count[j]))) %>%
				select(SUBJ_ID, one_of(clin_var_to_count)) %>%
				group_by(SUBJ_ID) %>%
				summarise_each(funs(sum(., na.rm=T))) %>%
				rename_(.dots = setNames(names(.)[-1], 
																 paste0(path_var_to_count[j],"_",
																 			 names(.)[-1], "_",k)))
			df <- merge(df, df2, by="SUBJ_ID", all=T)
			df[is.na(df)] <- 0
		}
	} else {
		df <- subset_df %>%
			select(SUBJ_ID, one_of(clin_var_to_count)) %>%
			group_by(SUBJ_ID) %>%
			summarise_each(funs(sum(., na.rm=T))) %>%
			rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_", k)))
	}
	assign(paste0("virus_resp_sum_noexcl_y",k,".df"), df)
}

# Year 1 = 425d
# ----

subset_df <- virus_resp.df %>%
	filter(age<=425 & EXCL_AGE>425)
df <- subset_df %>%
	select(SUBJ_ID, one_of(clin_var_to_count, path_var_to_count)) %>%
	group_by(SUBJ_ID) %>%
	summarise_each(funs(sum(., na.rm=T))) %>%
	rename_(.dots = setNames(names(.)[-1], paste0(names(.)[-1], "_425")))
for (j in 1:length(path_var_to_count)){
	df2 <- subset_df %>%
		filter_(interp(~ var == 1, var = as.name(path_var_to_count[j]))) %>%
		select(SUBJ_ID, one_of(clin_var_to_count)) %>%
		group_by(SUBJ_ID) %>%
		summarise_each(funs(sum(., na.rm=T))) %>%
		rename_(.dots = setNames(names(.)[-1], 
														 paste0(path_var_to_count[j],"_",
														 			 names(.)[-1], "_425")))
	df <- merge(df, df2, by="SUBJ_ID", all=T)
	df[is.na(df)] <- 0
}
virus_resp_sum_noexcl_y425.df <- df

# Exclusion ages
# ----
for (k in 1:3){
	df <- eval(parse(text=paste0("virus_resp_sum_",monthly_var[k],".df")))
	df <- merge(df, excl_age.df, by="SUBJ_ID", all=T)
	for (i in 1:nrow(df)) {
		if(df[i,]$EXCL_AGE <= 365/12*monthly_val[k+1]) {
			df[i,-which(colnames(df) %in% c("SUBJ_ID","EXCL_AGE"))] <- NA
		} else {
			df[i,which(is.na(df[i,]))] <- 0
		}
	}
	assign(paste0("virus_resp_sum_noexcl_",monthly_var[k],".df"), df)
}

for (k in 1:5){
	df <- eval(parse(text=paste0("virus_resp_sum_y",k,".df")))
	df <- merge(df, excl_age.df, by="SUBJ_ID", all=T)
	for (i in 1:nrow(df)) {
		if(df[i,]$EXCL_AGE <=365*k) {
			df[i,-which(colnames(df) %in% c("SUBJ_ID","EXCL_AGE"))] <- NA
		} else {
			df[i,which(is.na(df[i,]))] <- 0
		}
	}
	assign(paste0("virus_resp_sum_noexcl_y",k,".df"), df)
}

df <- virus_resp_sum_noexcl_y425.df
df <- merge(df, excl_age.df, by="SUBJ_ID", all=T)
for (i in 1:nrow(df)) {
	if(df[i,]$EXCL_AGE <=425) {
		df[i,-which(colnames(df) %in% c("SUBJ_ID","EXCL_AGE"))] <- NA
	} else {
		df[i,which(is.na(df[i,]))] <- 0
	}
}
virus_resp_sum_noexcl_y425.df <- df

# merge all
virus_resp_sum_noexcl.df <- Reduce(function(...) 
	merge(..., by=c("SUBJ_ID"), all=T),
	list(virus_resp_sum_noexcl_0_2m.df, virus_resp_sum_noexcl_2_6m.df, 
			 virus_resp_sum_noexcl_6_12m.df, virus_resp_sum_noexcl_y1.df,
			 virus_resp_sum_noexcl_y2.df, virus_resp_sum_noexcl_y3.df,
			 virus_resp_sum_noexcl_y4.df, virus_resp_sum_noexcl_y5.df,
			 virus_resp_sum_noexcl_y425.df))

virus_resp_sum_noexcl.df <- 
	merge(virus_resp_sum_noexcl.df, with(virus_resp_sum.df, data.frame(SUBJ_ID, EXCL_AGE)), 
				by="SUBJ_ID", all=T)
with(virus_resp_sum_noexcl.df[
	which(virus_resp_sum_noexcl.df$SUBJ_ID %in% virus_resp_sum.df[which(is.na(virus_resp_sum.df$lri_1)),]$SUBJ_ID),],
	data.frame(SUBJ_ID, lri_1, EXCL_AGE)) %>% filter(EXCL_AGE>=365)
with(virus_resp_sum_noexcl.df[
	which(virus_resp_sum_noexcl.df$SUBJ_ID %in% virus_resp_sum.df[which(is.na(virus_resp_sum.df$lri_2)),]$SUBJ_ID),],
	data.frame(SUBJ_ID, lri_2, EXCL_AGE)) %>% filter(EXCL_AGE>=365*2)
with(virus_resp_sum_noexcl.df[
	which(virus_resp_sum_noexcl.df$SUBJ_ID %in% virus_resp_sum.df[which(is.na(virus_resp_sum.df$lri_2)),]$SUBJ_ID),],
	data.frame(SUBJ_ID, lri_2, EXCL_AGE)) %>% filter(EXCL_AGE>=365*3)

new_excl_ari.df <- read.csv("data/new ari variables.csv", header=T)
virus_resp_sum_noexcl.df <-
	merge(virus_resp_sum_noexcl.df, with(new_excl_ari.df, data.frame(SUBJ_ID = SUBJID, Exclusion.age)), 
				by="SUBJ_ID", all=T)