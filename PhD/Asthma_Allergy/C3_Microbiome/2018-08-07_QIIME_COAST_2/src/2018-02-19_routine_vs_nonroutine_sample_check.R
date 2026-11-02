# Emergency crosstab

# Routine vs. non-routine
# ----------------------
# Define new variable routine vs. non-routine
log.table.df$routine.state <- factor(log.table.df$visittype, 
																		 levels=c("Well", "Well/Sick", "Sick"), 
																		 labels=c("Routine", "Routine Sick", "Non-routine"))
levels(log.table.df$routine.state) <- c("Routine", "Routine", "Non-routine")

with(log.table.df, prop.table(table(season, routine.state), 2))
with(log.table.df, prop.table(table(season2, routine.state), 2))

# chi-square test
with(log.table.df, chisq.test(table(season, routine.state))) # p < 2.2e-16
with(log.table.df, chisq.test(table(season2, routine.state))) # p = 3.06e-11

# Expected given that routine samples have mostly well, and non-routine mostly sick

# Only unwell samples - compare by routine.state  
# ----------------------
with(log.table.df[which(log.table.df$visittype2 != "Well"),], 
		 prop.table(table(season, routine.state), 2))

with(log.table.df[which(log.table.df$visittype2 != "Well"),], 
		 prop.table(table(season)))
		 
with(log.table.df[which(log.table.df$visittype2 != "Well"),], 
		 chisq.test(table(season, routine.state))) # p = 0.002768


# Only sick samples - compare by routine.state  
# ----------------------
with(log.table.df[which(log.table.df$visittype2 == "Sick"),], 
		 prop.table(table(season, routine.state), 2))

with(log.table.df[which(log.table.df$visittype2 == "Sick"),], 
		 prop.table(table(season)))

with(log.table.df[which(log.table.df$visittype2 == "Sick"),], 
		 chisq.test(table(season, routine.state))) # p = 0.002768

# Only mild samples - compare by routine.state  
# ----------------------
with(log.table.df[which(log.table.df$visittype2 == "Mild"),], 
		 prop.table(table(season, routine.state), 2))

with(log.table.df[which(log.table.df$visittype2 == "Mild"),], 
		 prop.table(table(season)))

with(log.table.df[which(log.table.df$visittype2 == "Mild"),], 
		 chisq.test(table(season, routine.state))) # p = 0.002768


# Mild vs. severe
# ----------------------
with(log.table.df, prop.table(table(season, visittype2), 2))
with(log.table.df, prop.table(table(season2, visittype2), 2))

df <- log.table.df %>% 
	filter(visittype2 != "Well") %>%
	mutate(visittype2 = factor(visittype2, levels=c("Mild", "Sick")))
with(df, chisq.test(table(season, visittype2))) # p =1e-6
with(df, chisq.test(table(season2, visittype2))) # p=0.1402

# Supports argument against combining severe/mild samples at the same time as
# combining routine/non-routine sick samples

# Routine x severity
with(log.table.df, table(season, visittype2, routine.state))

# For routine samples
with(log.table.df[which(log.table.df$routine.state == "Routine"),], 
		 prop.table(table(season, visittype2), 2))

# For non-routine samples
with(df[which(df$routine.state == "Non-routine"),], 
		 prop.table(table(season, visittype2), 2))


# Remove the milds and repeat routine vs. non-routine for unwell samples
# ----------------------

df <- log.table.df %>% 
	filter(visittype2 != "Mild") %>%
	mutate(visittype2 = factor(visittype2, levels=c("Well", "Sick")))
with(df[which(df$visittype2 == "Sick"),], 
		 table(season, routine.state))

with(df[which(df$visittype2 == "Sick"),], 
		 prop.table(table(season, routine.state), 2))
with(df[which(df$visittype2 == "Sick"),], 
		 chisq.test(table(season, routine.state))) # p = 0.8619


# CONSENSUS - remove the milds (or rather ignore milds for well vs. sick comparison)
