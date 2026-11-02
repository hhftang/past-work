# Get averaged GRS - use polar coordinates?!

load("output/2018-07-19_GRS_associations.Rdata")

master_GRS.df <- master_GRS.df %>% 
		 mutate(score_magnitude = (SCORESUM_allergic_std^2 + SCORESUM_asthma_std^2)^0.5,
		 			 score_angle = atan2(y = SCORESUM_allergic_std, x = SCORESUM_asthma_std))

# cluster
master_GRS.df %>% filter(!is.na(group)) %>%
	ggplot(aes(x = SCORESUM_asthma_std, y = SCORESUM_allergic_std)) + 
	geom_point(aes(color = group)) + 
	scale_color_manual(values = quaPalette.npem)

master_GRS.df %>% filter(!is.na(group)) %>%
	ggplot(aes(y = score_magnitude, x = score_angle)) + 
	geom_point(aes(color = group)) + 
	scale_color_manual(values = quaPalette.npem)

master_GRS.df %>% filter(!is.na(group)) %>%
	ggplot(aes(x = group, y = score_angle)) + 
	geom_boxplot(aes(color = group)) + 
	scale_color_manual(values = quaPalette.npem)


t.test(master_GRS.df$score_angle[master_GRS.df$group == 3], 
			 master_GRS.df$score_angle[master_GRS.df$group != 3])
# p = 0.060

# crwz5
master_GRS.df %>% filter(!is.na(group)) %>%
ggplot(aes(x = SCORESUM_asthma_std, y = SCORESUM_allergic_std)) + 
	geom_point(aes(color = factor(crwz5))) + 
	scale_color_manual(values = c("black", "red"))

ggplot(master_GRS.df, 
			 aes(y = score_magnitude, x = score_angle)) + 
	geom_point(aes(color = factor(crwz5))) + 
	scale_color_manual(values = c("black", "red"))

ggplot(master_GRS.df, 
			 aes(x = group, y = score_angle)) + 
	geom_boxplot(aes(color = factor(crwz5))) + 
	scale_color_manual(values = c("black", "red"))

glm.out <- glm(crwz5 ~ score_angle + group, data = master_GRS.df)
summary(glm.out) # p = 0.179 for score_angle, 1.82e-06 for cluster3

df <- master_GRS.df %>% filter(group == 3)
t.test(df$score_angle[df$crwz5 == 1], 
			 df$score_angle[df$crwz5 == 0])
# p = 0.51

# cras_5
master_GRS.df %>% filter(!is.na(group)) %>%
	ggplot(aes(x = SCORESUM_asthma_std, y = SCORESUM_allergic_std)) + 
	geom_point(aes(color = factor(cras_5))) + 
	scale_color_manual(values = c("black", "red"))

ggplot(master_GRS.df, 
			 aes(y = score_magnitude, x = score_angle)) + 
	geom_point(aes(color = factor(cras_5))) + 
	scale_color_manual(values = c("black", "red"))

ggplot(master_GRS.df, 
			 aes(x = group, y = score_angle)) + 
	geom_boxplot(aes(color = factor(cras_5))) + 
	scale_color_manual(values = c("black", "red"))

glm.out <- glm(cras_5 ~ score_angle + group, data = master_GRS.df)
summary(glm.out) # p = 0.73 for score_angle, 4.36e-05 for cluster3

df <- master_GRS.df %>% filter(group == 3)
t.test(df$score_angle[df$cras_5 == 1], 
			 df$score_angle[df$cras_5 == 0])
# p = 0.93
