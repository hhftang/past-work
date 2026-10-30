#' ---
#' title: "Dataset of immunological variables from CAS"
#' author: "Howard Tang"
#' date: "2016-09-23"
#' output: pdf_document
#' ---
#' 
#' This protocol loads all available CAS datasets relating to immune variables, 
#' and generates variables in the dataframe "immunol_cyto_all.df".
#' 
#' Libraries
#' ----
#' None
#'
#' Internal functions
#' ----
#' None
#'  
#' Immunology data
#' ----
#' 

# Load immunology tables
# ----

# Load immunol file for ages 0 to 5
immunol_5.df <- read.csv("data/CAS_Immunology_db_310516.csv", header=T)

# Load immunol file for age 10
immunol_10.df <- read.csv("data/asthma140828.csv", header=T)
# Merge
immunol.df <- merge(immunol_5.df, immunol_10.df, by="SUBJ_ID", all=T)

# Cytokine response data
# ----
# CMI simulation data
# Given values are deltas
# Negative delta values (i.e. case > control) are censored to some very small 
# positive value near zero (0.000001 for mRNA, 0.01 for protein).
# Is this actually appropriate?!

cmi_y0.df <- read.csv("data/2016-03-23_CMI_y0.csv",header=T,sep=",")
cmi_y0.5.df <- read.csv("data/2016-03-23_CMI_y0.5.csv",header=T,sep=",")
cmi_y1.df <- read.csv("data/2016-03-23_CMI_y1.csv",header=T,sep=",")
cmi_y2.df <- read.csv("data/2016-03-23_CMI_y2.csv",header=T,sep=",")
cmi_y3.df <- read.csv("data/2016-03-23_CMI_y3.csv",header=T,sep=",")
cmi_y5.df <- read.csv("data/2016-03-23_CMI_y5.csv",header=T,sep=",")

cmi.df <- merge(cmi_y0.df, cmi_y0.5.df, by="SUBJ_ID", all=T)
cmi.df <- merge(cmi.df, cmi_y1.df, by="SUBJ_ID", all=T)
cmi.df <- merge(cmi.df, cmi_y2.df, by="SUBJ_ID", all=T)
cmi.df <- merge(cmi.df, cmi_y3.df, by="SUBJ_ID", all=T)
cmi.df <- merge(cmi.df, cmi_y5.df, by="SUBJ_ID", all=T)

immunol_cmi.df <- merge(immunol.df, cmi.df, by="SUBJ_ID", all=T)

# Poly-IC stimulation data
# ------
# Stick to deltas only
polyic_y4_10.df <- read.csv("data/2016-03-23_polyic_y4-10.csv",header=T,sep=",")
immunol_cyto_all.df <- merge(immunol_cmi.df, polyic_y4_10.df, by="SUBJ_ID", all=T)

# dob vs. DOB.x vs DOB.y
immunol_cyto_all.df$dob <- as.Date(immunol_cyto_all.df$dob, format="%d/%m/%Y") 
immunol_cyto_all.df$DOB.x <- as.Date(immunol_cyto_all.df$DOB.x, format="%d/%m/%y") 
immunol_cyto_all.df$DOB.y <- as.Date(immunol_cyto_all.df$DOB.y, format="%d/%m/%Y")
with(immunol_cyto_all.df, data.frame(SUBJ_ID, dob, DOB.x, dob - DOB.x))
with(immunol_cyto_all.df, data.frame(SUBJ_ID, dob, DOB.x, dob - DOB.y))
with(immunol_cyto_all.df, sum(dob - DOB.x, na.rm=T))
with(immunol_cyto_all.df, sum(dob - DOB.y, na.rm=T))
# Use dob only
immunol_cyto_all.df$DOB.x <- NULL
immunol_cyto_all.df$DOB.y <- NULL

# Sex vs. gender
with(immunol_cyto_all.df, data.frame(SUBJ_ID, sex, GENDER.x, sex - GENDER.x))
with(immunol_cyto_all.df, data.frame(SUBJ_ID, sex, GENDER.y, sex - GENDER.y))
with(immunol_cyto_all.df, sum(sex - GENDER.x, na.rm=T))
with(immunol_cyto_all.df, sum(sex - GENDER.y, na.rm=T))
# Use sex only
immunol_cyto_all.df$GENDER.x <- NULL
immunol_cyto_all.df$GENDER.y <- NULL

save.image("output/2016-09-23_immuno_var_gen.Rdata")
