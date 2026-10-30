#' ---
#' title: "Generate master table for analysis"
#' author: "Howard Tang"
#' date: "2017-06-19"
#' output: pdf_document
#' ---
#' 
#' Libraries
#' ----
#'
desired_packages <- c("dplyr", "ggplot2", "reshape2")
new_packages <- desired_packages[!(desired_packages %in% 
																	 	installed.packages()[,"Package"])]
if(length(new_packages)>0) install.packages(new_packages)
for (i in 1:length(desired_packages)){
	library(desired_packages[i], character.only=T)
}

#' Internal functions
#' ----
#' None
#' 
#' 
#' Load data tables 
#' -----------------------------------------------------------------------------
#' 

load("output/2016-09-26_vir_respchart_var_gen.Rdata")

# Keep premaster.df only
rm(list=setdiff(ls(), "premaster.df"))

premaster.df <- rename(premaster.df, 
											 hdmigg_10 = hdmigg10,
											 catigg_10 = catigg10,
											 grassigg_10 = grassigg10,
											 
											 # cd_ic10, cd_ic1b, cd_ic6, cd_ic8, cd_icg, cd_icip, cd_icmp, cd_ict,
											 # y4_ic10, y4_ic1b, y4_ic6, y4_ic8, y4_icg, y4_icip, y4_icmp, y4_ict,
											 # y10_ic10, y10_ic1b, y10_ic6, y10_ic8, y10_icg, X.y10_icip, y10_icmp, y10_ict,
											 # Swap variable name and timepoint; change X.y10_icip
											 
											 ic10_cd = cd_ic10,
											 ic1b_cd = cd_ic1b,
											 ic6_cd = cd_ic6,
											 ic8_cd = cd_ic8,
											 icg_cd = cd_icg,
											 icip_cd = cd_icip, 
											 icmp_cd = cd_icmp, 
											 ict_cd = cd_ict, 
											 
											 ic10_4 = y4_ic10,
											 ic1b_4 = y4_ic1b,
											 ic6_4 = y4_ic6,
											 ic8_4 = y4_ic8,
											 icg_4 = y4_icg,
											 icip_4 = y4_icip, 
											 icmp_4 = y4_icmp, 
											 ict_4 = y4_ict, 
											 
											 ic10_10 = y10_ic10,
											 ic1b_10 = y10_ic1b,
											 ic6_10 = y10_ic6,
											 ic8_10 = y10_ic8,
											 icg_10 = y10_icg,
											 icip_10 = X.y10_icip, 
											 icmp_10 = y10_icmp, 
											 ict_10 = y10_ict,
											 
											 # currentwhz.10, currentasthma.10
											 # Change to crwz10, cras_10
											 crwz10 = currentwhz.10,
											 cras_10 = currentasthma.10
)

# Use new definitions of uri, lri from virus_respchart.df instead of immunol
# dataset variables

#' Variable quality control and generate old master table
#' -----------------------------------------------------------------------------
# Trim to non-redundant and useful variables
# Create table of column names for easy reference
col_premaster <- colnames(premaster.df)
write(col_premaster, "output/col_premaster.txt")

# Trust viral dataset for pathogen counts
# Trust respchart dataset for clinical infection counts
# Otherwise, use immunol dataset for all other variables

# Ignore new ari phenotypes for now

# Generate old master table
# ---
# Trim, re-order, and subset out variables 
master.df <- premaster.df %>%
	select(
		# SUBJ_ID, Sex (not GENDER),
		SUBJ_ID, sex,
		
		# General infection characteristics
		# ari_1, ari_2, ari_3, ari_4, ari_5,
		uri_1, uri_2, uri_3, uri_4, uri_5,
		lri_1, lri_2, lri_3, lri_4, lri_5,
		wlri_1, wlri_2, wlri_3, wlri_4, wlri_5,
		flri_1, flri_2, flri_3, flri_4, flri_5,
		
		# Number of pathogen-related infections, up to age 3
		# rsv_1, influ_1, hrv_1, hrva_1, hrvb_1, hrvc_1, 
		# rsv_ari_1, 
		rsv_uri_1, rsv_lri_1, rsv_wlri_1, rsv_flri_1,
		# influ_ari_1, 
		influ_uri_1, influ_lri_1, influ_wlri_1, influ_flri_1,
		# hrv_ari_1, 
		# hrv_uri_1, hrv_lri_1, hrv_wlri_1, hrv_flri_1,
		# hrva_ari_1, 
		hrva_uri_1, hrva_lri_1, hrva_wlri_1, hrva_flri_1,
		# hrvb_ari_1, 
		hrvb_uri_1, hrvb_lri_1, hrvb_wlri_1, hrvb_flri_1,
		# hrvc_ari_1, 
		hrvc_uri_1, hrvc_lri_1, hrvc_wlri_1, hrvc_flri_1,
		# rsv_2, influ_2, hrv_2, hrva_2, hrvb_2, hrvc_2, 
		# rsv_ari_2, 
		rsv_uri_2, rsv_lri_2, rsv_wlri_2, rsv_flri_2,
		# influ_ari_2, 
		influ_uri_2, influ_lri_2, influ_wlri_2, influ_flri_2,
		# hrv_ari_2, 
		# hrv_uri_2, hrv_lri_2, hrv_wlri_2, hrv_flri_2,
		# hrva_ari_2, 
		hrva_uri_2, hrva_lri_2, hrva_wlri_2, hrva_flri_2,
		# hrvb_ari_2, 
		hrvb_uri_2, hrvb_lri_2, hrvb_wlri_2, hrvb_flri_2,
		# hrvc_ari_2, 
		hrvc_uri_2, hrvc_lri_2, hrvc_wlri_2, hrvc_flri_2,
		# rsv_3, influ_3, hrv_3, hrva_3, hrvb_3, hrvc_3, 
		# rsv_ari_3, 
		rsv_uri_3, rsv_lri_3, rsv_wlri_3, rsv_flri_3,
		# influ_ari_3, 
		influ_uri_3, influ_lri_3, influ_wlri_3, influ_flri_3,
		# hrv_ari_3, 
		# hrv_uri_3, hrv_lri_3, hrv_wlri_3, hrv_flri_3,
		# hrva_ari_3, 
		hrva_uri_3, hrva_lri_3, hrva_wlri_3, hrva_flri_3,
		# hrvb_ari_3, 
		hrvb_uri_3, hrvb_lri_3, hrvb_wlri_3, hrvb_flri_3,
		# hrvc_ari_3, 
		hrvc_uri_3, hrvc_lri_3, hrvc_wlri_3, hrvc_flri_3,			
		
		# Allergen IgE (need to log10 later)
		# Year 10 IgE dataset also had variables with "r" prefix -> raw values
		# prior to applying limit of detection (<0.030 -> 0.015). These were ignored.
		totige6m, hdmige6m, catige6m, pntige6m, couige6m, ryeige6m, mouige6m, phinfige6m,
		totige1, hdmige1, catige1, pntige1, couige1, ryeige1, mouige1, phinfige1,
		totige2, hdmige2, catige2, pntige2, couige2, ryeige2, mouige2, phinfige2,
		totige3, hdmige3, catige3, pntige3, couige3, ryeige3, mouige3, phinfige3,
		totige4, hdmige4, catige4, pntige4, couige4, ryeige4, mouige4, phinfige4,
		totige5, hdmige5, catige5, pntige5, couige5, ryeige5, mouige5, phadige5, fmxige5, 
		totige10, hdmige10, catige10, pntige10, couige10, ryeige10, mouige10, phige10, fmxige10,
		
		# Allergen IgG4 (need to log10)
		hdmg4_6m, catg4_6m, pntg4_6m, coug4_6m, ryeg4_6m, moug4_6m, phinfg4_6m, 
		hdmg4_1, catg4_1, pntg4_1, coug4_1, ryeg4_1, moug4_1, phinfg4_1, 
		hdmg4_2, catg4_2, pntg4_2, coug4_2, ryeg4_2, moug4_2, phinfg4_2, 
		hdmg4_3, catg4_3, pntg4_3, coug4_3, ryeg4_3, moug4_3, phinfg4_3, 
		hdmg4_4, catg4_4, pntg4_4, coug4_4, ryeg4_4, moug4_4, phinfg4_4, 
		hdmg4_5, catg4_5, pntg4_5, coug4_5, ryeg4_5, moug4_5, phadg4_5, fmxg4_5,
		
		# Allergen IgG (need to log10)
		hdmigg_1, catigg_1, grassigg_1, 
		hdmigg_2, catigg_2, grassigg_2, 
		hdmigg_3, catigg_3, grassigg_3, 
		hdmigg_5, catigg_5, grassigg_5,
		hdmigg_10, catigg_10, grassigg_10,
		
		# Anti-peptide IgE (need to log10)--- removed on Kat's request
		# Anti-peptide IgG1 (need to log10) --- removed on Kat's request
		# Anti-peptide IgG4 (need to log10) --- removed on Kat's request
		
		# Allergen SPT wheal size
		# Ignore binary SPT results (derived variables based on size >1-2mm)
		hiswl6m, hdmwl6m, catwl6m, ryewl6m, altwl6m, aspwl6m, cmwl6m, ewwl6m,
		hiswl2, hdmwl2, catwl2, ryewl2, altwl2, aspwl2, cmwl2, ewwl2, 
		hiswl5, hdmwl5, catwl5, ryewl5, altwl5, aspwl5, cmwl5, ewwl5, pntwl5,
		hiswl10, hdmwl10, catwl10, couwl10, ryewl10, grwl10, mouwl10, cmwl10, ewwl10, pntwl10, 
		
		# FHx atopy defined by SPT
		mat_atopy, pat_atopy,
		
		# FHx atopy dx by dr
		mat_had, pat_had,
		
		# FHx asthma dx by dr
		mat_ast, pat_ast,
		
		# FHx non-asthmatic atopic disease dx by dr (derived)
		matfh_xast, patfh_xast, 
		# In newer tables, maybe I should replace mat_had, pat_had with 
		# matfh_xast, patfh_xast
		
		# Smoke exposure
		smoke1, smoke2, smoke3, smoke4, smoke5, smoke10,
		
		# Exposure to childcare
		care1, care2, care3, care4,
		
		# Exposure to cat
		# for years 1 to 5, 2 represents outdoor exposure; 1 indoor;
		# Need to change all 2s to 1s later, to simply get "exposure" (consistent
		# with age 10)
		crcatexp1, crcatexp2, crcatexp3, crcatexp4, crcatexp5, crcatexp10,
		
		# Exposure to other animals (age 10 only)
		crdogexp10, crfurpet10,
		
		# Number of children
		nochdn1, nochdn2, nochdn3, nochdn4, nochdn5, nochdn10,
		
		# Number of older children (not derived) in same house - should remain relatively
		# constant with time
		oldchdbirth, oldchd1, oldchd2, oldchd3, oldchd4, oldchd5, oldchd10,
		
		# Birth order
		bthodr,
		
		# FACs cell count and HLA-DR MFI
		bas6m, bas1, bas2, bas5, 
		pDC6m, pDC1, pDC2, pDC5,
		mDC6m, mDC1, mDC2, mDC5, 
		pDCmfi6m, pDCmfi1, pDCmfi2, pDCmfi5, 
		mDCmfi6m, mDCmfi1, mDCmfi2, mDCmfi5,
		
		# TRECS (need to log10)
		trecs6m, trecs1, trecs2, trecs5,
		
		# Body size
		# Ignore BMI (derived)
		ht3, ht4, ht5, ht10,
		wt3, wt4, wt5, wt10,
		
		# Maternal breast milk fatty acid content
		# Ignored fatty acid ratios (derived)
		n3_6w, n6_6w, n3_6m, n6_6m,
		
		# Maternal breast milk IgA and protein content
		iga6w, iga6m, 
		lys6w, lys6m, 
		prot6w, prot6m, 
		lact6w, lact6m, 
		
		# Maternal Vitamin D
		vitdmat,
		
		# Bacterial average abundance --- removed on Shumei's request
		# Bacterial MPG --- removed on Shumei's request
		
		# mRNA with HDM stimulation of CBMCs (need to log10)
		mhdm13_cd, mhdm4_cd, mhdm5_cd, mhdm9_cd, mhdmg_cd,
		mhdm13_6m, mhdm4_6m, mhdm4r_6m, mhdm5_6m, mhdm9_6m, mhdmg_6m,
		mhdm13_1, mhdm4_1, mhdm4r_1, mhdm5_1, mhdm9_1, mhdmg_1, 
		mhdm13_2, mhdm4_2, mhdm4r_2, mhdm5_2, mhdm9_2, mhdmg_2, 
		mhdm13_3, mhdm4_3, mhdm4r_3, mhdm5_3, mhdm9_3, mhdmg_3,
		mhdm13_5, mhdm4_5, mhdm4r_5, mhdm5_5, mhdm9_5, mhdmg_5,
		
		# Expression for CISH, DACT1 and PLXDC1
		# See http://www.jimmunol.org/content/176/8/4766.full for why these
		# genes are included
		mhdmcish_2, mhdmdac_2, mhdmplx_2,
		
		# Protein with CBMC stimulation (need to log10)
		# ignored ratios rpha, rseb etc.
		# Cord blood:
		hdm10_cd, hdm13_cd, hdm5_cd, hdmg_cd,
		cat10_cd, cat13_cd, cat5_cd, catg_cd, 
		# Fel d1 cat allergen
		fel10_cd, fel13_cd, fel5_cd, felg_cd, 
		# Lipopolysaccharide (LPS)
		lpsg10_cd, lpsg12_cd,
		# ovalbumin
		ova10_cd, ova13_cd, ova5_cd, ovag_cd, 
		# phytohaemagluttinin
		pha10_cd, pha13_cd, pha5_cd, phag_cd, 
		# PPD, purified protein derivative of TB
		ppd10_cd, ppd13_cd, ppd5_cd, ppdg_cd, 
		# Staphylococcal enterotoxin B
		seb10_cd, seb13_cd, seb5_cd, sebg_cd, 
		
		# Age 6m:
		hdm10_6m, hdm13_6m, hdm5_6m, hdmg_6m,
		# Arah2 peanut allergen
		ara10_6m, ara13_6m, ara5_6m, arag_6m,
		lpsg10_6m, lpsg12_6m, 
		ova10_6m, ova13_6m, ova5_6m, ovag_6m, 
		pha10_6m, pha13_6m, pha5_6m, phag_6m, 
		# RSV
		rsv10_6m, rsv13_6m, rsv5_6m, rsvg_6m, 
		seb10_6m, seb13_6m, seb5_6m, sebg_6m, 
		# Tetanus toxoid
		tt10_6m, tt13_6m, tt5_6m, ttg_6m, 
		
		# Age 1: 
		# ***t_1 etc. = TNF-alpha expression
		hdm10_1, hdm13_1, hdm5_1, hdmg_1, hdmt_1, 		
		ara10_1, ara13_1, ara5_1, arag_1, arat_1, 
		lpsg10_1, lpsg12_1,
		ova10_1, ova13_1, ova5_1, ovag_1, ovat_1, 
		pha10_1, pha13_1, pha5_1, phag_1, phat_1,
		rsv10_1, rsv13_1, rsv5_1, rsvg_1, rsvt_1, 
		seb10_1, seb13_1, seb5_1, sebg_1, sebt_1, 
		tt_t_1, tt10_1, tt13_1, tt5_1, ttg_1,
		
		# Age 2:
		hdm10_2, hdm13_2, hdm5_2, hdmg_2,
		# PolyIC response
		ic10_2, icg_2, ict_2,
		lpsg10_2, lpsg_g_2, lpsgt_2,
		# Lipoteichoic acid (LTA)
		lta10_2, ltag_2, ltat_2,
		ova10_2, ova13_2, ova5_2, ovag_2, 
		pha10_2, pha13_2, pha5_2, phag_2, 
		# PMA/ionomycin
		pi10_2, pi13_2, pi5_2, pig_2, 
		pnt10_2, pnt13_2, pnt5_2, pntg_2, 
		# Pansorbin		
		ps10_2, psg_2, pst_2,
		rye10_2, rye13_2, rye5_2, ryeg_2, 
		seb10_2, seb13_2, seb5_2, sebg_2, 
		tt10_2, tt13_2, tt5_2, ttg_2,
		
		# Age 3:
		hdm10_3, hdm13_3, hdm5_3, hdmg_3, 
		# CpG
		cpg10_3, cpgg_3, cpgt_3,
		# Flagellin
		fla10_3, flag_3, flat_3, 
		# Gardiquimod
		ga10_3, gag_3, gat_3, 
		ic10_3, icg_3, ict_3, 
		# Imiquimod
		imi10_3, imig_3, imit_3,
		lpsg10_3, lpsg_g_3, lpsgt_3, 
		lta10_3, ltag_3, ltat_3, 
		# Measles lysate
		msl10_3, msl13_3, msl5_3, mslg_3, 
		pha10_3, pha13_3, pha5_3, phag_3, 
		pi10_3, pi13_3, pi5_3, pig_3, 
		ps10_3, psg_3, pst_3, 
		rye10_3, rye13_3, rye5_3, ryeg_3, 
		seb10_3, seb13_3, seb5_3, sebg_3, 
		tt10_3, tt13_3, tt5_3, ttg_3,
		
		# Age 5:
		hdm10_5, hdm13_5, hdm5_5, hdmg_5, 
		fel10_5, fel13_5, fel5_5, felg_5, 
		ic10_5, icg_5, ict_5, 
		lpsg10_5, lpsg_g_5, lpsgt_5, 
		lta10_5, ltag_5, ltat_5, 
		msl10_5, msl13_5, msl5_5, mslg_5, 
		ova10_5, ova13_5, ova5_5, ovag_5, 
		pha10_5, pha13_5, pha5_5, phag_5, 
		pi10_5, pi13_5, pi5_5, pig_5, 
		pnt10_5, pnt13_5, pnt5_5, pntg_5,
		ps10_5, psg_5, pst_5, 
		rye10_5, rye13_5, rye5_5, ryeg_5, 
		seb10_5, seb13_5, seb5_5, sebg_5, 
		tt10_5, tt13_5, tt5_5, ttg_5,
		
		# PolyIC stimulation of CBMC (need to log10)
		ic10_cd, ic1b_cd, ic6_cd, ic8_cd, icg_cd, icip_cd, icmp_cd, ict_cd, 
		ic10_4, ic1b_4, ic6_4, ic8_4, icg_4, icip_4, icmp_4, ict_4,
		ic10_10, ic1b_10, ic6_10, ic8_10, icg_10, icip_10, icmp_10, ict_10,
		
		# Current wheeze
		crwz1, crwz2, crwz3, crwz4, crwz5, crwz10,
		
		# Current asthma
		# Ignored currentasthma10wm, currentasthma10worm, and DDX10 
		# I do not know what they mean!
		cras_3, cras_4, cras_5, cras_10, 
		
		# Current (atopic) rhinoconjunctivitis (note: semi-derived)
		arc5,
		
		# Current eczema
		crez6m, crez1, crez2, crez3, crez4, crez5
	)

#' Ignored variables
#' 
#' - derived variables for total number of events from age 1 to 3 
#' - frequency of infections within first year (<6m)
#' - binary infection variables (also derived)
#' - log Ig values (will recalculate manually)
#' - binary IgE results (positive = >0.35) (derived) 
#' - binary SPT results (derived variables based on size >1-2mm)
#' - log TREC values (will calculate manually)
#' - collection dates (not relevant)
#' - FHx from any one parent (derived; parents treated separately)
#' - non-asthmatic atopic disease dx by dr (derived)
#' - atopic disease as defined positive test (SPT, IgE) (derived)
#' - BMI (derived)
#' - fatty acid ratios (derived)
#' - transient wheeze (derived)
#'
#' Potential problematic variables
#' 
#' - total number of infections > number of specific infection
#' - ARI = URI + LRI or "URI OR LRI"???
#' - Total IgE > sum of measured IgE specifities, or value of one IgE specifity
#' - Number of children in household > number of older children in household
#' - MFI From flow cytometry samples may depend on cell count
#' - Total protein content in breast milk > content of individual protein 
#' constituents
#' - arc5 (atopic rhinoconjunctivitis at age 5) is based on SPT results for atopy
#' - Current furry pet exposure > current cat / dog exposure (age 10 only)
#' - Current treatment of protein/mRNA expression data is suboptimal 
#' (censoring)? - should get raw data
#'

# Create table of column names for easy reference
col_master <- colnames(master.df)
write(col_master, "output/col_master.txt")

# Apply Log10 to all relevant variables
# -----
# Log10 for the antibody assays
col_antibody <- which(grepl("ige",colnames(master.df)) |
												grepl("g4_", colnames(master.df)) |
												grepl("igg_", colnames(master.df)))
master.df[,col_antibody] <- log10(master.df[,col_antibody])

# Log10 for TRECS
col_trecs <- which(grepl("trecs", colnames(master.df)))
master.df[,col_trecs] <- log10(master.df[,col_trecs])

# Log10 for CMI stimulation mRNA and protein assays
# Convert all zero values to limit for now (0.000001 for mRNA, 0.01 for protein)
col_mhdm <- which(grepl("mhdm", colnames(master.df)))
for (i in col_mhdm){
	master.df[which(master.df[,i]==0),i] <- 0.000001
}
master.df[,col_mhdm] <- log10(master.df[,col_mhdm])

# Warning: for protein - indices not automated; need to be manually-keyed in
col_protein <- seq(354,622)
for (i in col_protein){
	master.df[which(master.df[,i]==0),i] <- 0.01
}
master.df[,col_protein] <- log10(master.df[,col_protein])

# Convert cat exposure to binary variable
# ----
# 0 (none) vs. 1 (indoor) vs. 2 (outdoor) collapsed to 0 (no exposure) vs. 1 (yes)
col_crcatexp <- which(grepl("crcatexp", colnames(master.df)))
for (i in col_crcatexp){
	master.df[which(master.df[,i]==2), i] <- 1
}

# Rename rownames as SUBJ_ID, and convert NaNs to NA
# ----
rownames(master.df) <- master.df$SUBJ_ID
master.df[is.na(master.df)] <- NA

# Make table without subj_id
# ----
# Remove pesky non-numeric columns (e.g. subject IDs)
master_nosub.df <- master.df
master_nosub.df$SUBJ_ID <- NULL

premaster_new.df <- premaster.df %>%
	select(-contains("fari"), -contains("cari"), -contains("uri_new"),
				 -uri1, -uri2, -uri3, -uri4, -uri5,
				 -lri1, -lri2, -lri3, -lri4, -lri5,
				 -flri1, -flri2, -flri3, -flri4, -flri5,
				 -wlri1, -wlri2, -wlri3, -wlri4, -wlri5,
				 -slri1, -slri2, -slri3, -slri4, -slri5,
				 -furi1, -furi2, -furi3, -furi4, -furi5,
				 -ari1, -ari2, -ari3, -ari4, -ari5,
				 -uri1x6m, -lri1x6m, -flri1x6m, -wlri1x6m, -flri1x6m, -slri1x6m, -ari1x6m, -furi1x6m,
				 -uri2x6m, -lri2x6m, -flri2x6m, -wlri2x6m, -flri2x6m, -slri2x6m, -ari2x6m, -furi2x6m,
				 -mlri1, -mlri2, -mlri3, -mlri4, -mlri5,
				 -lri24m, -wri24m, -flri24m, -slri24m, -mlri24m, 
				 -mlri1x6mo, -mlri2x6mo, 
				 -hrv1, -hrv2, -hrv3, -rsv1, -rsv2, -rsv3, 
				 -rsv_a2, -rsv_a3, -rsv_b2, -rsv_b3,
				 -contains("lrihrv"), -contains("lrirsv"), -contains("lriother"),
				 -contains("lrimv"))

# 1x6m variables can be rederived by summing 0_2m and 2_6m
write.table(premaster_new.df, "output/cas_immunology_db_131016.txt", sep="\t")
write.csv(premaster_new.df, "output/cas_immunology_db_131016.csv")

# Save output to .Rdata
save.image("output/2017-06-19_master_table_gen.Rdata")
