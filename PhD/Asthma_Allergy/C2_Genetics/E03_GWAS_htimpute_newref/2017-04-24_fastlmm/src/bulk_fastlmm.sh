#!/bin/bash

# phenotypes re-sorted in order of priority due to lengthy duration of analysis

pheno_list=( 
	# ari_1_bin ari_1_num ari_2_bin ari_2_num ari_3_bin ari_3_num \
	crwz_1_bin crwz_2_bin crwz_3_bin crwz_4_bin crwz_5_bin \
	cras_3_bin cras_4_bin cras_5_bin \
	# wlri_1_bin \
	wlri_1_num wlri_2_bin wlri_2_num wlri_3_bin wlri_3_num \
	lri_1_bin lri_1_num lri_2_bin lri_2_num lri_3_bin lri_3_num \
	hdmige_6m_bin hdmige_6m_num hdmige_1_bin hdmige_1_num hdmige_2_bin hdmige_2_num hdmige_3_bin hdmige_3_num \
	npem_clust_1_bin npem_clust_3_bin npem_clust_num \
	# the other phenotypes below are not as important: \
	flri_1_bin flri_1_num flri_2_bin flri_2_num flri_3_bin flri_3_num \
	uri_1_bin uri_1_num uri_2_bin uri_2_num uri_3_bin uri_3_num \
	hdmg4_6m_bin hdmg4_6m_num hdmg4_1_bin hdmg4_1_num hdmg4_2_bin hdmg4_2_num hdmg4_3_bin hdmg4_3_num \
	pntige_6m_bin pntige_6m_num pntige_1_bin pntige_1_num pntige_2_bin pntige_2_num pntige_3_bin pntige_3_num \
	pntg4_6m_bin pntg4_6m_num pntg4_1_bin pntg4_1_num pntg4_2_bin pntg4_2_num pntg4_3_bin pntg4_3_num \
	hdmwl_6m_bin hdmwl_6m_num hdmwl_2_bin hdmwl_2_num \
	hrva_lri_1_bin hrva_lri_1_num hrva_lri_2_bin hrva_lri_2_num hrva_lri_3_bin hrva_lri_3_num \
	hrvc_lri_1_bin hrvc_lri_1_num hrvc_lri_2_bin hrvc_lri_2_num hrvc_lri_3_bin hrvc_lri_3_num \
	influ_lri_1_bin influ_lri_1_num influ_lri_2_bin influ_lri_2_num influ_lri_3_bin influ_lri_3_num \
	rsv_lri_1_bin rsv_lri_1_num rsv_lri_2_bin rsv_lri_2_num rsv_lri_3_bin rsv_lri_3_num \
	crez_6m_bin crez_1_bin crez_2_bin crez_3_bin crez_4_bin crez_5_bin \
)

for i in "${pheno_list[@]}"
do 
	./fastlmm_task.sh cas_imp_3 "$i"
	./fastlmm_task_sex.sh cas_imp_3 "$i"
	Rscript man_qq_lamb_script.R "$i"
	Rscript man_qq_lamb_script.R "$i"_sex
done
