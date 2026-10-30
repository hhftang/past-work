write.table(
  with(premaster.df,
    data.frame(SUBJ_ID, slri_1, slri_2, slri_3, slri_4, slri_5)
  ), "output/slri.txt"
)