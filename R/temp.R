################################################################################
# #20220603 add network data ---------------------------------------------------
#
#
# etwork_emrn <- emrn_rpair$step2;rm(emrn_rpair);gc()
# load('D:/project/00_zhulab/01_metdna2/00_code/PackageFile/mrn/210329/emrn_rpair_210329.RData')
#
# network_emrn_step2 <- emrn_rpair$step2
# network_emrn_step0 <- emrn_rpair$step0
# usethis::use_data(network_emrn_step0, overwrite = TRUE)
# usethis::use_data(network_emrn_step2, overwrite = TRUE)


################################################################################
# 20220604 add compound info --------------------------------------------------
#
# load('D:/project/00_zhulab/01_metdna2/00_code/PackageFile/library/20210629_update_emrn/cpd_emrn_210629.RData')
# usethis::use_data(cpd_emrn, overwrite = TRUE)
#
# load('D:/project/00_zhulab/01_metdna2/00_code/PackageFile/library/20210615_update_cpd_info/zhuMetlib_210615.RData')
# cpd_lib <- zhuMetlib$meta$compound
#
# usethis::use_data(cpd_lib, overwrite = TRUE)

################################################################################
