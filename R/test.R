# dir.create('d:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization',
#   showWarnings = FALSE, recursive = TRUE)
# file.copy(from = 'D:/project/00_zhulab/01_metdna2/00_data/20220412_compare_CAMERA_KGMN/KGMN/MetDNA2_pos/00_annotation_table/00_intermediate_data/table_identification',
#   to = 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization', overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
# file.copy(from = 'D:/project/00_zhulab/01_metdna2/00_data/20220412_compare_CAMERA_KGMN/KGMN/MetDNA2_pos/03_annotation_credential/ms2_data.RData',
#   to = 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization', overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
# file.copy(from = 'D:/project/00_zhulab/01_metdna2/00_data/20220412_compare_CAMERA_KGMN/KGMN/MetDNA2_pos/03_annotation_credential/ms2_data.msp',
#   to = 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization', overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
#
# network_include_all_level <- constructNetwork2(seed_table = seed_table_pos,
#   annotation_table = table_identification_pos,
#   network_emrn = network_emrn,
#   ms2 = ms2_pos,
#   max_reaction_step = 3,
#   ms2_link = 'hybrid')


# file.copy(from = 'D:/project/00_zhulab/01_metdna2/00_data/20220412_compare_CAMERA_KGMN/KGMN/MetDNA2_pos/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData',
#   to = 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization', overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
#
# file.copy(from = 'D:/project/00_zhulab/01_metdna2/00_data/20220412_compare_CAMERA_KGMN/KGMN/MetDNA2_pos/00_annotation_table/00_intermediate_data/list_peak_group',
#   to = 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization', overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
#
# file.copy(from = 'D:/project/00_zhulab/01_metdna2/00_data/20220412_compare_CAMERA_KGMN/KGMN/MetDNA2_pos/00_annotation_table/00_intermediate_data/peak_group_id_table',
#   to = 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization', overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)


# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/00_annotation_table/00_intermediate_data/list_peak_group')
# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/00_annotation_table/00_intermediate_data/peak_group_id_table')
# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/00_annotation_table/00_intermediate_data/table_identification')

################################################################################
# test 46std -------------------------------------------------------------------

# # load packages
# library(CHNOSZ)
# library(dplyr)
# library(MetDNA2Vis)
#
# # set working directory
# # setwd('D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/Demo_MetDNA2_NIST_urine_pos/06_visualization/')

# # Export global networks
# construct network 1
# reconstructNetwork1(dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220613_S9_fraction_ms2/pos/06_visualization/',
#                     is_unknown_annotation = TRUE,
#                     mode = '46std')

# # construct network 2
# annotation_table <- reformatTable1(dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220613_S9_fraction_ms2/pos/06_visualization/')
# reconstructNetwork2(annotation_table = annotation_table,
#                     dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220613_S9_fraction_ms2/pos/06_visualization/',
#                     is_unknown_annotation = TRUE,
#                     mode = '46std')
#
# construct network 3
# reconstructNetwork3(dir_path = 'G:/00_projects/03_MetDNA2/10_project/MetDNA2_project/Data/20220613_S9_fraction_ms2/pos/06_visualization/')

# # Export subnetworks -----------------------------------------------------------
# # network 1 of unknown peak subnetwork
# # Note: the folder_output should keep same among different layer subnetworks
# retrieveSubNetwork1(centric_met = c('C00082', 'KeggExd000923'),
#                     is_unknown_annotation = TRUE,
#                     folder_output = c('M182T541_M262T526'))
#
#
# # network 2 of unknown peak subnetwork
# retrieveSubNetwork2(from_peak = 'M182T541',
#                     end_peak = 'M262T526',
#                     folder_output = c('M182T541_M262T526'))
#
# # network 3 of unknown peak subnetwork
# retrieveSubNetwork3(base_peaks = c('M182T541', 'M262T526'),
#                     base_adducts = c('[M+H]+', '[M+H]+'),
#                     folder_output = c('M182T541_M262T526'))
#
#
# # merge subnetwork
# mergeSubnetwork(from_peak = 'M182T541',
#                 end_peak = 'M262T526',
#                 folder_output = 'M182T541_M262T526')


