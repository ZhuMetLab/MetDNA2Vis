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
