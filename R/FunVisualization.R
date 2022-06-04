################################################################################
# prepare functions ------------------------------------------------------------

  # copyFiles4Vis --------------------------------------------------------------
#' @title copyFiles4Vis
#' @author Zhiwei Zhou
#' @description copy intermidate data files for visualization
#' @param dir_path the path of KGMN result folder
#' @export
copyFiles4Vis <- function(dir_path = '.') {

  to_folder <- file.path(dir_path, '06_visualization')
  dir.create(to_folder, showWarnings = FALSE, recursive = TRUE)
  file.copy(from = file.path(dir_path,
    '00_annotation_table/00_intermediate_data/table_identification'),
    to = to_folder,
    overwrite = TRUE,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE)

  file.copy(from = file.path(dir_path,
    '03_annotation_credential/ms2_data.RData'),
    to = to_folder,
    overwrite = TRUE,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE)

  file.copy(from = file.path(dir_path,
    '03_annotation_credential/ms2_data.msp'),
    to = to_folder,
    overwrite = TRUE,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE)

  file.copy(from = file.path(dir_path,
    '03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData'),
    to = to_folder,
    overwrite = TRUE,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE)

  file.copy(from = file.path(dir_path,
    '00_annotation_table/00_intermediate_data/list_peak_group'),
    to = to_folder,
    overwrite = TRUE,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE)

  file.copy(from = file.path(dir_path,
    '00_annotation_table/00_intermediate_data/peak_group_id_table'),
    to = to_folder,
    overwrite = TRUE,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE)

}


  # checkFiles4Vis -------------------------------------------------------------
#' @title checkFiles4Vis
#' @author Zhiwei Zhou
#' @description Check required files for visualization
#' @param dir_path '.'
#' @export
checkFiles4Vis <- function(dir_path = '.'){
  cat('Check required files ...\n')

  file_require <- c('table_identification',
    'ms2_data.RData',
    'ms2_data.msp',
    'list_peak_group_annotation_concised.RData',
    'list_peak_group',
    'peak_group_id_table')

  temp <- file_require %in% list.files(dir_path)

  if (sum(temp) < 6) {
    error <- which(!temp) %>% file_require[.]
    stop('Intermidate data files:', paste(error, collapse = ','), ' are missed!\n')
  }

  cat('Check required files: done! \n\n')

}



################################################################################
# reformatTable1 ---------------------------------------------------------------

#' @title reformatTable1
#' @description
#' @author Zhiwei Zhou
#' @export
#' @importFrom magrittr '%>%' '%$%'

# dir_path <- 'd:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization'
reformatTable1 <- function(dir_path = '.'){
  load(file.path(dir_path, 'table_identification'))
  load(file.path(dir_path, 'ms2_data.RData'))
  table_identification <- table_identification %>%
    # mutate(peak_name = paste(peak_name, 'POS', sep = '_')) %>%
    dplyr::mutate(with_ms2 = case_when(peak_name %in% names(raw_msms) ~ TRUE,
      !(peak_name %in% names(raw_msms)) ~ FALSE)) %>%
    dplyr::mutate(id_kegg = case_when(is.na(id_kegg) ~ id_zhulab,
      !is.na(id_kegg) ~ id_kegg)) %>%
    dplyr::mutate(label_seed = paste(peak_name, id_kegg, sep = '_'))

  return(table_identification)
}



################################################################################
# reconstructNetwork1 ----------------------------------------------------------

#' @title reconstructNetwork1
#' @author Zhiwei Zhou
#' @param dir_path path of working directory. Default: '.'
#' @param is_unknown_annotation whether used unknown annotation. Default: TRUE
#' @export

reconstructNetwork1 <- function(dir_path = '.',
  is_unknown_annotation = TRUE) {
  cat('Extract network1 from KMRN... \n')

  if (is_unknown_annotation) {
    data("network_emrn_step2", envir = environment())
    network_emrn <- network_emrn_step2
    rm(network_emrn_step2);gc()
  } else {
    data("network_emrn_step0", envir = environment())
    network_emrn <- network_emrn_step0
    rm(network_emrn_step0);gc()
  }

  network_emrn <- igraph::get.data.frame(network_emrn, what = 'both')

  # modify node table
  node_table <- network_emrn$vertices
  data("cpd_emrn", envir = environment())
  temp_idx <- match(node_table$name, cpd_emrn$id)
  node_table <- node_table %>%
    dplyr::left_join(cpd_emrn, by = c('name' = 'id')) %>%
    dplyr::rename(cpd_name = name.y)

  rm(cpd_emrn);gc()

  # modify edge table
  edge_table <- network_emrn$edges
  idx_from <- match(edge_table$from, node_table$name)
  idx_to <- match(edge_table$to, node_table$name)
  diff_formula <- retrieveDiffFormula(from_formula = node_table$formula[idx_from],
    to_formula = node_table$formula[idx_to],
    from_mass = node_table$monoisotopic_mass[idx_from],
    to_mass = node_table$monoisotopic_mass[idx_to])
  diff_formula <- diff_formula %>% stringr::str_replace(pattern = '\\+0', '')
  diff_formula[is.na(diff_formula)] <- 'isomer'
  edge_table <- edge_table %>% dplyr::mutate(diff_formula = diff_formula)

  # export
  path_export <- file.path(dir_path, '00_network1')
  dir.create(path_export, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(node_table, file = file.path(path_export, 'node_table.tsv'))
  readr::write_tsv(edge_table, file = file.path(path_export, 'edge_table.tsv'))

  cat('Done!\n')
}

################################################################################
# reconstructNetwork2 ----------------------------------------------------------
  # extractNeighbor ------------------------------------------------------------

#' @title extractNeighbor
#' @description extract neighbor
#' @param node_target
#' @param network_emrn
#' @param reaction_step
#' @examples

# load('./PackageFile/mrn/210329/emrn_rpair_210329.RData')
# network_emrn <- emrn_rpair$step2;rm(emrn_rpair);gc()
extractNeighbor <- function(node_target,
  network_emrn,
  reaction_step = 1) {
  subnetwork_obj <- try(igraph::make_ego_graph(network_emrn, order = reaction_step, nodes = node_target), silent = TRUE)
  if (class(subnetwork_obj) == 'try-error') return(NULL)

  cpd_names <- igraph::V(subnetwork_obj[[1]])$name
  return(cpd_names)
}



  # extractNetwork2FromNode ----------------------------------------------------

#' @title extractNetwork2FromNode
#' @author Zhiwei Zhou
#' @param id_kegg
#' @param peak
#' @param annotation_table
#' @param network_emrn
#' @param ms2
#' @param reaction_step

# id_kegg <- 'C01762'
# peak <- 'M285T555_POS'
# annotation_table <- table_id_recursive
# ms2 <- ms2_pos

# id_kegg <- 'C00082'
# peak <- 'M182T541_POS'
# annotation_table <- table_id_recursive
# ms2 <- ms2_pos

# extractNetwork2FromNode(id_kegg = 'C03618',
#                         peak = 'M130T725_POS',
#                         network_emrn = network_emrn,
#                         annotation_table = table_id_final_pos,
#                         ms2 = ms2_pos,
#                         reaction_step = 1,
#                         ms2_link = 'hybrid')


extractNetwork2FromNode <- function(id_kegg,
  peak,
  annotation_table,
  network_emrn,
  ms2,
  reaction_step = 1, # exact n_step compound
  ms2_link = c('smilarity', 'hybrid')) {
  # browser()

  ms2_link <- match.arg(ms2_link)

  # retrieve linked peaks
  # seed peak: M347T826 L0019
  # temp_subgraph <- try(extractSubgraph(id = id_kegg, step = reaction_step), silent = TRUE)
  #
  # if (class(temp_subgraph) == 'try-error') {return(NULL)}
  if (reaction_step > 1) {
    linked_cpds <- extractNeighbor(node_target = id_kegg, network_emrn = network_emrn, reaction_step = reaction_step)
    linked_cpds2 <- extractNeighbor(node_target = id_kegg, network_emrn = network_emrn, reaction_step = reaction_step-1)
    linked_cpds <- setdiff(linked_cpds, linked_cpds2)

    if (length(linked_cpds) == 0) {
      return(NULL)
    }

  } else {
    linked_cpds <- extractNeighbor(node_target = id_kegg, network_emrn = network_emrn, reaction_step = reaction_step)
    if (length(linked_cpds) == 0) {
      return(NULL)
    }

    linked_cpds <- setdiff(linked_cpds, id_kegg)
    if (length(linked_cpds) == 0) {
      return(NULL)
    }
  }

  # meta info
  result_meta <- c('peak' = peak, 'id_kegg' = id_kegg)

  # retrieve annotated linked neighbor peaks
  temp_id_peaks <- annotation_table %>%
    dplyr::filter(id_kegg %in% linked_cpds) %>%
    dplyr::filter(with_ms2) %>%
    dplyr::group_by(peak_name) %>%
    dplyr::summarise(peak_name = peak_name[1],
      mz = mz[1],
      rt = rt[1],
      id_kegg = paste(id_kegg, collapse = ';'),
      isotope = paste(isotope, collapse = ';'),
      adduct = paste(adduct, collapse = ';'),
      formula = paste(unique(formula), collapse = ';'),
      with_ms2 = paste(unique(with_ms2), collapse = ';')) %>%
    dplyr::arrange(mz)

  # if no neighbor metabolite with ms2
  if (nrow(temp_id_peaks) == 0) {
    return(NULL)
  }

  # calculate MS2 similarity centric to M267T372_NEG
  peaks_msms_similarity <- calcuPairSpecSim(annot_a = 'centric_node',
    annot_a_peak = peak,
    annot_b = 'neighbor_node',
    annot_b_peak = temp_id_peaks$peak_name,
    ms2_data = ms2)

  if (ms2_link == 'hybrid') {
    peak_ms2_linked <- peaks_msms_similarity$pair_matrix %>%
      tibble::as_tibble() %>%
      dplyr::mutate(peak_a = as.character(peak_a),
        peak_b = as.character(peak_b)) %>%
      dplyr::filter(ms2_score >= 0.5 | n_frag_match >= 4)
  } else {
    peak_ms2_linked <- peaks_msms_similarity$pair_matrix %>%
      as_tibble() %>%
      dplyr::mutate(peak_a = as.character(peak_a),
        peak_b = as.character(peak_b)) %>%
      dplyr::filter(ms2_score >= 0.5)
  }

  if (nrow(peak_ms2_linked) == 0) {
    return(NULL)
  }

  # node table
  temp_peak <- c(peak_ms2_linked$peak_a, peak_ms2_linked$peak_b) %>% unique()
  # annotation_table_included <- annotation_table %>%
  #   filter(peak_name %in% temp_peak) %>%
  #   mutate(formula = match(id_kegg, MetDNA2::cpd_emrn$id) %>% MetDNA2::cpd_emrn$formula[.])
  # annotation_table_included %>%
  #   distinct(peak_name, adduct, formula) %>% arrange(peak_name)

  node_table <- annotation_table %>%
    dplyr::filter(peak_name %in% temp_peak) %>%
    dplyr::mutate(formula = match(id_kegg, MetDNA2::cpd_emrn$id) %>% MetDNA2::cpd_emrn$formula[.]) %>%
    dplyr::group_by(peak_name) %>%
    dplyr::summarise(peak_name = peak_name[1],
      mz = mz[1],
      rt = rt[1],
      id_kegg = paste(id_kegg, collapse = ';'),
      isotope = paste(isotope, collapse = ';'),
      adduct = paste(adduct, collapse = ';'),
      formula = paste(unique(formula), collapse = ';'),
      with_ms2 = paste(unique(with_ms2), collapse = ';'))

  # edge table
  # ms2 edge
  edge_table_ms2_sim <- peak_ms2_linked %>%
    dplyr::mutate(reaction_step = reaction_step) %>%
    dplyr::select(peak_a, peak_b, ms2_score, n_frag_match, reaction_step) %>%
    dplyr::rename(from = peak_a, to = peak_b) %>%
    dplyr::mutate(edge_value = as.character(round(ms2_score, 4)),
      edge_label = 'ms2_sim')

  # rp edge
  idx_from <- match(edge_table_ms2_sim$from, node_table$peak_name)
  idx_to <- match(edge_table_ms2_sim$to, node_table$peak_name)
  edge_table_rp <- edge_table_ms2_sim %>%
    dplyr::mutate(from_formula = node_table$formula[idx_from],
      to_formula = node_table$formula[idx_to],
      from_mass = node_table$mz[idx_from],
      to_mass = node_table$mz[idx_to])

  # if one peak consisted by multiple formula, separate them
  edge_table_rp <- edge_table_rp %>%
    tidyr::separate_rows(from_formula, sep = ';') %>%
    tidyr::separate_rows(to_formula, sep = ';')

  temp_diff_formula <- retrieveDiffFormula(from_formula = edge_table_rp$from_formula,
    to_formula = edge_table_rp$to_formula,
    from_mass = edge_table_rp$from_mass,
    to_mass = edge_table_rp$to_mass)


  edge_table_rp <- edge_table_rp %>%
    dplyr::mutate(diff_formula = temp_diff_formula) %>%
    dplyr::select(from, to, ms2_score, n_frag_match, reaction_step, diff_formula) %>%
    dplyr::rename(edge_value = diff_formula) %>%
    dplyr::mutate(edge_label = 'rp')

  edge_table <- edge_table_rp %>%
    dplyr::bind_rows(edge_table_ms2_sim)


  result <- list(meta = result_meta, node_table = node_table, edge_table = edge_table)

  return(result)
}


  # reconstructNetwork2 --------------------------------------------------------

#' @title reconstructNetwork2
#' @author Zhiwei Zhou
#' @param annotation_table
#' @param dir_path
#' @param is_unknown_annotation
#' @param whether_show_single_node Default: FALSE
#' @param max_reaction_step Default: 3
#' @param ms2_link 'hybrid' or 'smilarity'; 'hybrid' means linking nodes through MS/MS similarity (>= 0.5) or matched fragments (>=4); 'smilarity' means linking nodes through MS/MS similarity (>= 0.5) only; Default: hybrid
#' @examples
#' @export


# annotation_table <- reformatTable1(dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')
#
# test <- reconstructNetwork2(annotation_table = annotation_table,
#   dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')

reconstructNetwork2 <- function(
  annotation_table,
  dir_path = '.',
  is_unknown_annotation = TRUE,
  whether_show_single_node = FALSE,
  max_reaction_step = 3,
  ms2_link = c('hybrid', 'smilarity')) {

  ms2_link <- match.arg(ms2_link)
  seed_table <- annotation_table %>% dplyr::filter(with_ms2)
  load(file.path(dir_path, 'ms2_data.RData'))
  ms2 <- raw_msms; rm(raw_msms); gc()

  if (is_unknown_annotation) {
    data("network_emrn_step2", envir = environment())
    network_emrn <- network_emrn_step2
    rm(network_emrn_step2);gc()
  } else {
    data("network_emrn_step0", envir = environment())
    network_emrn <- network_emrn_step0
    rm(network_emrn_step0);gc()
  }

  # browser()
  list_unfinished_seed <- seed_table$label_seed

  global_node_table <- NULL
  global_edge_table <- NULL
  list_finished_peak <- NULL
  list_finished_id <- NULL

  while (length(list_unfinished_seed) > 0) {
    cat('Number of unfinished seed:', length(list_unfinished_seed));cat('\n')
    idx_seed <- list_unfinished_seed[1] %>% match(seed_table$label_seed)

    new_seed_peak <-  seed_table$peak_name[idx_seed]
    new_seed_id <- seed_table$id_kegg[idx_seed]
    # cat('Seed peak:', new_seed_peak, 'Seed id:', new_seed_id);cat('\n')

    # cat('Start recursive: ')
    i <- 1
    node_table <- NULL
    edge_table <- NULL

    while (length(new_seed_peak) > 0) {
      # cat(i, ' ')

      temp_step <- 1
      # retrieve neighbor nodes via rp and ms2 similarity in annotation pool
      # if no neighbor retrieved, try step2, step3
      while (temp_step <= max_reaction_step) {
        # cat('Step', temp_step)
        result1 <- mapply(function(x, y){
          extractNetwork2FromNode(id_kegg = x,
            peak = y,
            annotation_table = annotation_table,
            network_emrn = network_emrn,
            ms2 = ms2,
            ms2_link = ms2_link,
            reaction_step = temp_step)
        },
          x = new_seed_id,
          y = new_seed_peak,
          SIMPLIFY = FALSE)

        temp <- sapply(result1, length)

        if (sum(temp) == 0) {
          temp_step <- temp_step + 1
        } else {
          break
        }
      }

      temp <- sapply(result1, length)

      if (sum(temp) == 0) {

        # whether remove node witnout any neighbor
        if (!whether_show_single_node) {
          node_table <- NULL
          edge_table <- NULL
          list_finished_peak <- c(list_finished_peak, new_seed_peak)
          list_finished_id <- c(list_finished_id, new_seed_id)
          break
        } else {
          temp_idx <- which(annotation_table$peak_name == new_seed_peak)
          node_table <- tibble::tibble(peak_name = new_seed_peak,
            mz = annotation_table$mz[temp_idx[1]],
            rt = annotation_table$rt[temp_idx[1]],
            id_kegg = annotation_table$id_kegg[temp_idx] %>% paste(collapse = ';'),
            isotope = annotation_table$isotope[temp_idx] %>% paste(collapse = ';'),
            adduct = annotation_table$adduct[temp_idx] %>% paste(collapse = ';'),
            formula = annotation_table$formula[temp_idx] %>% unique() %>% paste(collapse = ';'),
            with_ms2 = annotation_table$with_ms2[temp_idx[1]])

          edge_table <- tibble::tibble(from = c(new_seed_peak, new_seed_peak),
            to = c(new_seed_peak, new_seed_peak),
            ms2_score = c(1, 1),
            n_frag_match = c(0, 0),
            reaction_step = c(0, 0),
            edge_value = c('', ''),
            edge_label = c('rp', 'ms2_sim'))

          list_finished_peak <- c(list_finished_peak, new_seed_peak)
          list_finished_id <- c(list_finished_id, new_seed_id)
          break
        }
      }

      temp_node_table <- lapply(result1, function(x){x$node_table}) %>%
        dplyr::bind_rows()

      temp_edge_table <- lapply(result1, function(x){x$edge_table}) %>%
        dplyr::bind_rows()

      node_table <- node_table %>%
        dplyr::bind_rows(temp_node_table) %>%
        dplyr::distinct(peak_name, .keep_all = TRUE)

      edge_table <- edge_table %>%
        dplyr::bind_rows(temp_edge_table) %>%
        dplyr::mutate(temp1 = pmin(from, to), temp2 = pmax(from, to)) %>%
        dplyr::distinct(temp1, temp2, edge_label, edge_value, .keep_all = TRUE) %>%
        dplyr::filter(from != to) %>%
        dplyr::select(-c(temp1, temp2))

      list_finished_peak <- c(list_finished_peak, new_seed_peak)
      list_finished_id <- c(list_finished_id, new_seed_id)

      temp <- annotation_table %>%
        dplyr::filter(peak_name %in% node_table$peak_name) %>%
        dplyr::filter(!(peak_name %in% list_finished_peak) & !(id_kegg %in% list_finished_id)) %>%
        dplyr::distinct(peak_name, id_kegg, .keep_all = TRUE)

      new_seed_peak <- temp$peak_name
      new_seed_id <- temp$id_kegg

      i <- i + 1

    }

    temp_finished <- paste(list_finished_peak, list_finished_id, sep = '_')

    global_node_table <- global_node_table %>% dplyr::bind_rows(node_table)
    global_edge_table <- global_edge_table %>% dplyr::bind_rows(edge_table)

    # select new seed not used as seed
    list_unfinished_seed <- setdiff(list_unfinished_seed, temp_finished)

    if (length(list_unfinished_seed) == 0) break
    # cat('\n')
  }

  # save result
  # remove redundancy
  global_node_table <- global_node_table %>%
    dplyr::distinct(peak_name, .keep_all = TRUE)

  global_edge_table <- global_edge_table %>%
    dplyr::mutate(temp1 = pmin(from, to), temp2 = pmax(from, to)) %>%
    dplyr::distinct(temp1, temp2, edge_label, edge_value, .keep_all = TRUE) %>%
    dplyr::select(-c(temp1, temp2))

  # browser()

  # add peak label
  # if any candidates have unknown annotation, this peak is labeled as unknown_peak
  peak_label <- sapply(unique(annotation_table$peak_name), function(x){
    temp <- annotation_table %>% dplyr::filter(peak_name == x) %>%
      dplyr::pull(confidence_level)
    if (any(temp %in% c('level1', 'level2'))) return('seed_peak')
    if (any(temp %in% c('level3.2'))) return('unknown_peak')
    return('known_peak')
  })

  global_node_table$label <- match(global_node_table$peak_name, names(peak_label)) %>%
    peak_label[.]

  network2_result <- list(node_table = global_node_table,
    edge_table = global_edge_table)

  dir.create(file.path(dir_path, '01_files_network2'),
    showWarnings = FALSE, recursive = TRUE)

  readr::write_tsv(network2_result$node_table,
    file = file.path(dir_path, '01_files_network2', 'node_table_network2.tsv'))

  readr::write_tsv(network2_result$edge_table,
    file = file.path(dir_path, '01_files_network2', 'edge_table_network2.tsv'))

  save(network2_result,
    file = file.path(dir_path, '01_files_network2', 'network2_result.RData'))

  # export network2 igraph object
  network2_result$node_table <- network2_result$node_table %>% dplyr::rename(name = peak_name)
  network2_obj <- igraph::graph.data.frame(d = network2_result$edge_table,
    vertices = network2_result$node_table)

  save(network2_obj,
    file = file.path(dir_path, '01_files_network2', 'network2_obj.RData'))
}


################################################################################
# reconstructNetwork3 --------------------------------------------------------

#' @title reconstructNetwork3
#' @author Zhiwei Zhou
#' @description retrieve network3 for visualization
#' @param dir_path path of working directory
#' @return a list of network3, includes 2 objects (node_table & edge_table)
#' @importFrom magrittr '%>%' '%$%'
#' @export
#' @examples
#' setwd('D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')
#' reconstructNetwork3()


# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/03_annotation_credential/00_intermediate_data/list_peak_group_annotation_concised.RData')
# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/00_annotation_table/00_intermediate_data/list_peak_group')
# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/00_annotation_table/00_intermediate_data/peak_group_id_table')
# load('./Data/20220118_biological_sample_analysis/nist_urine_pos/00_annotation_table/00_intermediate_data/table_identification')
#
# reconstructNetwork3(dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')

reconstructNetwork3 <- function(dir_path = '.') {
  cat('Start reconstructing Network3 for visualization...\n')
  # browser()
  load(file.path(dir_path, 'list_peak_group'))
  load(file.path(dir_path, 'peak_group_id_table'))
  load(file.path(dir_path, 'table_identification'))
  load(file.path(dir_path, 'list_peak_group_annotation_concised.RData'))

  peak_group_name <- match(names(list_peak_group), peak_group_id_table$peak_group_id) %>%
    peak_group_id_table$included_peak_group[.]

  # node table
  node_table <- list_peak_group %>%
    dplyr::bind_rows() %>%
    dplyr::select(peak_name:rt, adduct, isotope, peak_group_id:rela_peaks)

  # assign node attribute: seed, unknown, abiotic peaks
  temp_ids <- table_identification %>%
    dplyr::group_by(peak_name) %>%
    dplyr::summarise(peak_name = peak_name[1],
      mz = mz[1],
      rt = rt[1],
      id_kegg = paste(id_kegg, collapse = ';'),
      id_zhulab = paste(id_zhulab, collapse = ';'),
      confidence_level = paste(unique(confidence_level), collapse = ';'),
      isotope = paste(isotope, collapse = ';'),
      adduct = paste(adduct, collapse = ';'),
      formula = paste(unique(formula), collapse = ';'),
      peak_group_id = paste(unique(peak_group_id), collapse = ';'),
      num_peaks = paste(unique(num_peaks), collapse = ';'))

  # add label
  idx_seed_peak <- temp_ids %>%
    dplyr::filter(stringr::str_detect(confidence_level, 'level1|level2')) %>%
    dplyr::pull(peak_name) %>%
    match(node_table$peak_name)

  idx_known_peak <- temp_ids %>%
    dplyr::filter(confidence_level == 'level3.1') %>%
    dplyr::pull(peak_name) %>%
    match(node_table$peak_name)

  idx_unknown_peak <- temp_ids %>%
    dplyr::filter(stringr::str_detect(confidence_level, 'level3.2')) %>%
    dplyr::pull(peak_name) %>%
    match(node_table$peak_name)

  node_table$label <- 'abiotic_peak'
  node_table$label[idx_seed_peak] <- 'seed_peak'
  node_table$label[idx_known_peak] <- 'known_peak'
  node_table$label[idx_unknown_peak] <- 'unknown_peak'

  # edge table
  list_pg_final <- match(peak_group_name, names(list_peak_group_annotation_concised)) %>%
    list_peak_group_annotation_concised[.]

  edge_table <- lapply(list_pg_final, function(x){
    x@peak_list_annotated %>%
      dplyr::mutate(base_peak = paste(x@base_peak_name, x@base_peak_adduct, sep = '_'))
  }) %>% dplyr::bind_rows()

  edge_label <- apply(edge_table[, 13:16], 1, function(x){
    # browser()
    if (x[1] != '[M]') {
      return(x[1])
    } else {
      idx <- which(!(is.na(x[2:4])))
      return(x[idx+1])
    }
  })

  edge_table <- edge_table %>%
    dplyr::select(query_peak, peak_name, dplyr::everything()) %>%
    dplyr::rename(from = query_peak, to = peak_name) %>%
    dplyr::mutate(ms2_score = -1,
      n_frag_match = -1,
      reaction_step = -1,
      edge_value = edge_label,
      edge_label = 'peak_correlation') %>%
    dplyr::select(from, to, ms2_score, n_frag_match, reaction_step, base_peak, edge_value, edge_label)

  network3_result <- list(node_table = node_table, edge_table = edge_table)

  dir.create(file.path(dir_path, '02_files_network3'), showWarnings = FALSE, recursive = TRUE)

  readr::write_tsv(network3_result$node_table,
    file = file.path(dir_path, '02_files_network3', 'node_table_network3.tsv'))

  readr::write_tsv(network3_result$edge_table,
    file = file.path(dir_path, '02_files_network3', 'edge_table_network3.tsv'))

  save(network3_result,
    file = file.path(dir_path, '02_files_network3', 'network3_result.RData'))

  # export network3 igraph object
  network3_result$node_table <- network3_result$node_table %>%
    dplyr::rename(name = peak_name) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(name = name[1],
      mz = mz[1],
      rt = rt[1],
      adduct = paste(adduct, collapse = ';'),
      isotope = paste(isotope, collapse = ';'),
      peak_group_id = paste(peak_group_id, collapse = ';'),
      base_peak = paste(base_peak, collapse = ';'),
      rela_peaks = paste(rela_peaks, collapse = ';'),
      label = paste(label, collapse = ';'))
  network3_obj <- igraph::graph.data.frame(d = network3_result$edge_table,
    vertices = network3_result$node_table)

  save(network3_obj,
    file = file.path(dir_path, '02_files_network3', 'network3_obj.RData'))

  cat('Done!\n')
}







################################################################################
# retrieve subnetworks ------------------------------------------------------------
  # retrieveSubNetwork1 --------------------------------------------------------

#' @title retrieveSubNetwork1
#' @description retrieve subnetwork1 from constructed metabolic reaction network
#' @param centric_met charecter vector. If the length of centric_met equals 1, then retrieve subnetwork with defined step. If the length of centric_met >= 2, then only retrieve subnetwork with nodes from input centric_met.
#' @param step neigbor step. Default: 1
#' @param is_unknown_annotation whether used unknown annotation. Default: TRUE
#' @param dir_path path of working directory. Default: '.'
#' @param show_plot whether show subnetwork plot in R. Default: TRUE
#' @export
#' @examples
#' # single node input
#' retrieveSubNetwork1()

# retrieveSubNetwork1(centric_met = 'C00082',
#   dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/')
#
# retrieveSubNetwork1(centric_met = c('C00082', 'KeggExd000923'),
#   dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/')


  # centric_met <- 'C00082'
retrieveSubNetwork1 <- function(centric_met,
    step = 1,
    is_unknown_annotation = TRUE,
    dir_path = '.',
    show_plot = TRUE) {

    if (missing(centric_met)) {
      stop('Please input centric_met\n')
    }

    if (is_unknown_annotation) {
      data("network_emrn_step2", envir = environment())
      network_emrn <- network_emrn_step2
      rm(network_emrn_step2);gc()
    } else {
      data("network_emrn_step0", envir = environment())
      network_emrn <- network_emrn_step0
      rm(network_emrn_step0);gc()
    }

    if (!all(centric_met %in% names(igraph::V(network_emrn)))) {
      error <- which(centric_met %in% names(igraph::V(network_emrn))) %>%
        centric_met[.]
      stop('Nodes: ', paste(error, collapse = ';'), 'are not included in KMRN!\n')
    }

    if (length(centric_met) == 1) {
      sub_graph <- igraph::make_ego_graph(network_emrn,
        order = step,
        nodes = centric_met)
      sub_graph <- sub_graph[[1]]
    } else {
      sub_graph <- igraph::subgraph(network_emrn, v = c(centric_met))
    }


    # extract node table and edge table
    sub_graph_result <- sub_graph %>% igraph::get.data.frame(what = 'both')

    # modify node table
    node_table <- sub_graph_result$vertices
    data("cpd_emrn", envir = environment())
    temp_idx <- match(node_table$name, cpd_emrn$id)
    node_table <- node_table %>%
      dplyr::left_join(cpd_emrn, by = c('name' = 'id')) %>%
      dplyr::rename(cpd_name = name.y)
    rm(cpd_emrn);gc()

    # modify edge table
    edge_table <- sub_graph_result$edges
    idx_from <- match(edge_table$from, node_table$name)
    idx_to <- match(edge_table$to, node_table$name)
    diff_formula <- retrieveDiffFormula(from_formula = node_table$formula[idx_from],
      to_formula = node_table$formula[idx_to],
      from_mass = node_table$monoisotopic_mass[idx_from],
      to_mass = node_table$monoisotopic_mass[idx_to])
    diff_formula <- diff_formula %>% stringr::str_replace(pattern = '\\+0', '')
    diff_formula[is.na(diff_formula)] <- 'isomer'
    edge_table <- edge_table %>% dplyr::mutate(diff_formula = diff_formula)

    path_export <- file.path(dir_path,
      '03_subnetworks',
      paste(centric_met, collapse = '_'),
      'network1')

    dir.create(path_export, showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(node_table, file = file.path(path_export, 'node_table.tsv'))
    readr::write_tsv(edge_table, file = file.path(path_export, 'edge_table.tsv'))

    sub_graph <- igraph::graph.data.frame(d = edge_table,
      vertices = node_table,
      directed = FALSE)
    save(sub_graph, file = file.path(path_export, 'sub_graph.RData'))

    subnetwork1_result <- list(node_table = node_table, edge_table = edge_table)
    save(subnetwork1_result, file = file.path(path_export, 'subnetwork1_result.RData'))

    # plot example
    if (show_plot) {
      subgraph_tbl <- tidygraph::tbl_graph(nodes = node_table,
        edges = edge_table)

      temp_plot <- ggraph::ggraph(subgraph_tbl, layout = 'nicely') +
        ggraph::geom_edge_fan(ggplot2::aes(label = diff_formula), colour = 'black') +
        ggraph::geom_node_point(ggplot2::aes(colour = type, shape = type), size = 5) +
        ggplot2::scale_colour_manual(values = c('known_known' = 'dodgerblue',
          'known_unknown' = 'tomato',
          'unknown_unknown' = 'tomato')) +
        ggplot2::scale_shape_manual(values = c('known_known' = 16,
          'known_unknown' = 17,
          'unknown_unknown' = 17)) +
        ggraph::geom_node_text(ggplot2::aes(label = name)) +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = 'top')

      return(temp_plot)
    }
  }

  # retrieveSubNetwork2 ----------------------------------------------------

  #' @title retrieveSubNetwork2
  #' @author Zhiwei Zhou
  #' @param from_peak
  #' @param end_peak Character. Default: NULL
  #' @param step Numeric. Default: step = 1
  #' @param dir_path path of working directory. Default: '.'
  #' @export
  #' @examples
  #' # subnetwork from centric peak
  #'   retrieveSubNetwork2(from_peak = 'M196T420', step = 1, dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')
  #'
  #' # subnetwork with defined from_peak and end_peak
  #'    retrieveSubNetwork2(from_peak = 'M262T526', end_peak = 'M182T541', dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')
  #'

  # retrieveSubNetwork2(from_peak = 'M196T420',
  #   step = 1,
  #   dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')
  #
  # retrieveSubNetwork2(from_peak = 'M262T526',
  #   end_peak = 'M182T541',
  #   dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')

  # load('D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/01_files_network2/network2_result.RData')
  # network2_result$node_table <- network2_result$node_table %>% dplyr::rename(name = peak_name)
  # dir_path <- 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization'
  # from_peak <- 'M262T526'
  # end_peak <- 'M182T541'

  retrieveSubNetwork2 <- function(from_peak,
    end_peak = NULL,
    step = 1,
    dir_path = '.',
    show_plot = TRUE) {
    # browser()
    if (missing(from_peak)) {
      stop('Please input from_peak')
    }

    if (length(from_peak) > 1 | length(end_peak) > 1) {
      stop('The input from_peak or to_peak > 1\n')
    }

    if (!('network2_obj.RData' %in% list.files(file.path(dir_path, '01_files_network2')))) {
      stop('No network2_result existed!')
    } else {
      load(file.path(dir_path, '01_files_network2', 'network2_obj.RData'))
    }

    temp_node <- c(from_peak, end_peak)

    if (!all(temp_node %in% names(igraph::V(network2_obj)))) {
      error <- which(temp_node %in% names(igraph::V(network2_obj))) %>% temp_node[.]
      stop('Nodes: ', paste(error, collapse = ';'), 'are not included in network2!\n')
    }

    if (length(end_peak) == 0) {
      sub_graph <- igraph::make_ego_graph(network2_obj,
        order = step,
        nodes = from_peak)
      sub_graph <- sub_graph[[1]]
    } else {
      sub_graph <- igraph::subgraph(network2_obj, v = c(from_peak, end_peak))
    }

    # export subgraph examples
    ifelse(length(end_peak) > 0,
      path_export <- file.path(dir_path,
        '03_subnetworks',
        paste(from_peak, end_peak, sep = '_'),
        'network2'),
      path_export <- file.path(dir_path,
        '03_subnetworks',
        from_peak,
        'network2'))

    dir.create(path_export,
      showWarnings = FALSE, recursive = TRUE)

    sub_graph_result <- igraph::get.data.frame(sub_graph, what = 'both')
    node_table <- sub_graph_result$vertices
    readr::write_tsv(node_table,
      file = file.path(path_export, 'node_table.tsv'))

    edge_table <- sub_graph_result$edges
    readr::write_tsv(edge_table, file.path(path_export, 'edge_table.tsv'))

    subnetwork2_result <- list(node_table = node_table,
      edge_table = edge_table)
    save(subnetwork2_result, file = file.path(path_export, 'subnetwork2_result.RData'))

    save(sub_graph, file = file.path(path_export, 'sub_graph.RData'))

    # plot example
    if (show_plot) {
      subgraph_tbl <- tidygraph::tbl_graph(nodes = node_table,
        edges = edge_table)

      temp_plot <- ggraph::ggraph(subgraph_tbl) +
        ggraph::geom_edge_fan(ggplot2::aes(colour = edge_label)) +
        ggraph::geom_node_point(ggplot2::aes(colour = label, shape = label), size = 5) +
        ggplot2::scale_colour_manual(values = c('seed_peak' = 'orange',
          'known_peak' = 'dodgerblue',
          'unknown_peak' = 'tomato')) +
        ggplot2::scale_shape_manual(values = c('seed_peak' = 16,
          'known_peak' = 16,
          'unknown_peak' = 17)) +
        ggraph::scale_edge_color_manual(values = c('ms2_sim' = 'tomato',
          'rp' = 'black')) +
        ggraph::geom_node_text(ggplot2::aes(label = name),
          # angle = 45,
          repel = TRUE) +
        ggplot2::theme_void()

      return(temp_plot)
    }

  }

  # retrieveSubNetwork3 --------------------------------------------------------

#' @title retrieveSubNetwork3
#' @author Zhiwei Zhou
#' @param base_peaks a vector of base_peaks
#' @param base_adducts a vector of base_adducts. Note: it should keep same length with base_peaks
#' @param dir_path path of working directory.
#' @param show_plot Whether show subnetwork plot. Default: TRUE
#' @export
#' @examples
#' retrieveSubNetwork3(base_peaks = c('M182T541', 'M262T526'),
#' base_adducts = c('[M+H]+', '[M+H]+'))

  # load('D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/02_files_network3/network3_result.RData')
  # base_peaks <- c('M182T541', 'M262T526')
  # base_adducts <- c('[M+H]+', '[M+H]+')

# retrieveSubNetwork3(base_peaks = c('M262T526', 'M182T541'),
#   base_adducts = c('[M+H]+', '[M+H]+'),
#   dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization')

  retrieveSubNetwork3 <- function(base_peaks,
    base_adducts,
    dir_path = '.',
    show_plot = TRUE) {
    # browser()
    if (missing(base_peaks)) {
      stop('Please input base_peaks\n')
    }

    if (missing(base_adducts)) {
      stop('Please input base_adducts\n')
    }

    if (length(base_peaks) != length(base_adducts)) {
      stop('The length of base_peaks should equal with the length of base_adducts\n')
    }

    if (!('network3_result.RData' %in% list.files(file.path(dir_path, '02_files_network3')))) {
      stop('No network3_result existed!\n')
    } else {
      load(file.path(dir_path, '02_files_network3', 'network3_result.RData'))
    }

    base_peak_group <- paste(base_peaks, base_adducts, sep = '_')
    temp_peak_group_list <- unique(network3_result$node_table$base_peak)

    if (!all(base_peak_group %in% temp_peak_group_list)) {
      error <- which(base_peak_group %in% temp_peak_group_list) %>% temp_peak_group[.]
      stop('Nodes: ', paste(error, collapse = ';'), 'are not included in network2!\n')
    }
    rm(list = c('temp_peak_group_list'));gc()

    # extract node_table
    node_table <- network3_result$node_table %>%
      dplyr::filter(base_peak %in% base_peak_group)

    # if redundancy exist, then merge them
    redun_peak <- table(node_table$peak_name)
    redun_peak <- names(redun_peak)[redun_peak>1]
    if (length(redun_peak) > 0) {
      label_redun_peak <- sapply(redun_peak, function(x){
        temp_idx <- which(node_table$peak_name == x)
        temp_label <- node_table$label[temp_idx]
        temp_label <- match(temp_label,
          c('seed_peak', 'known_peak', 'unknown_peak', 'abiotic_peak')) %>%
          which.min() %>%
          temp_label[.]
      })
    }

    node_table <- node_table %>%
      dplyr::rename(name = peak_name) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(name = name[1],
        mz = mz[1],
        rt = rt[1],
        adduct = paste(adduct, collapse = ';'),
        isotope = paste(isotope, collapse = ';'),
        peak_group_id = paste(peak_group_id, collapse = ';'),
        base_peak = paste(base_peak, collapse = ';'),
        rela_peaks = paste(rela_peaks, collapse = ';'),
        label = paste(label, collapse = ';')) %>%
      dplyr::ungroup()

    if (length(redun_peak) > 0) {
      node_table$label[match(redun_peak, node_table$name)] <- label_redun_peak
    }


    # # replace_as base peaks
    # node_table$label[node_table$label != 'abiotic_peak'] <- 'base_peak'


    # extract edge_table
    edge_table <- network3_result$edge_table %>%
      dplyr::filter(base_peak %in% base_peak_group)
    # temp_idx <- which(edge_table$from == edge_table$to)
    # edge_table$edge_label[temp_idx] <- 'base_peak'


    # export results
    if (length(base_peaks) <= 2) {
      path_export <- file.path(dir_path,
        '03_subnetworks',
        paste(base_peaks, collapse = '_'),
        'network3')
    } else {
      path_export <- file.path(dir_path,
        '03_subnetworks',
        paste(base_peaks[1:2], collapse = '_'),
        'network3')
    }


    dir.create(path_export, showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(node_table,
      file = file.path(path_export, 'node_table.tsv'))
    readr::write_tsv(edge_table,
      file = file.path(path_export, 'edge_table.tsv'))

    subnetwork3_result <- list(node_table = node_table,
      edge_table = edge_table)
    save(subnetwork3_result, file = file.path(path_export, 'subnetwork3_result.RData'))


    sub_graph <- igraph::graph.data.frame(d = edge_table, vertices = node_table)
    save(sub_graph, file = file.path(path_export, 'sub_graph.RData'))

    # plot example
    if (show_plot) {
      # browser()
      subgraph_tbl <- tidygraph::tbl_graph(nodes = node_table,
        edges = edge_table)

      # options(ggrepel.max.overlaps = Inf)
      temp_plot <- ggraph::ggraph(subgraph_tbl) +
        ggraph::geom_edge_fan(ggplot2::aes(label = edge_value), colour = 'limegreen') +
        ggraph::geom_node_point(ggplot2::aes(colour = label, shape = label),
          size = 5) +
        ggplot2::scale_colour_manual(values = c('seed_peak' = 'orange',
          'known_peak' = 'dodgerblue',
          'unknown_peak' = 'tomato',
          'abiotic_peak' = 'limegreen')) +
        ggplot2::scale_shape_manual(values = c('seed_peak' = 16,
          'known_peak' = 16,
          'unknown_peak' = 17,
          'abiotic_peak' = 16)) +
        # ggraph::scale_edge_color_manual(values = c('peak_correlation' = 'limegreen')) +
        ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
        ggplot2::theme_void()

      return(temp_plot)
    }

  }


################################################################################
# mergeSubNetworks -------------------------------------------------------------

#' @title mergeSubnetwork
#' @description merge subnetwork2 and subnetwork3 from KGMN. Note: this function only effective after performing subnetwork2 and subnetwork3 extraction
#' @author Zhiwei Zhou
#' @param from_peak
#' @param end_peak
#' @param dir_path path of working directory. Default: '.'
#' @param show_plot whether show merged network. Default: TRUE
#' @export
#' @examples
#' mergeSubnetwork(from_peak = 'M262T526', end_peak = 'M182T541', dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/')

# dir_path <- 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/'
# from_peak <- 'M262T526'
# end_peak <- 'M182T541'
# mergeSubnetwork(from_peak = 'M262T526', end_peak = 'M182T541', dir_path = 'D:/project/00_zhulab/01_metdna2/00_data/20220602_visualization_kgmn/06_visualization/')

mergeSubnetwork <- function(from_peak,
  end_peak,
  dir_path = '.',
  show_plot = TRUE) {

  if (missing(from_peak)) {
    stop('Please input from_peak')
  }

  if (missing(end_peak)) {
    stop('Please input end_peak')
  }

  if (!(paste(from_peak, end_peak, sep = '_') %in%
      list.files(file.path(dir_path, '03_subnetworks')))) {
    stop('Please run subnetwork extraction first for network2 and network3\n')
  }

  temp_folder <- file.path(dir_path,
    '03_subnetworks',
    paste(from_peak, end_peak, sep = '_')) %>%
    list.files()
  if (!(all(c('network2', 'network3') %in% temp_folder))) {
    error <- which(!(c('network2', 'network3') %in% temp_folder))
    stop('Nodes: ', paste(error, collapse = ';'), 'are not retrieved!\n')
  }

  # load subnetworks of each layer network
  load(file.path(dir_path,
    '03_subnetworks',
    paste(from_peak, end_peak, sep = '_'),
    'network2',
    'subnetwork2_result.RData'))

  load(file.path(dir_path,
    '03_subnetworks',
    paste(from_peak, end_peak, sep = '_'),
    'network3',
    'subnetwork3_result.RData'))

  # merge subnetworks
  node_table_subnetwork2 <- subnetwork2_result$node_table %>%
    dplyr::select(name:rt, id_kegg, isotope, adduct, formula, label)
  node_table_subnetwork3 <- subnetwork3_result$node_table %>%
    dplyr::mutate(id_kegg = '', formula = '') %>%
    dplyr::select(name:rt, id_kegg, isotope, adduct, formula, label)
  node_table <- node_table_subnetwork2 %>%
    dplyr::bind_rows(node_table_subnetwork3) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  rownames(node_table) <- NULL

  edge_table_subnetwork2 <- subnetwork2_result$edge_table
  edge_table_subnetwork3 <- subnetwork3_result$edge_table %>%
    dplyr::select(-base_peak)

  edge_table <- edge_table_subnetwork2 %>%
    dplyr::bind_rows(edge_table_subnetwork3)

  # export merged networks
  path_export <- file.path(dir_path,
    '03_subnetworks',
    paste(from_peak, end_peak, sep = '_'),
    'network_merge')

  dir.create(path_export, showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(node_table,
    file = file.path(path_export, 'node_table.tsv'))
  readr::write_tsv(edge_table,
    file = file.path(path_export, 'edge_table.tsv'))

  network_merge_result <- list(node_table = node_table, edge_table = edge_table)
  save(network_merge_result,
    file = file.path(path_export, 'network_merge_result.RData'))

  sub_graph <- igraph::graph.data.frame(d = edge_table, vertices = node_table)
  save(sub_graph, file = file.path(path_export, 'sub_graph.RData'))


  # plot example
  if (show_plot) {
    # browser()
    subgraph_tbl <- tidygraph::tbl_graph(nodes = node_table,
      edges = edge_table)

    # options(ggrepel.max.overlaps = Inf)
    temp_plot <- ggraph::ggraph(subgraph_tbl) +
      ggraph::geom_edge_fan(ggplot2::aes(label = edge_value, colour = edge_label)) +
      ggraph::geom_node_point(ggplot2::aes(colour = label, shape = label),
        size = 5) +
      ggplot2::scale_colour_manual(values = c('seed_peak' = 'orange',
        'known_peak' = 'dodgerblue',
        'unknown_peak' = 'tomato',
        'abiotic_peak' = 'limegreen')) +
      ggplot2::scale_shape_manual(values = c('seed_peak' = 16,
        'known_peak' = 16,
        'unknown_peak' = 17,
        'abiotic_peak' = 16)) +
      ggraph::scale_edge_color_manual(values = c('ms2_sim' = 'tomato',
        'rp' = 'black',
        'peak_correlation' = 'limegreen')) +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE) +
      ggplot2::theme_void()

    return(temp_plot)
  }

}



################################################################################
# small tools ------------------------------------------------------------------
  # getCompoundInfo ------------------------------------------------------------

#' @title getCompoundInfo
#' @description query compound information in KGMN
#' @author Zhiwei Zhou
#' @param query_id The compound id in annotation_table one. Single compound query is supported
#' @export
#' @examples
#' # compounds in KMRN
#' getCompoundInfo('C01802')
#' getCompoundInfo('KeggExd000923')
#'
#' # compounds in library
#' getCompoundInfo('L0449')

# getCompoundInfo('C01802')
# getCompoundInfo('L0449')

getCompoundInfo <- function(query_id) {
  if (length(query_id) > 1) {
    stop('This function is only support query one compound per time\n')
  }

  data("cpd_emrn", envir = environment())
  data("cpd_lib", envir = environment())

  temp_ids <- c(cpd_emrn$id, cpd_lib$id, cpd_lib$id_kegg) %>% unique()
  if (!(query_id %in% temp_ids)) {
    stop('Queryed compound not included in KGMN\n')
  }

  if (query_id %in% cpd_emrn$id) {
    cpd_emrn %>%
      dplyr::filter(id == query_id) %>%
      showDfLine()
  } else if (query_id %in% cpd_lib$id) {
    cpd_lib %>%
      dplyr::filter(id == query_id) %>%
      showDfLine()
  } else {
    cpd_lib %>%
      dplyr::filter(id_kegg == query_id) %>%
      showDfLine()
  }

}


################################################################################


################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
If you have any questions, please send email to zhouzw@sioc.ac.cn or jiangzhu@sioc.ac.cn.
Authors: Dr. Zhiwei Zhou, and Dr. Zhengjiang Zhu (jiangzhu@sioc.ac.cn).
Maintainer: Zhiwei Zhou

Version 0.1.0 (20220604)
-------------
o Initial version for visualize multi-layer network results form KGMN
")
}

