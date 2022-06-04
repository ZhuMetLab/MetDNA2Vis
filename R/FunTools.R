# convertSpectraData ---------------------------------------------------------
#' @title convertSpectraData
#' @param ms2_data
#' @importClassesFrom SpectraTools 'SpectraData'
#' @export

convertSpectraData <- function(ms2_data) {
  options(readr.num_columns = 0)
  temp_info <- ms2_data$info %>%
    dplyr::rename(name = NAME,
      mz = PRECURSORMZ) %>%
    dplyr::select(name:mz) %>%
    readr::type_convert()

  temp_ms2_data <- ms2_data$spec

  result <- new('SpectraData',
    info = temp_info,
    spectra = list(temp_ms2_data))

  return(result)
}







# runSpecMatch -----------------------------------------------------------------

#' @title runSpecMatch
#' @description a interphace of runing SpectraTools
#' @author Zhiwei Zhou
#' @param obj_ms2_cpd1 experimental ms2 object. Note: info is a data.frame, spec is a list formed by matrix
#' @param obj_ms2_cpd2 library ms2 object.
#' @param mz_tol_ms2 Default: 35 ppm
#' @param scoring_approach 'dp', 'bonanza', 'hybrid', 'gnps'

# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd1_for_metdna2.RData')
# load('/home/zhouzw/Data_processing/20210224_metdna2_development_test/obj_ms2_cpd2_for_metdna2.RData')
# score_dp <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'dp')
# score_bonanza <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'bonanza')
# score_hybrid <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'hybrid')
# score_gnps <- runSpecMatch(obj_ms2_cpd1 = obj_ms2_cpd1, obj_ms2_cpd2 = obj_ms2_cpd2, mz_tol_ms2 = 35, scoring_approach = 'gnps')

runSpecMatch <- function(
    obj_ms2_cpd1,
  obj_ms2_cpd2,
  mz_tol_ms2 = 35,
  scoring_approach = c('dp', 'bonanza', 'hybrid', 'gnps'),
  ...
) {
  # browser()
  match.arg(scoring_approach)

  switch (scoring_approach,
    'dp' = {
      intensityNormedMethod <- 'maximum'
      methodScore <- 'dp'
    },
    'bonanza' = {
      intensityNormedMethod <- 'bonanza'
      methodScore <- 'bonanza'
    },
    'hybrid' = {
      intensityNormedMethod <- 'maximum'
      methodScore <- 'hybrid'
    },
    'gnps' = {
      intensityNormedMethod <- 'gnps'
      methodScore <- 'gnps'
    }
  )

  matchParam <- SpectraTools::MatchParam(ppm = mz_tol_ms2,
    cutoff = 0,
    weightIntensity = 1,
    weightMZ = 0,
    normIntensity = TRUE,
    tuneLibSpectra = TRUE,
    intensityExpNormed = TRUE,
    intensityLibNormed = TRUE,
    includePrecursor = TRUE,
    ppmPrecursorFilter = 30,
    thrIntensityAbs = 0,
    thrIntensityRel = 0,
    intensityNormedMethod = intensityNormedMethod,
    methodMatch = 'direct',
    methodScore = methodScore) %>%
    new(Class = 'MatchParam')


  result <- try(SpectraTools::MatchSpectra(dataExp = obj_ms2_cpd1,
    dataRef = obj_ms2_cpd2,
    matchParam),
    silent = TRUE)




  # # if the spectra has only one fragment, and it larger than precursor, it was removed
  # if (length(result) == 0) {
  #   n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
  #   n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[1]])
  #   n_frag_match <- 0
  #   n_nl_match <- 0
  #
  #   if (scoring_approach == 'dp') {
  #
  #   }
  #   result@info <- obj_ms2_cpd2@info %>%
  #     dplyr::mutate(scoreReverse = 0,
  #                   scoreForward = 0) %>%
  #     dplyr::mutate(n_frag_cpd1 = n_frag_cpd1,
  #                   n_frag_cpd2 = n_frag_cpd2,
  #                   n_frag_match = n_frag_match,
  #                   n_frag_nl = n_nl_match) %>%
  #     dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)
  # }

  # add matched_frag and matched_nl into the result table
  stat_matched_frag <- lapply(seq_along(result@matchedFragments), function(i){
    temp_matchedFragments <- result@matchedFragments[[i]]

    if (length(temp_matchedFragments) > 0) {
      n_frag_cpd1 <- temp_matchedFragments %>%
        dplyr::filter(intensity > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

      n_frag_cpd2 <- temp_matchedFragments %>%
        dplyr::filter(intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

      n_frag_match <- temp_matchedFragments %>%
        dplyr::filter(intensity > 0 & intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()

    } else {
      n_frag_cpd1 <- nrow(obj_ms2_cpd1@spectra[[1]])
      n_frag_cpd2 <- nrow(obj_ms2_cpd2@spectra[[i]])
      n_frag_match <- 0
    }

    if (scoring_approach == 'dp') {
      n_nl_match <- 0
    } else {
      n_nl_match <- result@nlFragments[[1]] %>%
        dplyr::filter(intensity > 0 & intensityExp > 0) %>%
        dplyr::count() %>%
        dplyr::pull()
    }

    temp_result <- tibble::tibble(n_frag_cpd1 = n_frag_cpd1,
      n_frag_cpd2 = n_frag_cpd2,
      n_frag_match = n_frag_match,
      n_frag_nl = n_nl_match) %>%
      dplyr::mutate(n_frag_total = n_frag_match + n_frag_nl)

  })

  stat_matched_frag <- stat_matched_frag %>% dplyr::bind_rows()

  result@info <- result@info %>%
    dplyr::bind_cols(stat_matched_frag)

  return(result)
}




# calcuPairSpecSim -------------------------------------------------------------

#' @title calcuPairSpecSim
#' @author Zhiwei Zhou
#' @param calculate
calcuPairSpecSim <- function(annot_a = 'C00204',
  annot_a_peak = c('M177T553'),
  annot_b = 'C04349',
  annot_b_peak = c('M235T596', 'M157T548'),
  ms2_data = raw_msms) {
  # browser()

  if (length(annot_a) > 1) {stop('Only one metabolite is required for a')}
  if (length(annot_b) > 1) {stop('Only one metabolite is required for b')}

  pair_matrix <- expand.grid(annot_a_peak, annot_b_peak) %>%
    dplyr::rename(peak_a = Var1,
      peak_b = Var2) %>%
    dplyr::mutate(annot_a = annot_a,
      annot_b = annot_b) %>%
    dplyr::select(annot_a, annot_b, peak_a, peak_b)

  idx_a <- match(as.character(pair_matrix$peak_a), names(ms2_data))
  if (length(idx_a) == 0) {stop('a_peak has no ms2!')}
  ms2_a <- ms2_data[idx_a]

  idx_b <- match(as.character(pair_matrix$peak_b), names(ms2_data))
  if (length(idx_b) == 0) {stop('b_peak has no ms2!')}
  ms2_b <- ms2_data[idx_b]


  match_list <- mapply(function(x, y){
    temp_ms2_a <- x %>% convertSpectraData()
    temp_ms2_b <- y %>% convertSpectraData()

    if (temp_ms2_a@info$mz >= temp_ms2_b@info$mz) {
      lib_ms2 <- temp_ms2_b
      exp_ms2 <- temp_ms2_a
    } else {
      lib_ms2 <- temp_ms2_a
      exp_ms2 <- temp_ms2_b
    }

    runSpecMatch(obj_ms2_cpd1 = exp_ms2,
      obj_ms2_cpd2 = lib_ms2,
      mz_tol_ms2 = 25,
      scoring_approach = 'dp')
  },
    x = ms2_a,
    y = ms2_b,
    SIMPLIFY = FALSE)


  match_score <- lapply(match_list, function(x)x@info) %>% bind_rows() %>% pull(scoreReverse)
  n_frag <- lapply(match_list, function(x)x@info) %>% bind_rows() %>% pull(n_frag_match)

  pair_matrix <- pair_matrix %>% mutate(ms2_score = match_score, n_frag_match = n_frag)
  names(match_list) <- paste(pair_matrix$peak_a, pair_matrix$peak_b, sep = '@')

  result <- list(pair_matrix = pair_matrix,
    match_list = match_list)

  return(result)
}


# retrieveDiffFormula ----------------------------------------------------------
retrieveDiffFormula <- function(
  from_formula,
  to_formula,
  from_mass,
  to_mass
) {
  diff_formula <- mapply(
    function(from_formula,
      to_formula,
      from_mass,
      to_mass){
      if (from_mass > to_mass) {
        result <- CHNOSZ::makeup(c(from_formula, to_formula), c(1, -1), sum = TRUE)
      } else {
        result <- CHNOSZ::makeup(c(to_formula, from_formula), c(1, -1), sum = TRUE)
      }

      return(result)
    },
    from_formula,
    to_formula,
    from_mass,
    to_mass,
    SIMPLIFY = FALSE)


  diff_formula2 <- sapply(seq_along(diff_formula), function(i){
    # cat(i, ' ')
    result <- try({diff_formula[[i]] %>% as.chemical.formula()},
      silent = TRUE)

    if (class(result) == 'try-error') {
      return(NA)
    } else {
      return(result)
    }
  })

  return(diff_formula2)
}

# showDfLine -------------------------------------------------------------------
showDfLine <- function(x){
  if (nrow(x) != 1) {stop('Only one row is supported')}

  x %>%
    dplyr::mutate_all(as.character) %>%
    tidyr::pivot_longer(everything()) %>%
    print(n = Inf)
}
