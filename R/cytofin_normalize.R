#' Batch normalize CyTOF plates from heterogeneous sources using external anchors
#'
#' This function batch normalizes CyTOF data from multiple plates (from one or more 
#' experimental cohorts) using external (i.e. "generalized") anchors.
#'
#' @param metadata_path A filepath leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include
#' `filename`, `cohort`, `plate_number`, `patient_id`, `condition`, `is_anchor`, 
#' and `validation`. 
#' 
#' See the vignette for details: \code{vignette("help", package = "cytofin")} 
#' 
#' @param panel_path A file path leading to an .xlsx or .csv file containing 
#' a table of standardized antigen panel information. Columns should include 
#' `metal_name`, `antigen_name`, `antigen_pattern`, 
#' `lineage`, `functional`, and `general`. 
#' 
#' See the vignette for details: \code{vignette("help", package = "cytofin")}
#' 
#' @param anchor_statistics a list produced by the `cytofin_prep_anchors`
#' function or the file path to an .rds object containing anchor reference statistics.
#' 
#' @param input_data_path A folder directory containing the input CyTOF files
#' to be normalized. In most cases, this will be the directory to which the output
#' .fcs files from `cytofin_homogenize` were written.
#' 
#' @param output_data_path A folder directory to which the output (i.e. 
#' batch normalized/batch corrected) .fcs files should be written.
#' 
#' @param mode A string indicating which transformation function should be used 
#' for batch normalization ("meanshift", "meanshift_bulk", "variance", "z_score", 
#' or "beadlike").
#' 
#' @param input_prefix The string that was appended to the name of the input files 
#' of `cytofin_homogenize` to create their corresponding output file names. 
#' Defaults to "homogenized_".
#' 
#' @param output_prefix A string to be appended to the name of each input file 
#' to create the name of the corresponding output file (post-homogenization). 
#' Defaults to "normalized_" (e.g. an input file named "file1.fcs" will correspond to 
#' the output file "normalized_file1.fcs" saved in `output_data_path`).
#' 
#' @param shift_factor The scalar value `a` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arc-sine function:  
#'    
#'    `new_x <- asinh(a + b * x)`.
#' 
#' Defaults to 0. 
#' 
#' @param scale_factor The scalar value `b` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arc-sine function:  
#'    
#'    `new_x <- asinh(a + b * x)`.
#' 
#' Defaults to 0.2. 
#'
#' @return Batch-normalized .fcs files are saved in the directory specified by 
#' `output_data_path`. 
#' 
#' In addition, a data.frame containing information about 
#' each input .fcs file (that can be used for plotting with `cytofin_make_plots`)
#' is returned with the following columns: 
#'    * All of the columns in the input metadata table (located at `metadata_path`)
#'    * __universal_mean__: the universal mean vector to which all files are adjusted 
#'    (will be identical for all input .fcs files)
#'    * __universal_var__: the universal mean vector to which all files are adjusted
#'    (will be identical for all input .fcs files)
#'    * __anchor_mean__: the mean (across all cells) vector for the anchor file associated
#'    with each input .fcs file (i.e. the anchor located on the same plate as the 
#'    input .fcs file) before batch normalization.
#'    * __anchor_var__: the variance (across all cells) vector for the anchor file associated
#'    with each input .fcs file (i.e. the anchor located on the same plate as the 
#'    input .fcs file)
#'    * __mean_b4norm__: the mean (across all cells) vector of the input .fcs file 
#'    before batch normalization. 
#'    * __var_b4norm__: the variance (across all cells) vector of the input .fcs file 
#'    before batch normalization. 
#'    * __mean_norm__: the mean (across all cells) vector of the input .fcs file 
#'    after batch normalization. 
#'    * __var_norm__: the variance (across all cells) vector of the input .fcs file 
#'    after batch normalization.
#'    * __anchor_mean_norm__: the mean (across all cells) vector for the anchor file associated
#'    with each input .fcs file (i.e. the anchor located on the same plate as the 
#'    input .fcs file) after batch normalization.
#'    * __anchor_var_norm__: the variance (across all cells) vector for the anchor file associated
#'    with each input .fcs file (i.e. the anchor located on the same plate as the 
#'    input .fcs file) after batch normalization.
#' 
#' @export
#' 
cytofin_normalize <-
  function(
    metadata_path,
    panel_path,
    anchor_statistics, 
    input_data_path,
    output_data_path,
    mode = c("meanshift", "meanshift_bulk", "variance", "z_score", "beadlike"),
    input_prefix = "homogenized_", 
    output_prefix = "normalized_", 
    shift_factor = 0, 
    scale_factor = 0.2
  ) {
    
    # create output directory
    dir.create(output_data_path)
    
    #read metadata table
    md <- cytofin_read_metadata(metadata_path)
    
    # separate metadata for anchor samples
    md_control <- dplyr::filter(md, is_anchor == 1)

    # if anchor_statistics is a file path
    if (is.character(anchor_statistics)) { 
      anchor_statistics_list <- readr::read_rds(anchor_statistics)
      # else if anchor_statistics is a list 
    } else if (is.list(anchor_statistics)) { 
      anchor_statistics_list <- anchor_statistics
    } else { 
      stop("anchor_statistics must be either a character vector (file path) or a list")
    }
    
    # extract needed values from anchor_statistics_list
    universal_var <- anchor_statistics_list$universal_var
    universal_mean <- anchor_statistics_list$universal_mean
    bulk_var <- anchor_statistics_list$bulk_var
    bulk_mean <- anchor_statistics_list$bulk_mean
    
    # read in standardized panel
    ref_panel <- cytofin_read_panel_info(panel_path = panel_path)
    
    # compile list of all markers to keep during analysis
    lineage_markers <- 
      as.character(ref_panel$metal_name[ref_panel$lineage == 1])
    
    functional_markers <- 
      as.character(ref_panel$metal_name[ref_panel$functional == 1])
    
    all_markers <- c(lineage_markers, functional_markers)
    
    # create transformation functions
    norm_1 <- function(x) {
      y <- universal_mean[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m] - anchor_mean[m] + universal_mean[m]
      return(z)
    } #meanshift
    
    norm_2 <- function(x) {
      y <- universal_mean[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m] - mean(anchor_mean[m]) + mean(universal_mean[m])
      return(z)
    } #meanshift bulk
    
    norm_3 <- function(x) {
      y <- universal_mean[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- 
        (z[m] - anchor_mean[m] + universal_mean[m]) * sqrt(universal_var[m])/sqrt(anchor_var[m])
      return(z)
    } #variance
    
    norm_4 <- function(x) {
      y <- universal_mean[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- 
        (z[m] - anchor_mean[m]) * sqrt(universal_var[m])/sqrt(anchor_var[m]) + universal_mean[m]
      return(z)
    } #z-score
    
    norm_5 <- function(x) {
      y <- universal_mean[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m] * lm(universal_mean[m] ~ anchor_mean[m])$coefficient[[2]]
      return(z)
    } #beadlike

    # find the user-specified normalization function
    if (mode == "meanshift") {
      norm <- norm_1 
    } else if (mode == "meanshift_bulk") {
      norm <- norm_2
    } else if (mode == "variance") {
      norm <- norm_3
    } else if (mode == "z_score") {
      norm <- norm_4
    } else if (mode == "beadlike") {
      norm <- norm_5
    }
    
    # create final data structure 
    result <- 
      dplyr::mutate(
        md, 
        universal_mean = list(0), 
        universal_var = list(0), 
        anchor_mean = list(0), 
        anchor_var = list(0), 
        mean_b4norm = list(0), 
        var_b4norm = list(0), 
        mean_norm = list(0), 
        var_norm = list(0), 
        anchor_mean_norm = list(0), 
        anchor_var_norm = list(0)
      )
    
    # for each file being batch-normalized...
    for (i in 1:length(md$filename)) {
      # calculate adjustment parameters from control plate
      
      # find the anchor file corresponding to the same plate and cohort 
      # as the file being batch normalized
      filename_anchor <- 
        md_control$filename[
          which(
            (md_control$plate_number == md$plate_number[i]) & 
              (md_control$cohort == md$cohort[i])
          )
        ]
      
      # read in the anchor file
      fcs <- 
        flowCore::read.FCS(
          file.path(input_data_path, paste0(input_prefix, filename_anchor)), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      
      # arcsinh transform all columns of the anchor file
      asinh_transform <- 
        flowCore::arcsinhTransform(a = shift_factor, b = scale_factor)
      col_names <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = col_names, tfun = asinh_transform)
      fcs_asinh <- flowCore::transform(fcs, tlist)
    
      # find the mean and variance vector of the anchor file
      anchor_expr <- flowCore::exprs(fcs_asinh)
      anchor_mean <- apply(anchor_expr, 2, mean)
      anchor_var <- apply(anchor_expr, 2, var)
      
      # find the bulk mean and bulk variance of the anchor file
      anchor_bulk_mean <- mean(anchor_mean)
      anchor_bulk_var <- mean(anchor_var)
      
      # normalize the anchor file
      anchor_expr_norm <- 
        t(apply(anchor_expr, 1, norm))
      
      anchor_mean_norm <- apply(anchor_expr_norm, 2, mean)
      anchor_var_norm <- apply(anchor_expr_norm, 2, var)
      
      # normalize the target file
      
      ## read in target file  
      filename <- md$filename[i]
      fcs <- 
        flowCore::read.FCS(
          file.path(input_data_path, paste0(input_prefix, filename)), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      
      # arcsinh-transform all columns of the target file
      col_names <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = col_names, tfun = asinh_transform)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      
      # extract target file's expression matrix before normalization
      expr_b4norm <- flowCore::exprs(fcs_asinh)
      
      # find the target file's un-normalized mean and variance vectors
      mean_b4norm <- apply(expr_b4norm, 2, mean)
      var_b4norm <- apply(expr_b4norm, 2, var)
      
      ## normalize the target file
      expr_norm <- 
        t(apply(expr_b4norm, 1, norm))
      
      # find the mean and variance vectors of the normalized target file
      mean_norm <- apply(expr_norm, 2, mean)
      var_norm <- apply(expr_norm, 2, var)
      
      # create flowFrame to be written as the output .fcs file for this sample
      fcs_norm <- flowCore::flowFrame(expr_norm)
      
      # normalization completed, reverse asinh transformation for final output
      my_rev_asinh <- 
        function(x) {
          rev_asinh(x, shift_factor = shift_factor, scale_factor = scale_factor)
        }
      tlist2 <- flowCore::transformList(from = col_names, tfun = my_rev_asinh)
      fcs_asinh_rev <- flowCore::transform(fcs_norm, tlist2)
      
      # prepare and write out final .fcs file
      flowCore::pData(flowCore::parameters(fcs_asinh_rev))$desc <- 
        flowCore::pData(flowCore::parameters(fcs_asinh))$desc
      fcs_name <- file.path(output_data_path, paste0(output_prefix, filename))
      flowCore::write.FCS(x = fcs_asinh_rev, filename = fcs_name)
      
      # update final data structure 
      result$universal_mean[[i]] <- universal_mean 
      result$universal_var[[i]] <- universal_var
      result$anchor_mean[[i]] <- anchor_mean
      result$anchor_var[[i]] <- anchor_var
      result$mean_b4norm[[i]] <- mean_b4norm
      result$var_b4norm[[i]] <- var_b4norm
      result$mean_norm[[i]] <- mean_norm
      result$var_norm[[i]] <- var_norm
      result$anchor_mean_norm[[i]] <- anchor_mean_norm
      result$anchor_var_norm[[i]] <- anchor_var_norm
    }
    
    # add marker list and arcsinh transformation parameters to the final data structure
    attr(result, which = "shift_factor") <- shift_factor
    attr(result, which = "scale_factor") <- scale_factor
    attr(result, which = "all_markers") <- all_markers
    
    # return result
    return(result)
    
  }
