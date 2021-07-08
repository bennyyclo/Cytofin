#' Batch normalize CyTOF plates from heterogeneous sources using stable channels
#'
#' This function batch normalizes CyTOF data from multiple plates (from one or more 
#' experimental cohorts) by computing the non-redundancy score (NRS) for each 
#' channel in the dataset, then using the most redundant (i.e. the "most stable") 
#' channels as a reference for batch normalization.
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
#' @param input_data_path A folder directory containing the input CyTOF files
#' to be normalized. In most cases, this will be the directory to which the output
#' .fcs files from `cytofin_homogenize` were written.
#' 
#' @param output_data_path A folder directory to which the output (i.e. 
#' batch normalized/batch corrected) .fcs files should be written.
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
#' transform CyTOF raw data ion counts using the hyperbolic arcsinh function:
#'   
#'    `new_x <- asinh(a + b * x)`.
#'    
#' Defaults to 0. 
#' 
#' @param scale_factor The scalar value `b` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arcsinh function: 
#'  
#'    `new_x <- asinh(a + b * x)`.
#'    
#' Defaults to 0.2. 
#' 
#' @param nchannels An integer representing the number of most stable channels to
#' use during batch normalization. Defaults to 3.
#' 
#' @param make_plot A boolean value indicating if a plot depicting the non-
#' redundancy scores of each marker in each .fcs file being batch normalized
#' should be plotted as a side-effect of the function call. Defaults to FALSE.
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
#'    input .fcs file)
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
cytofin_normalize_nrs <-
  function(
    metadata_path,
    panel_path,
    input_data_path,
    output_data_path,
    input_prefix = "homogenized_", 
    output_prefix = "normalized_", 
    shift_factor = 0, 
    scale_factor = 0.2,
    nchannels = 3, 
    make_plot = FALSE
  ) {
    
    # create output directory
    dir.create(output_data_path, showWarnings = FALSE, recursive = TRUE)
  
    #read metadata table
    md <- cytofin_read_metadata(metadata_path)
    
    # read in standardized panel
    ref_panel <- cytofin_read_panel_info(panel_path = panel_path)
    
    # compile list of all markers to keep during analysis
    lineage_markers <- 
      as.character(ref_panel$metal_name[ref_panel$lineage == 1])
    
    functional_markers <- 
      as.character(ref_panel$metal_name[ref_panel$functional == 1])
    
    all_markers <- c(lineage_markers, functional_markers)
    
    # transformation function
    norm <- function(x) {
      y <- universal_mean[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- 
        z[m] - 
        mean(mean_ctr[selected_markers]) + 
        mean(universal_mean[selected_markers])
      return(z)
    } # meanshift bulk
    
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
    
    ## create function to compute non-redundancy score for all channels
    NRS <- function(x, ncomp = 3) {
      pr <- prcomp(x, center = TRUE, scale. = FALSE)
      score <- 
        rowSums(
          outer(
            rep(1, ncol(x)),
            pr$sdev[1:ncomp]^2
          ) * 
            abs(pr$rotation[, 1:ncomp]))
      return(score)
    }
    
    # read in all .fcs files to be normalized
    fcs <- 
      flowCore::read.flowSet(
        file.path(input_data_path, paste0(input_prefix, md$filename)), 
        transformation = FALSE, 
        truncate_max_range = FALSE
      )
    
    # arcsinh transform all channels of the input .fcs files
    asinh_transform <- 
      flowCore::arcsinhTransform(a = shift_factor, b = scale_factor)
    col_names <- flowCore::colnames(fcs)
    tlist <- flowCore::transformList(from = col_names, tfun = asinh_transform)
    fcs_asinh <- flowCore::transform(fcs, tlist)
    
    # find the mean and variance vector of all cells in the dataset
    expr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
    universal_mean <- apply(expr, 2, mean)
    universal_var <- apply(expr, 2, var)
    
    # calculate non-redundancy scores for each antigen in each .fcs file
    nrs_sample <-
      flowCore::fsApply(fcs_asinh[, all_markers], NRS, use.exprs = TRUE)
    
    # find mean non-redundancy scores for each antigen across all samples 
    colnames(nrs_sample) <-
      as.character(ref_panel$antigen_name[match((colnames(nrs_sample)), ref_panel$metal_name)])
    nrs <- colMeans(nrs_sample, na.rm = TRUE)
    
    nrs_sample <- data.frame(nrs_sample)
    markers_ord <- names(sort(nrs, decreasing = TRUE))
    nrs_sample <- data.frame(nrs_sample)
    nrs_sample$sample_id <- rownames(nrs_sample)
    
    if (make_plot) {
      # make data.frame for plotting
      ggdf <-
        reshape2::melt(
          nrs_sample, 
          id.var = "sample_id", 
          value.name = "nrs", 
          variable.name = "antigen"
        )
      
      ggdf$antigen <- 
        factor(ggdf$antigen, levels = markers_ord)
      
      # make plot
      p <- 
        ggplot2::ggplot(ggdf, ggplot2::aes(x = antigen, y = nrs)) +
        ggplot2::geom_point(
          ggplot2::aes(color = sample_id),
          alpha = 0.9,
          position = ggplot2::position_jitter(width = 0.3, height = 0)
        ) +
        ggplot2::geom_boxplot(outlier.color = NA, fill = NA) +
        ggplot2::stat_summary(fun = "mean", geom = "point", shape = 21, fill = "white") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
        )
      
      print(p)
      
    }
    
    ####--------#####
    
    # select nchannels antigens with the lowest NRS for calibration
    selected_markers <- names(sort(nrs, decreasing = FALSE))[1:nchannels]
    
    # find the metal names corresponding to the chosen antigens
    selected_markers <- as.character(ref_panel$metal_name[match(selected_markers, ref_panel$antigen_name)])
    
    for (i in 1:length(md$filename)) {
      # calculate adjustment parameters from control plate
      
      # read in .fcs file
      filename_ctr <- md$filename[i]
      fcs <- 
        flowCore::read.FCS(
          file.path(input_data_path, paste0(input_prefix, filename_ctr)), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      
      # asinh-transform .fcs file and subset out only its selected channels
      col_names <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = col_names, tfun = asinh_transform)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      expr_ctr <- flowCore::exprs(fcs_asinh[, selected_markers])
      
      # find the mean and variance of the nchannels selected channels in the .fcs file
      mean_ctr <- apply(expr_ctr, 2, mean)
      mean_ctr_mean <- mean(mean_ctr)
      var_ctr <- apply(expr_ctr, 2, var)
      var_ctr_mean <- mean(var_ctr)

      # batch normalize the .fcs channels
      expr_ctr_norm <- 
        t(apply(flowCore::exprs(fcs_asinh), 1, norm))
      
      # find the mean and variance vectors of the normalized input file 
      mean_ctr_norm <- apply(expr_ctr_norm[, selected_markers], 2, mean)
      var_ctr_norm <- apply(expr_ctr_norm[, selected_markers], 2, var)
      
      # find the bulk mean and variance of the normalized input file
      mean_ctr_norm_mean <- mean(mean_ctr_norm)
      var_ctr_norm_mean <- mean(var_ctr_norm)
      
      # normalize the target plate
      ## before
      
      # read in input .fcs file
      filename <- md$filename[i]
      fcs <- 
        flowCore::read.FCS(
          file.path(input_data_path, paste0(input_prefix, filename)), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      
      # asinh-transform input .fcs file
      col_names <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = col_names, tfun = asinh_transform)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      expr_b4norm <- flowCore::exprs(fcs_asinh)
      
      # find the mean and variance vectors before batch correction
      mean_b4norm <- apply(expr_b4norm, 2, mean)
      var_b4norm <- apply(expr_b4norm, 2, var)
      
      # find the bulk mean and bulk variance before batch correction
      mean_b4norm_mean <- mean(mean_b4norm)
      var_b4norm_mean <- mean(var_b4norm)
      
      ## after
      expr_norm <- 
        t(apply(flowCore::exprs(fcs_asinh), 1, norm))
      
      # calculate mean and variance vectors for the normalized input .fcs file
      mean_norm <- apply(expr_norm, 2, mean)
      var_norm <- apply(expr_norm, 2, var)
      
      # calculate bulk mean and variance values for the normalized input .fcs file
      mean_norm_mean <- mean(mean_norm)
      var_norm_mean <- mean(var_norm)
      
      # create output flowFrame 
      fcs_norm <- flowCore::flowFrame(expr_norm)
      
      # normalization completed, reverse transformation
      my_rev_asinh <- 
        function(x) {
          rev_asinh(x, shift_factor = shift_factor, scale_factor = scale_factor)
        }
      tlist2 <- flowCore::transformList(from = col_names, tfun = my_rev_asinh)
      fcs_asinh_rev <- flowCore::transform(fcs_norm, tlist2)
      
      
      flowCore::pData(flowCore::parameters(fcs_asinh_rev))$desc <- 
        flowCore::pData(flowCore::parameters(fcs_asinh))$desc
      
      # write out output .fcs file
      fcs_name <- file.path(output_data_path, paste0(output_prefix, filename))
      flowCore::write.FCS(fcs_asinh_rev, fcs_name)
      
      # update final data structure
      result$universal_mean[[i]] <- universal_mean 
      result$universal_var[[i]] <- universal_var
      result$anchor_mean[[i]] <- mean_ctr
      result$anchor_var[[i]] <- var_ctr
      result$mean_b4norm[[i]] <- mean_b4norm
      result$var_b4norm[[i]] <- var_b4norm
      result$mean_norm[[i]] <- mean_norm
      result$var_norm[[i]] <- var_norm
      result$anchor_mean_norm[[i]] <- mean_ctr_norm
      result$anchor_var_norm[[i]] <- var_ctr_norm
    }
    
    # add marker list and arcsinh transformation parameters to the final data structure
    attr(result, which = "shift_factor") <- shift_factor
    attr(result, which = "scale_factor") <- scale_factor
    attr(result, which = "all_markers") <- all_markers
    
    # return result
    return(result)
    
    
  }
