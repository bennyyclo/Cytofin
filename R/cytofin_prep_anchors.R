#' Prepare CyTOF controls for batch normalization
#'
#' This function calculates reference statistics for CyTOF data normalization by 
#' finding the universal means and universal variances of the generalized 
#' anchors identified in the metadata file at `metadata_path`.
#'
#' @param metadata_path A file path leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include
#' `filename`, `cohort`, `plate_number`, `patient_id`, `condition`, `population`, 
#' and `validation`. TO DO: Change the names of these columns to more descriptive
#' names and make sure that they are all actually needed. 
#' TO DO: Add to the vignette how users are supposed to identify which rows from 
#' the large metadata table should be used as anchors.
#' See the vignette for details: \code{vignette("help", package = "cytofin")} 
#' 
#' @param panel_path A file path leading to an .xlsx or .csv file containing 
#' a table of standardized antigen panel information. Columns should include 
#' `desc`, `range`, `metal_pattern`, `antigen_pattern`, `Lineage`, `Functional`, 
#' and `General`. TO DO: Change the names of these columns to more descriptive
#' names and make sure that they are all actually needed. 
#' See the vignette for details: \code{vignette("help", package = "cytofin")}
#' 
#' @param input_data_path A folder directory containing the input CyTOF files
#' to be prepped for normalization
#' 
#' @param output_path A file path specifying where to save the output .rds
#' file containing the statistics calculated from this step and the concatenated 
#' .FCS files containing all cells from the generalized anchor samples. Defaults
#' to "none", in which case no files are saved.
#' 
#' @return a `list()` of summary statistics with the following elements:
#' * var_uni:  a named numeric vector in which each entry corresponds to the 
#' universal variance of an antigen channel in the homogenized dataset
#' * mean_uni:  a named numeric vector in which each entry corresponds to the 
#' universal mean of an antigen channel in the homogenized dataset
#' * var_uni_mean:  The mean of all the channel-specific universal variances 
#' in `var_uni` (a scalar value)
#' * mean_uni_mean:  The mean of all the channel-specific universal means 
#' in `mean_uni` (a scalar value)
#' 
#' 
#' @export
#' 

cytofin_prep_anchors <- function(
  metadata_path, 
  panel_path, 
  input_data_path, 
  output_path = "none"
) {
  
  # read metadata table
  md_control <- cytofin_read_metadata(metadata_path)
  
  # read reference panel information
  ref_panel <- cytofin_read_panel_info(panel_path)

  # extract character vectors of the lineage markers' metals and 
  # the functional markers' metals
  lineage_markers <- ref_panel$desc[ref_panel$Lineage == 1]
  functional_markers <- ref_panel$desc[ref_panel$Functional == 1]
  all_markers <- c(lineage_markers, functional_markers)
  
  # read in the input data as a flowSet
  # TO DO change the hard-coding of the "homogenized_" tag
  fcs_control <- 
    flowCore::read.flowSet(
      file.path(input_data_path, paste0("homogenized_", md_control$filename)), 
      transformation = FALSE, 
      truncate_max_range = FALSE
    )
  
  # arcsinh-transform all data
  asinh_transform <- flowCore::arcsinhTransform(a = 0, b = 0.2)
  col_names <- flowCore::colnames(fcs_control)
  expr_untransformed <- flowCore::fsApply(fcs_control, flowCore::exprs)
  transform_list <- flowCore::transformList(from = col_names, tfun = asinh_transform)
  fcs_asinh <- flowCore::transform(fcs_control, transform_list)
  expr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
  
  # calculate universal mean and variance
  # TO DO: This universal mean and variance are not calculated in the way that 
  #        They are reported in the manuscript. (This is just a mean of all cells
  #        in the dataset, but the universal mean that we report calculating is
  #        supposed to be the mean of the means of each sample independently).
  mean_uni <- apply(expr, 2, mean)
  var_uni <- apply(expr, 2, var)
  
  # calculate the mean and variance of all the channel-specific universal means
  # and variances, respectively
  var_uni_mean <- mean(var_uni[all_markers])
  mean_uni_mean <- mean(mean_uni[all_markers])
  
  # collate all reference statistics into a list
  result <-
    list(
      var_uni = var_uni, 
      mean_uni = mean_uni,
      var_uni_mean = var_uni_mean, 
      mean_uni_mean = mean_uni_mean
    )
  
  # if the user wants to store intermediate files
  if (output_path != "none") {
    
    # save universal mean and variance information
    readr::write_rds(
      x = result, 
      file = file.path(output_path, "anchor_statistics.rds")
    )
    
    # write concatenated control file (asinh-transformed)
    gc()
    filename <- file.path(out_path, "concatenated_control.fcs")
    ff <- flowCore::flowFrame(expr)
    data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
    flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
    flowCore::write.FCS(ff, filename)
    
    # write concatenated control file (untransformed)
    gc()
    filename <- file.path(out_path, "concatenated_control_untransformed.fcs")
    ff <- flowCore::flowFrame(expr_untransformed)
    data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
    flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
    flowCore::write.FCS(ff, filename)
  }
    
    return(result)
  }
  
