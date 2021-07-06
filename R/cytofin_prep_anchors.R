#' Prepare CyTOF controls for batch normalization across plates
#'
#' This function calculates reference statistics needed for CytofIn batch normalization. 
#' Specifically, it calculates the universal mean and universal variance vectors
#' of the generalized anchors identified in the metadata file at `metadata_path`; 
#' in addition, it calculates the non-channel-specific bulk mean and bulk variance
#' of the generalized anchors.
#'
#' @param metadata_path A file path leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include
#' `filename`, `cohort`, `plate_number`, `patient_id`, `condition`, `population`, 
#' and `validation`.
#' 
#' See the vignette for details: \code{vignette("help", package = "cytofin")} 
#' 
#' @param panel_path A file path leading to an .xlsx or .csv file containing 
#' a table of standardized antigen panel information. Columns should include 
#' `metal_name`, `antigen_name`, `antigen_pattern`, `lineage`, `functional`, 
#' and `general`.
#'  
#' See the vignette for details: \code{vignette("help", package = "cytofin")}
#' 
#' @param input_data_path A folder directory containing the input CyTOF files
#' to be prepped for normalization. These files should already be homogenized,  
#' and in most cases this will be the directory to which the output
#' .fcs files from `cytofin_homogenize` were written.
#' 
#' @param input_prefix The string that was appended to the name of the input files 
#' of `cytofin_homogenize` to create their corresponding output file names. 
#' Defaults to "homogenized_".
#' 
#' @param output_path A file path specifying where to save the output .rds
#' file containing the statistics calculated from this step and the concatenated 
#' .FCS files containing all cells from the generalized anchor samples. Defaults
#' to "none", in which case no files are saved.
#'
#' @param shift_factor The scalar value `a` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arc-sine function:  
#'    
#'    `new_x <- asinh(a + b*x)`.
#' 
#' Defaults to 0. 
#' 
#' @param scale_factor The scalar value `b` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arc-sine function:  
#'    
#'    `new_x <- asinh(a + b*x)`.
#' 
#' Defaults to 0.2.
#' 
#' @return a `list()` of summary statistics with the following elements:
#' * __universal_var__:  a named numeric vector in which each entry corresponds to the 
#' universal variance of an antigen channel in the homogenized dataset
#' * __universal_mean__:  a named numeric vector in which each entry corresponds to the 
#' universal mean of an antigen channel in the homogenized dataset
#' * __bulk_var__:  The mean of all the channel-specific universal variances 
#' in `universal_var` (a scalar value)
#' * __bulk_mean__:  The mean of all the channel-specific universal means 
#' in `universal_mean` (a scalar value)
#' 
#' 
#' @export
#' 

cytofin_prep_anchors <- function(
  metadata_path, 
  panel_path, 
  input_data_path, 
  input_prefix = "homogenized_",
  output_path = "none", 
  shift_factor = 0, 
  scale_factor = 0.2
) {
  
  # create output directory if needed
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  # read metadata table and select only the anchor samples
  md_control <- 
    dplyr::filter(cytofin_read_metadata(metadata_path), is_anchor == 1)
  
  # read reference panel information
  ref_panel <- cytofin_read_panel_info(panel_path)

  # extract character vectors of the lineage markers' metals and 
  # the functional markers' metals
  lineage_markers <- ref_panel$metal_name[ref_panel$lineage == 1]
  functional_markers <- ref_panel$metal_name[ref_panel$functional == 1]
  all_markers <- c(lineage_markers, functional_markers)
  
  # read in the input data as a flowSet
  fcs_control <- 
    flowCore::read.flowSet(
      file.path(input_data_path, paste0(input_prefix, md_control$filename)), 
      transformation = FALSE, 
      truncate_max_range = FALSE
    )
  
  # arcsinh-transform all data
  asinh_transform <- flowCore::arcsinhTransform(a = shift_factor, b = scale_factor)
  col_names <- flowCore::colnames(fcs_control)
  expr_untransformed <- flowCore::fsApply(fcs_control, flowCore::exprs)
  transform_list <- flowCore::transformList(from = col_names, tfun = asinh_transform)
  fcs_asinh <- flowCore::transform(fcs_control, transform_list)
  expr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
  
  # calculate universal mean and variance
  universal_mean <- apply(expr, 2, mean)
  universal_var <- apply(expr, 2, var)
  
  # calculate the mean and variance of all the channel-specific universal means
  # and variances, respectively
  bulk_var <- mean(universal_var[all_markers])
  bulk_mean <- mean(universal_mean[all_markers])
  
  # collate all reference statistics into a list
  result <-
    list(
      universal_var = universal_var, 
      universal_mean = universal_mean,
      bulk_var = bulk_var, 
      bulk_mean = bulk_mean
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
    filename <- file.path(output_path, "concatenated_control.fcs")
    ff <- flowCore::flowFrame(expr)
    data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
    flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
    flowCore::write.FCS(ff, filename)
    
    # write concatenated control file (untransformed)
    gc()
    filename <- file.path(output_path, "concatenated_control_untransformed.fcs")
    ff <- flowCore::flowFrame(expr_untransformed)
    data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
    flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
    flowCore::write.FCS(ff, filename)
  }
    
    return(result)
  }
  
