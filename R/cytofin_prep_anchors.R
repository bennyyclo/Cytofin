#' Prepare CyTOF controls for batch normalization
#'
#' This function establishes control reference for CyTOF data normalization.
#'
#' @param metadata_path A filepath leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include
#' `filename`, `cohort`, `plate_number`, `patient_id`, `condition`, `population`, 
#' and `validation`. TO DO: Change the names of these columns to more descriptive
#' names and make sure that they are all actually needed. 
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
#' @return Data of summary statistics and concatenated control file
#' 
#' @export
#' 

cytofin_prep_anchors <- function(metadata_path, panel_path, input_data_path) {

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
  fcs_control <- 
    flowCore::read.flowSet( # change to read.FCS
      file.path(input_data_path, paste0("homogenized_", md_control$filename)), 
      transformation = FALSE, 
      truncate_max_range = FALSE
    )
  
  #calculate universal mean and variance
  asinhNik <- flowCore::arcsinhTransform(a=0, b=0.2)
  colname <- flowCore::colnames(fcs_control)
  expr_untransformed <- flowCore::fsApply(fcs_control, flowCore::exprs)
  tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
  fcs_asinh <- flowCore::transform(fcs_control, tlist)
  expr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
  mean_uni <- apply(expr, 2, mean)
  var_uni <- apply(expr, 2, var)
  var_uni_mean <- mean(var_uni[all_markers])
  mean_uni_mean <- mean(mean_uni[all_markers])
  
  save(var_uni, mean_uni, var_uni_mean, mean_uni_mean, file=paste0("Prep_control.RData"))
  
  gc()
  filename <- paste0("concatenated_control.fcs")
  ff <- flowCore::flowFrame(expr)
  data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
  flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
  flowCore::write.FCS(ff, filename)
  
  gc()
  filename <- paste0("concatenated_control_untransformed.fcs")
  ff <- flowCore::flowFrame(expr_untransformed)
  data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
  flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
  flowCore::write.FCS(ff, filename)
}
