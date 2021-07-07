#' Find the extension for a file
#'
#' @param filename A string representing the name of a file in its local directory
#'
#' @return The the file extension of `filename` 
#'
#' @examples
#' \dontrun{
#' # example file name
#' my_filename <- "my_file.txt"
#' 
#' # find and print the extension
#' my_extension <- getExtension(my_filename)
#' print(my_extension)
#' }
get_extension <- function(filename) {
  ex <- strsplit(basename(filename), split="\\.")[[1]]
  return(ex[[-1]])
}


#' Alter a flowFrame to only include data from channels in a reference panel
#'
#' @param fcs_raw A flowFrame containing unprocessed CyTOF data
#' @param ref_panel A data.frame representing the reference panel data for a 
#' cytofin analysis.
#'
#' @return a homogenized flowFrame
#'
homogenize_flowFrame <- function(fcs_raw, ref_panel) {
  
  #extract some needed values from the raw fcs data and the reference panel
  ref_markers <- ref_panel$antigen_name
  ref_metals <- ref_panel$metal_name
  
  panel_fcs <- flowCore::pData(flowCore::parameters(fcs_raw))
  panel_markers <- panel_fcs$desc
  panel_metals <- as.character(panel_fcs$name)
  panel_rownames <- row.names(panel_fcs)
  
  # create new flowFrame to be modified
  fcs <- fcs_raw
  
  # only keep the markers/metals that are present in the reference marker list
  panel_markers_to_keep <- intersect(panel_markers, ref_markers)
  panel_metals_to_keep <- panel_metals[panel_markers %in% panel_markers_to_keep]
  
  # create dictionary to look up which metals from the reference panel 
  # correspond to shared antigens with the FCS file's panel (which may be on 
  # different metals)
  names(ref_metals) <- ref_markers
  
  # perform lookup to "rename" metals in the FCS file's panel to the standard
  # metal name in the reference
  new_panel_metals <- as.character(ref_metals[panel_markers_to_keep])
  
  # remove columns not present in ref_panel from final fcs file
  expr <- flowCore::exprs(fcs)
  new_expr <- expr[, panel_metals_to_keep] 
  
  # rename metals using the looked-up values
  colnames(new_expr) <- new_panel_metals
  
  # sort columns into the order in the reference panel
  final_expr <- new_expr[ , ref_metals]
  flowCore::exprs(fcs) <- final_expr
  
  # return result 
  return(fcs)
  
}


#' Read in a cytofin metadata file
#' 
#' This function reads a cytofin metadata file from a connection 
#' that points to a .csv or a .xlsx file
#'
#' @param metadata_path A filepath leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include
#' `filename`, `cohort`, `plate_number`, `patient_id`, `condition`, `population`, 
#' and `validation`. TO DO: Change the names of these columns to more descriptive
#' names and make sure that they are all actually needed. 
#' See the vignette for details: \code{vignette("help", package = "cytofin")}
#'
#' @return A data.frame containing the metadata information in the 
#' file stored at `metadata_path`. 
#'
#' @examples
#' \dontrun{
#' my_path <- file.path("~", "foo", "bar", "metadata.csv")
#' my_metadata <- cytofin:::cytofin_read_metadata(my_path)
#' }
cytofin_read_metadata <- function(metadata_path) {
  
  if (get_extension(metadata_path) == "xlsx") {
    md <- readxl::read_excel(metadata_path)
  } else if (get_extension(metadata_path) == "csv") {
    md <- read.csv(metadata_path)
  } else { 
    # throw error if the wrong kind of file is given
    stop("metadata_path must point to an .xlsx or .csv file")
  }
  
  # trim whitespace from all strings in metadata
  md <- data.frame(lapply(md, trimws), stringsAsFactors = FALSE)
  
  return(md)
}


#' Read in a cytofin reference panel information
#' 
#' This function reads cytofin reference panel information from a connection 
#' that points to a .csv or a .xlsx file
#'
#' @param panel_path A file path leading to an .xlsx or .csv file containing 
#' a table of standardized antigen panel information. Columns should include 
#' `desc`, `range`, `metal_pattern`, `antigen_pattern`, `Lineage`, `Functional`, 
#' and `General`. TO DO: Change the names of these columns to more descriptive
#' names and make sure that they are all actually needed. 
#' See the vignette for details: \code{vignette("help", package = "cytofin")}
#'
#' @return A data.frame containing the reference panel information in the 
#' file stored at `panel_path`. 
#'
#' @examples
#' \dontrun{
#' my_path <- file.path("~", "foo", "bar", "panel.csv")
#' my_metadata <- cytofin:::cytofin_read_panel_info(my_path)
#' }
cytofin_read_panel_info <- function(panel_path) {
  
  if (get_extension(panel_path) == "xlsx") {
    ref_panel <- readxl::read_excel(panel_path)
  } else if (get_extension(panel_path) == "csv") {
    ref_panel <- read.csv(panel_path)
  } else { 
    # throw error if the wrong kind of file is given
    stop("panel_path must point to an .xlsx or .csv file")
  }
  
  # trim whitespace from all strings in reference panel
  ref_panel <- data.frame(lapply(ref_panel, trimws), stringsAsFactors = FALSE)
  
  return(ref_panel)
}


#' Reverses arcsinh transformation with cofactor `scale_factor` and a shift of `shift_factor`.
#'
#' @param x A numeric vector.
#' 
#' @param shift_factor The scalar value `a` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arcsinh function:  
#'    `new_x <- asinh(a + b * x)`.
#' 
#' @param scale_factor The scalar value `b` in the following equation used to 
#' transform CyTOF raw data ion counts using the hyperbolic arcsinh function:  
#'    `new_x <- asinh(a + b * x)`.
#'
#' @return A numeric vector after undergoing reverse 
#' arcsinh transformation 
#' 
#'
rev_asinh <- function(x, shift_factor, scale_factor) {
  
  new_x <- (sinh(x) - shift_factor) / scale_factor
  return(new_x)
  
}
