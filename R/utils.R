#' Find the extension for a file
#'
#' @param filename A string representing the name of a file in its local directory
#'
#' @return The the file extension of `filename` 
#'
#' @examples
#' # example file name
#' my_filename <- "my_file.txt"
#' 
#' # find and print the extension
#' my_extension <- getExtension(my_filename)
#' print(my_extension)
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
#' @return
#'
#' @examples
#' TO DO: Insert examples
homogenize_flowFrame <- function(fcs_raw, ref_panel) {
  
  #extract some needed values from the raw fcs data and the reference panel
  ref_markers <- ref_panel$range
  ref_metals <- ref_panel$desc
  
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
