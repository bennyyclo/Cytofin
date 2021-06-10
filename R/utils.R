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


#' Title
#'
#' @param fcs_raw 
#' @param ref_panel 
#'
#' @return
#' @export
#'
#' @examples
#' TO DO: Insert examples
homogenize_flowFrame <- function(fcs_raw, ref_panel) {
  
  #extract some needed values from the raw fcs data and the reference panel
  ref_markers <- ref_panel$range
  ref_metals <- ref_panel$desc
  
  panel_fcs <- flowCore::pData(flowCore::parameters(fcs_raw))
  fcs <- fcs_raw
  
  # why is this being done???
  # what does the prefix ori_ mean???
  ori_all_metals <- panel_fcs$name[panel_fcs$desc %in% ref_markers]
  ori_vals <- row.names(panel_fcs)[panel_fcs$desc %in% ref_markers]
  
  # sort
  ori_all_metals_sorted <- 
    ori_all_metals[order(as.numeric(gsub("[^[:digit:]]", "", ori_vals)))]
  ref_metals_sorted <- 
    ref_metals[order(as.numeric(gsub("[^[:digit:]]", "", ori_vals)))]
  
  # remove columns not present in ref_panel from final fcs file
  expr <- flowCore::exprs(fcs)
  new_expr <- expr[, ori_all_metals_sorted]
  final_expr <- new_expr[, colnames(new_expr) %in% ref_metals]
  flowCore::exprs(fcs) <- final_expr
  
  # return result 
  return(fcs)
  
}
