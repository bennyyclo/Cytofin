#' Prep CyTOF control for normalization
#'
#' This function establish control reference for CyTOF data normalization.
#'
#' @param metadata_filename panel_filename input_file_dir output_file_dir
#' @return data of summary statistics and concatenated control file
#' @export

anprep <- function(metadata_filename, panel_filename, input_file_dir, output_file_dir){

#read metadata table
getExtension <- function(file){ 
    ex <- strsplit(basename(file), split="\\.")[[1]]
    return(ex[-1])
} 


#parse standard panel
if (getExtension(panel_filename)=="xlsx"){
ref_panel <- readxl::read_excel(panel_filename)
}else if (getExtension(panel_filename)=="csv"){
ref_panel <- read.csv(panel_filename)
}

(lineage_markers <- ref_panel$desc[ref_panel$Lineage == 1])
(functional_markers <- ref_panel$desc[ref_panel$Functional == 1])
all_markers <- c(lineage_markers, functional_markers)

#md_control <- readxl::read_excel(metadata_filename, col_names = TRUE)

if (getExtension(metadata_filename)=="xlsx"){
md_control <- readxl::read_excel(metadata_filename)
}else if (getExtension(metadata_filename)=="csv"){
md_control <- read.csv(metadata_filename)
}

fcs_control <- flowCore::read.flowSet(paste0(input_file_dir,paste0("homogenized_", md_control$filename)), transformation = FALSE, truncate_max_range = FALSE)

#calculate universal mean and variance
asinhNik <- flowCore::arcsinhTransform(a=0,b=0.2)
colname <- flowCore::colnames(fcs_control)
expr_untransformed <- flowCore::fsApply(fcs_control, flowCore::exprs)
tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
fcs_asinh <- flowCore::transform(fcs_control, tlist)
expr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
mean_uni <- apply(expr, 2, mean)
var_uni <- apply(expr, 2, var)
var_uni_mean <- mean(var_uni[all_markers])
mean_uni_mean <- mean(mean_uni[all_markers])

save(var_uni, mean_uni, var_uni_mean, mean_uni_mean, file=paste0(output_file_dir,"Prep_control.RData"))

gc()
filename <- paste0(output_file_dir, "concatenated_control.fcs")
ff <- flowCore::flowFrame(expr)
data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
flowCore::write.FCS(ff, filename)

gc()
filename <- paste0(output_file_dir, "concatenated_control_untransformed.fcs")
ff <- flowCore::flowFrame(expr_untransformed)
data_panel_name <- flowCore::pData(flowCore::parameters(fcs_control[[1]]))$desc
flowCore::pData(flowCore::parameters(ff))$desc <- data_panel_name  
flowCore::write.FCS(ff, filename)
}
