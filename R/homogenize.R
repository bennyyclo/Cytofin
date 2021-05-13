#' Homogenize CyTOF data panel
#'
#' This function homogenizes CyTOF data (.fcs files) from heterogeneous sources 
#' according to the standard panel in panel_filename.
#'
#' @param metadata_filename A filepath leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include ____. 
#' Must be in the current directory. 
#' 
#' @param panel_filename A filepath leading to an .xlsx or .csv file containing 
#' a table of standardized antigen panel information. Columns should include ____.  
#' Must be in the current directory.
#' 
#' @param input_file_dir A folder directory containing the input CyTOF files
#' to be homogenized.
#' 
#' @param output_file_dir A folder directory containing output homogenized files.
#' 
#' @return Homogenized files (in .fcs format) with user-defined channels according
#' to the file named panel_filename.
#' 
#' @export

homogenize <- 
  function(metadata_filename, panel_filename, input_file_dir, output_file_dir) {
    
    # create output directory for homogenized .fcs files
    dir.create(output_file_dir)
    
    # read metadata table
    getExtension <- function(file){ 
      ex <- strsplit(basename(file), split="\\.")[[1]]
      return(ex[-1])
    } 
    if (getExtension(metadata_filename)=="xlsx") {
      md <- readxl::read_excel(metadata_filename)
    } else if (getExtension(metadata_filename)=="csv") {
      md <- read.csv(metadata_filename)
    } else { 
      #throw error if the wrong kind of file is given
      stop("metadata_filename must point to an .xlsx or .csv file")
    }
    
    md <- data.frame(lapply(md, trimws), stringsAsFactors = FALSE)
    filename <- md$filename
    cohort <- md$cohort
    plate_number <- md$plate_number
    patient_id <- md$patient_id
    condition <- md$condition
    population <- md$population
    
    print(md)
    
    #read excel file
    if (getExtension(panel_filename)=="xlsx"){
      ref_panel <- readxl::read_excel(panel_filename)
    }else if (getExtension(panel_filename)=="csv"){
      ref_panel <- read.csv(panel_filename)
    }
    
    
    
    ref_panel <- data.frame(lapply(ref_panel, trimws), stringsAsFactors = FALSE)
    
    print(ref_panel)
    
    
    for (file in md$filename){
      #read in flowSet
      fcs_raw <- flowCore::read.flowSet(paste0(input_file_dir,file), transformation = FALSE, truncate_max_range = FALSE)
      cat("filename:",file,"\n")
      
      #parse standard panel
      data_panel <- flowCore::pData(flowCore::parameters(fcs_raw[[1]]))$desc
      data_panel_name <- flowCore::pData(flowCore::parameters(fcs_raw[[1]]))$name 
      
      #for each channel in the reference panel
      for (i in 1:length(ref_panel$range)) {
        cat(i, "\n")
        tryCatch({
          #correct antigen name
          ref_antigen <- ref_panel$range[i]
          ref_antigen_pattern <- ref_panel$antigen_pattern[i] #regex
          ref_metal <- ref_panel$desc[i]
          ref_metal_pattern <- ref_panel$metal_pattern[i]
          data_antigen <- data_panel[stringr::str_detect(tidyr::replace_na(data_panel,''),ref_antigen_pattern)]
          
          if (max(stringr::str_detect(tidyr::replace_na(data_panel,''),ref_antigen_pattern))==1){
            flowCore::pData(flowCore::parameters(fcs_raw[[1]]))$desc[stringr::str_detect(tidyr::replace_na(data_panel,''),ref_antigen_pattern)] <- ref_antigen
          }else{
            #pData(parameters(fcs_raw[[1]]))$desc[str_detect(replace_na(data_panel_name,''),ref_metal_pattern)] <- ref_antigen 
          }
          
          cat("matched data_antigen:",data_antigen,"ref_antigen:",ref_antigen,"ref_antigen_pattern",ref_antigen_pattern,"\n")  
        }, error = function(e) {
          txt <- paste(md$filename,"item",i,"data_antigen",data_antigen,"ref_antigen",ref_antigen,"ref_antigen_pattern",ref_antigen_pattern,as.character(e))
          cat(txt,"\n")
        })
        fcs_raw[[1]]
      }
      
      
      all_markers <- ref_panel$range
      panel_fcs <- flowCore::pData(flowCore::parameters(fcs_raw[[1]]))
      all_metals <- ref_panel$desc
      panel_metals <- ref_panel$desc
      ori_all_metals <- panel_fcs$name[match(all_markers, panel_fcs$desc)]
      ori_vals <- row.names(panel_fcs)[match(all_markers, panel_fcs$desc)]
      ori_all_metals_sorted <- ori_all_metals[order(as.numeric(gsub("[^[:digit:]]", "", ori_vals)))]
      all_metals_sorted <- all_metals[order(as.numeric(gsub("[^[:digit:]]", "", ori_vals)))]
      
      fcs <- flowCore::fsApply(fcs_raw, function(x){
        
        expr <- flowCore::exprs(x)
        expr <- expr[, ori_all_metals_sorted]
        flowCore::exprs(x) <- expr
        x
      })
      flowCore::colnames(fcs) <- all_metals_sorted
      
      fcs2 <- flowCore::fsApply(fcs, function(x){
        
        expr <- flowCore::exprs(x)
        expr <- expr[, panel_metals]
        flowCore::exprs(x) <- expr
        x
      })
      
      
      filename <- paste0(output_file_dir,"homogenized_",file)
      flowCore::write.FCS(fcs2[[1]], filename)
      fcs2[[1]]
      ##
    }
    
  }
