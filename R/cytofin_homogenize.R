#' Homogenize CyTOF data panel
#'
#' This function homogenizes CyTOF data (.fcs files) from heterogeneous sources 
#' according to the standard panel in `panel_path.`
#'
#' @param metadata_path A filepath leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names. Columns should include ____. 
#' Must be in the current directory. 
#' 
#' @param panel_path A filepath leading to an .xlsx or .csv file containing 
#' a table of standardized antigen panel information. Columns should include ____.  
#' Must be in the current directory.
#' 
#' @param input_data_path A folder directory containing the input CyTOF files
#' to be homogenized.
#' 
#' @param output_data_path A folder directory containing output homogenized files.
#' 
#' @return `cytofin_homogenize` doesn't return anything. Instead, it has the 
#' side-effect of saving homogenized files (in .fcs format) to the directory 
#' specified with `output_data_path`. Each of the saved files will contain 
#' homogenized, user-defined channels according to details specified in the 
#' file at `panel_path.`
#' 
#' @export
#' 
#' @examples
#' 
#' # For a complete example of the `cytofin` workflow, 
#' # see the packages vignette by running `vignette(package = "cytofin")`
#' 
#' # download and store the metadata file
#' TO DO 
#' 
#' # download and store the panel file
#' TO DO 
#' 
#' # download and store the input data in a temporary directory
#' TO DO 
#' 
#' # specify a temporary output directory
#' TO DO 
#' 
#' # run homogenization
#' TO DO 

cytofin_homogenize <- 
  function(metadata_path, panel_path, input_data_path, output_data_path) {
    
    # create output directory for homogenized .fcs files
    dir.create(output_data_path)
    
    # read metadata table
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
    
    # extract components of the metadata (why???)
    filename <- md$filename
    cohort <- md$cohort
    plate_number <- md$plate_number
    patient_id <- md$patient_id
    condition <- md$condition
    population <- md$population
    
    print(md)
    
    # read reference panel information
    if (get_extension(panel_path) == "xlsx") {
      ref_panel <- readxl::read_excel(panel_path)
    } else if (get_extension(panel_path) == "csv") {
      ref_panel <- read.csv(panel_path)
    }
    ref_panel <- data.frame(lapply(ref_panel, trimws), stringsAsFactors = FALSE)
    
    # for all files in the input directory
    for (file in md$filename) {
      # read in FCS file
      fcs_raw <- 
        flowCore::read.FCS(
          filename = file.path(input_data_path, file), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      cat("filename:", file, "\n")
      
      # parse panel in FCS files
      data_panel_antigens <- 
        flowCore::pData(flowCore::parameters(fcs_raw))$desc
      
      data_panel_metals <- 
        flowCore::pData(flowCore::parameters(fcs_raw))$name
      
      # for each channel in the reference panel
      for (i in 1:length(ref_panel$range)) {
        cat(i, "\n")
        tryCatch(
          {
            # extract the antigen name in the reference and its corresponding regex
            ref_antigen <- ref_panel$range[[i]]
            ref_antigen_regex <- ref_panel$antigen_pattern[[i]] #regex
            
            ### are neither of these used?
            ref_metal <- ref_panel$desc[[i]]
            ref_metal_regex <- ref_panel$metal_pattern[[i]]
            
            # Find the index of the data antigen corresponding to the reference antigen
            data_antigen_index <- 
              stringr::str_detect(
                string = tidyr::replace_na(data_panel_antigens,''), 
                pattern = ref_antigen_regex
              )
            # store the name of the data antigen for reporting
            data_antigen <- data_panel_antigens[data_antigen_index]
            
            # if there was a match with the reference antigen's regex
            if (max(data_antigen_index) == 1) {
              # rename the data antigen in the flowFrame using the reference name
              flowCore::pData(flowCore::parameters(fcs_raw))$desc[data_antigen_index] <- 
                ref_antigen
              # otherwise
            } else {
              # do nothing???
            }
            
            cat(
              "matched data antigen: ",
              data_antigen,
              "\nwith the reference antigen: ",
              ref_antigen,
              "\nusing the regex: ",
              ref_antigen_regex,
              "\n"
            )
          }, 
          # what's going on with this error function???
          error = 
            function(e) {
              txt <- 
                paste(
                  md$filename, "item", i , 
                  "data_antigen", data_antigen, "ref_antigen",
                  ref_antigen, "ref_antigen_pattern", ref_antigen_pattern,
                  as.character(e)
                )
              cat(txt,"\n")
            }
        )
        
        fcs_raw
      }
      
      # finalize the fcs file to write as output
      fcs <- homogenize_flowFrame(fcs_raw, ref_panel)
      
      # write output fcs file to the specified directory
      filename <- paste0(output_data_path,"homogenized_", file)
      flowCore::write.FCS(fcs, filename)
      
    }
    
  }
