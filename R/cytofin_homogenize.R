#' Homogenize CyTOF channels names using a consensus antigen panel
#'
#' This function homogenizes CyTOF data (.fcs files) from heterogeneous sources 
#' according to the standard panel in a .csv file located at `panel_path.`
#'
#' @param metadata_path A file path leading to an .xlsx or .csv file 
#' containing a table of CyTOF file (.fcs file) names in the first column (`filename`)
#' and additional information about each .fcs file in subsequent columns. 
#' Columns should include `filename`, `cohort`, `plate_number`, `patient_id`, 
#' `condition`, `is_anchor`, and `validation`. 
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
#' @param input_data_path A folder directory containing the input .fcs files
#' to be homogenized.
#' 
#' @param output_data_path A folder directory to which the output (i.e. 
#' homogenized) .fcs files should be written.
#' 
#' @param prefix A string appended to the name of each input file to create the 
#' name of the corresponding output file (post-homogenization). Defaults to 
#' "homogenized_" (e.g. an input file named "file1.fcs" will correspond to 
#' the output file "homogenized_file1.fcs" saved in `output_data_path`). 
#' 
#' @param verbose A boolean value indicating whether progress message should be 
#' printed to the console during homogenization. Defaults to FALSE.
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
#' # see the packages vignette by running the following: 
#' `vignette("help", package = "cytofin")`
#' 
cytofin_homogenize <- 
  function(
    metadata_path, 
    panel_path, 
    input_data_path, 
    output_data_path, 
    prefix = "homogenized_",
    verbose = FALSE
  ) {
    
    # create output directory for homogenized .fcs files
    dir.create(output_data_path, showWarnings = FALSE, recursive = TRUE)
    
    # read metadata table
    md <- cytofin_read_metadata(metadata_path)
    
    # read reference panel information
    ref_panel <- cytofin_read_panel_info(panel_path)
    
    # for all files in the input directory
    for (file in md$filename) {
      # read in FCS file
      sink(file = "/dev/null")
      fcs_raw <- 
        flowCore::read.FCS(
          filename = file.path(input_data_path, file), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      sink()
      if(verbose) {
        cat("filename:", file, "\n")
      }
      
      # parse panel in FCS files
      data_panel_antigens <- 
        flowCore::pData(flowCore::parameters(fcs_raw))$desc
      
      data_panel_metals <- 
        flowCore::pData(flowCore::parameters(fcs_raw))$name
      
      # for each channel in the reference panel
      for (i in 1:length(ref_panel$antigen_name)) {
        tryCatch(
          {
            # extract the antigen name in the reference and its corresponding regex
            ref_antigen <- ref_panel$antigen_name[[i]]
            ref_antigen_regex <- ref_panel$antigen_pattern[[i]]
            
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
              # rename the data antigen in the flowFrame using the reference antigen name
              flowCore::pData(flowCore::parameters(fcs_raw))$desc[data_antigen_index] <- 
                ref_antigen
              # otherwise
            }
            
            # report what was matched if verbose
            if(verbose) {
              cat(
                "matched data antigen: ",
                data_antigen,
                "\nwith the reference antigen: ",
                ref_antigen,
                "\nusing the regex: ",
                ref_antigen_regex,
                "\n"
              )
            }
          }, 
          # if an error is encountered, print some information
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
      }
      
      # finalize the fcs file to write as output
      fcs <- homogenize_flowFrame(fcs_raw, ref_panel)
      
      # write output fcs file to the specified directory
      filename <- file.path(output_data_path, paste0(prefix, file))
      flowCore::write.FCS(fcs, filename)
      
    }
    
  }
