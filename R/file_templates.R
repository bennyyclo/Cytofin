#' Generate a template for a cytofin metadata file
#' 
#' `cytofin_generate_metadata_template` creates a template metadata .csv file 
#' (with the correct columns and dummy example data) in a specified location. 
#'
#' @param file_name A string representing the name of the .csv file to be 
#' saved in the directory specified by `template_path`. Defaults to 
#' "template_metadata.csv"
#'
#' @param template_path A file path or connection where the template file should be 
#' written. Defaults to the current working directory
#'
#' 
#' @export
#'
#' @examples
#' # specify the path where you'd like to store the template file
#' my_name <- "metadata_template.csv"
#' my_path <- file.path("~", "Desktop", "template_folder")
#'     
#' 
#' # generate the template file, which then can be edited manually 
#' cytofin_generate_metadata_template(
#'    file_name = my_name, 
#'    template_path = my_path
#' )
#' 
cytofin_generate_metadata_template <- 
  function(
    file_name = "template_metadata.csv", 
    template_path = getwd()
  ) { 
  
  #create output data.frame
  output_frame <- 
    data.frame(
      filename = c("file_1.fcs", "file_2.fcs", "file_3.fcs", "file_4.fcs"), 
      cohort = c("cohort_1", "cohort_1", "cohort_2", "cohort_2"),
      plate_number = c("plate_1", "plate_1", "plate_2", "plate_2"), 
      patient_id = c("patient_1", "patient_2", "patient_a", "patient_b"), 
      condition = c("basal", "basal", "stimulation_1", "stimulation_2"), 
      is_anchor = c(0, 1, 0, 1), 
      validation = 
        paste0(
          "validation_", 
          c("file_1.fcs", "file_2.fcs", "file_3.fcs", "file_4.fcs")
        )
    )
  
  readr::write_csv(
    x = output_frame, 
    file = file.path(template_path, file_name)
  )
  
}


#' Generate a template for a cytofin reference panel file
#' 
#' `cytofin_generate_panel_template` creates a template reference panel .csv file 
#' (with the correct columns and dummy example data) in a specified location. 
#'
#' @param file_name A string representing the name of the .csv file to be 
#' saved in the directory specified by `template_path`. Defaults to 
#' "template_panel_info.csv"
#'
#' @param template_path File path or connection where the template file should be 
#' written. Defaults to the current working directory
#'
#' 
#' @export
#'
#' @examples
#' 
#' # specify the path where you'd like to store the template file
#' my_name <- "panel_template.csv"
#' my_path <- file.path("~", "Desktop", "template_folder")
#' 
#' # generate the template file, which then can be edited manually 
#' cytofin_generate_panel_template(
#'    file_name = my_name, 
#'    template_path = my_path
#' )
#' 
cytofin_generate_panel_template <- 
  function(
    file_name = "template_panel_info.csv", 
    template_path = getwd()
  ) { 
    
  # TO DO: Give a more compelling example of markers and their regex
  #create output data.frame
  output_frame <- 
    data.frame(
      metal_name = c("Time", "Event_length", "(Pd102)Di", "(Pd104)Di"), 
      antigen_name = c("Time", "Event_length", "marker_name_1", "marker_name_2"),
      lineage = c(0, 0, 1, 1), 
      functional = c(0, 0, 0, 1), 
      general = c(0, 1, 1, 1)
    )
  
  readr::write_csv(
    x = output_frame, 
    file = file.path(template_path, file_name)
  )
  
}
