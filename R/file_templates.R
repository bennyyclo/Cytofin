#' Generate a template for a cytofin metadata file
#' 
#' `cytofin_generate_metadata_template` creates a template metadata .csv file 
#' (with the correct columns and dummy example data) in a specified location. 
#'
#' @param template_path File path or connection where the template file should be 
#' written. Defaults to the current working directory
#'
#' 
#' @export
#'
#' @examples
#' # specify the path where you'd like to store the template file
#' my_path <- 
#'     file.path("~", "Desktop", "template_folder", 
#'               "metadata_template.csv")
#' 
#' # generate the template file, which then can be edited manually 
#' cytofin_generate_metadata_template(template_path = my_path)
#' 
cytofin_generate_metadata_template <- function(template_path = getwd()) { 
  
  #create output data.frame
  output_frame <- 
    data.frame(
      filename = c("file_1.fcs", "file_2.fcs", "file_3.fcs", "file_4.fcs"), 
      cohort = c("cohort_1", "cohort_1", "cohort_2", "cohort_2"),
      plate_number = c("plate_1", "plate_2", "plate_3", "plate_4"), 
      patient_id = c("patient_1", "patient_2", "patient_a", "patient_b"), 
      condition = c("basal", "basal", "stimulation_1", "stimulation_2"), 
      population = c("B-cells", "B-cells", "B-cells", "T-cells"), 
      validation = 
        paste0(
          "homogenized_", 
          c("file_1.fcs", "file_2.fcs", "file_3.fcs", "file_4.fcs")
        )
    )
  
  readr::write_csv(
    x = output_frame, 
    file = template_path
  )
  
}



#' Generate a template for a cytofin reference panel file
#' 
#' `cytofin_generate_panel_template` creates a template reference panel .csv file 
#' (with the correct columns and dummy example data) in a specified location. 
#'
#' @param template_path File path or connection where the template file should be 
#' written. Defaults to the current working directory
#'
#' 
#' @export
#'
#' @examples
#' # specify the path where you'd like to store the template file
#' my_path <- 
#'     file.path("~", "Desktop", "template_folder", 
#'               "panel_template.csv")
#' 
#' # generate the template file, which then can be edited manually 
#' cytofin_generate_panel_template(template_path = my_path)
#' 
cytofin_generate_panel_template <- function(template_path = getwd()) { 
  
  # TO DO: Give a more compelling example of markers and their regex
  #create output data.frame
  output_frame <- 
    data.frame(
      desc = c("Time", "Event_length", "(Pd102)Di", "(Pd104)Di"), 
      range = c("Time", "Event_length", "marker_name_1", "marker_name_2"),
      metal_pattern = c("[Tt]ime", "ength", "Pd102", "Pd104"), 
      antigen_pattern = c("[Tt]ime", "ength", "regex_1", "regex_2"), 
      Lineage = c(0, 0, 1, 1), 
      Functional = c(0, 0, 0, 1), 
      General = c(0, 1, 1, 1)
    )
  
  readr::write_csv(
    x = output_frame, 
    file = template_path
  )
  
}
