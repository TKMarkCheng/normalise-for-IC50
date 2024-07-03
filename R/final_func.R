#' Strings all functions together.
#'
#' @param neut_xlsx_path a file to normalise, where by default negative control is left most column, and positive control is second left column.
#' @param sheetname the sheet of the .xlsx file to read, default for promega excel file is "Results"
#' @param control_neg_column The column which negative/cell control are located (column number from left, integer,default=1)
#' @param control_pos_column The column which positive/cell control are located (column number from left, integer,default=2)
#' @param rotation_deg_needed degrees that dataframe needs to be rotated, default=0. Feeds into rotate.R
#' @param cluster if the k-clusters dependent anomaly detection is invoked, default=FALSE. Feed into normalise.R
#' @return A rotated, normalised dataframe where negative and positive control columns are removed.
#' @export

#final function
final_func <- function(neut_xlsx_path, 
                       control_neg_column=c(1), control_pos_column=c(2),
                       control_plating = "column",
                       control_neg_wells=c("A1","B1","C1","D1","E1","F1","G1","H1"),
                       control_pos_wells=c("A2","B2","C2","D2","E2","F2","G2","H2"),
                       sheetname="Results", rotation_deg_needed=0, cluster=FALSE){
  message("processing input file: ",neut_xlsx_path)
  input_name = basename(neut_xlsx_path)
  # 1. read in and select relevant 8x12 area of the excel Results sheet
  df <- read_promega_plate_excel(neut_xlsx_path)
  # 2. rotate by n degrees, if needed
  df <- rotate(df,rotation_deg_needed = rotation_deg_needed)
  # 3. normalize the rotated df
  output <- normalise(df = df, control_plating = control_plating,
                      control_neg_column = control_neg_column,control_pos_column = control_pos_column,
                      control_neg_wells = control_neg_wells, control_pos_wells = control_pos_wells, 
                      cluster = cluster,
                      filename = input_name)
  # 4. Return and announce completion
  message("Finished processing: ",neut_xlsx_path)
  return(output)
}
