#' Ask User to Continue with Yes/No Input
#'
#' This function prompts the user with a message asking whether they want to continue, 
#' and accepts either 'Y' (for Yes) or 'N' (for No) as valid inputs.
#' The function will repeatedly ask the user for input until they enter a valid response.
#'
#' @param msg A character string containing the message to display to the user before the prompt.
#' @return A logical value: `TRUE` if the user inputs 'Y' (Yes) and `FALSE` if the user inputs 'N' (No).
#' @export
ask_user_continue <- function(msg) {
  repeat {
    user_input <- readline(prompt = sprintf("%s. Do you want to continue? (Y/N): ", msg))

    if (toupper(user_input) == "Y") {
      return(TRUE)
    } else if (toupper(user_input) == "N") {
      return(FALSE)
    } else {
      cat("Invalid input. Please enter 'Y' for Yes or 'N' for No.\n")
    }
  }
}





# # Function to compress a CSV file using bzip2
# #' @importFrom utils write.csv
# compress_csv_bzip2 <- function(input_csv) {
#   output_csv_bz <- paste0(input_csv, ".bz2")
#   data <- read.csv(input_csv, check.names = FALSE)
#   write.csv(data, bzfile(output_csv_bz), row.names = FALSE)
# }

# Decompress and load a CSV file using bzip2
#' @importFrom utils read.csv
read.csv_bzip2 <- function(input_csv_bz, ...) {
  data <- read.csv(bzfile(input_csv_bz), ...)
  return(data)
}

# # # Directory containing your CSV files
# csv_directory <- "E:/xhmhc/BrainEnrich/inst/extdata/not_compressed"
# # List all CSV files in the directory
# csv_files <- list.files(path = csv_directory, pattern = "\\.csv$", full.names = TRUE)
# # Compress each CSV file using bzip2
# lapply(csv_files, compress_csv_bzip2)

# test1=read.csv_bzip2("E:/xhmhc/BrainEnrich/inst/extdata/not_compressed/desikan_r0.6.csv.bz2", check.names=FALSE)
# test2=read.csv("E:/xhmhc/BrainEnrich/inst/extdata/not_compressed/desikan_r0.6.csv",check.names=FALSE)

# # Compare column names
# if (!identical(names(test1), names(test2))) {
#   cat("Column names are different:\n")
#   print(setdiff(names(test1), names(test2)))
#   print(setdiff(names(test2), names(test1)))
# }

# # Compare dimensions
# if (!identical(dim(test1), dim(test2))) {
#   cat("Dimensions are different:\n")
#   print(dim(test1))
#   print(dim(test2))
# }

# # Compare data frames using all.equal
# differences <- all.equal(test1, test2)

# # Print the differences
# if (isTRUE(differences)) {
#   cat("Data frames are identical.\n")
# } else {
#   cat("Differences found:\n")
#   print(differences)
# }
