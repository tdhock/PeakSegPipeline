Rscript <- function
### Run sprintf and paste Rscript -e before -- useful for constructing
### command lines.
(...
### Passed to sprintf -- there should be no single quotes.
){
  code <- sprintf(...)
  if(any(grepl("'", code))){
    print(code)
    stop("there can not be any ' in code")
  }
  sprintf("Rscript -e '%s'", code)
### Command line.
}
