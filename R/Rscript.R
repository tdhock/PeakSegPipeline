Rscript <- function(...){
  code <- sprintf(...)
  if(any(grepl("'", code))){
    print(code)
    stop("there can not be any ' in code")
  }
  sprintf("Rscript -e '%s'", code)
}
