system.or.stop <- function
### Run a command line or stop with an error.
(cmd,
### Command line, passed to system.
  verbose=0
### print output?
){
  if(verbose)cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
### Nothing.
}
