system.or.stop <- function
### Run a command line or stop with an error.
(cmd
### Command line, passed to system.
){
  cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
### Nothing.
}
