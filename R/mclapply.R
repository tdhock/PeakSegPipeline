### Run mclapply inside of a for loop, ensuring that it never receives
### a first argument with a length more than maxjobs. This avoids some
### memory problems (swapping, or getting jobs killed on the cluster)
### when using mclapply(1:N, FUN) where N is large.
maxjobs.mclapply <- function(X, FUN, maxjobs=getOption("mc.cores", 1L)){
  maxjobs <- if(is.numeric(maxjobs) && length(maxjobs)==1){
    as.integer(maxjobs)
  }else{
    1L
  }
  if(maxjobs == 1L)return(lapply(X, FUN))
  N <- length(X)
  i.list <- splitIndices(N, N/maxjobs)
  result.list <- list()
  for(i in seq_along(i.list)){
    i.vec <- i.list[[i]]
    gc()
    result.list[i.vec] <- mclapply(X[i.vec], FUN)
    gc()
  }
  result.list
}

### mclapply with error checking.
mclapply.or.stop <- function(...){
  result.list <- maxjobs.mclapply(...)
  is.error <- sapply(result.list, inherits, "try-error")
  if(any(is.error)){
    print(result.list[is.error])
    stop("errors in mclapply")
  }
  result.list
}

### Set mc.cores option from the PBS_NUM_PPN environment variable.
PPN.cores <- function(variable="PBS_NUM_PPN"){
  stopifnot(is.character(variable))
  stopifnot(length(variable) == 1)
  ppn.txt <- Sys.getenv(variable)
  ppn <- as.integer(ppn.txt)
  if(is.finite(ppn)){
    options(mc.cores=ppn)
  }
  ppn
}

