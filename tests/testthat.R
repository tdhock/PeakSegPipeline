library(testthat)

download.to <- function
(u, f, writeFun=if(grepl("bigWig", f))writeBin else writeLines){
  if(!file.exists(f)){
    require(httr)
    f.dir <- dirname(f)
    dir.create(f.dir, showWarnings=FALSE, recursive=TRUE)
    request <- GET(u)
    stop_for_status(request)
    writeFun(content(request), f)
  }
}

suite <- Sys.getenv("TEST_SUITE")
if(suite==""){
  PeakSegPipeline::system.or.stop("cd .. && bash build.sh")
}else{
  test_check("PeakSegPipeline", filter=suite)
}
