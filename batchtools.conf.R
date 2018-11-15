cluster.functions = makeClusterFunctionsSlurm(system.file(
  file.path("templates", "slurm-afterok.tmpl"),
  package="PeakSegPipeline",
  mustWork=TRUE))
