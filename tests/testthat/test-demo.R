library(testthat)
context("demo")

## Download bigWig files from github.
bigWig.part.vec <- c(
  "Input/MS010302",
  "bcell/MS010302",
  "Input/MS002202",
  "kidney/MS002202",
  "Input/MS026601",
  "bcell/MS026601",
  "Input/MS002201",
  "kidney/MS002201"
    )
label.txt <- "
chr10:33,061,897-33,162,814 noPeaks
chr10:33,456,000-33,484,755 peakStart kidney
chr10:33,597,317-33,635,209 peakEnd kidney
chr10:33,662,034-33,974,942 noPeaks

chr10:35,182,820-35,261,001 noPeaks
chr10:35,261,418-35,314,654 peakStart bcell kidney
chr10:35,343,031-35,398,459 peakEnd bcell kidney

chr10:38,041,023-38,102,554 noPeaks
chr10:38,296,008-38,307,179 peakStart bcell kidney
chr10:38,379,045-38,391,967 peakStart bcell kidney
chr10:38,404,899-38,412,089 peakEnd bcell kidney
chr10:38,413,073-38,444,133 noPeaks

chr10:38,585,584-38,643,190 noPeaks
chr10:38,643,191-38,650,766 peakStart bcell kidney
chr10:38,731,066-38,750,574 peakEnd bcell kidney
chr10:38,750,960-38,790,663 noPeaks

chr10:38,807,475-38,815,200 noPeaks
chr10:38,815,201-38,816,355 peakStart bcell kidney Input
chr10:38,818,377-38,819,342 peakEnd bcell kidney Input

chr10:39,098,319-39,111,384 noPeaks
chr10:39,125,134-39,125,550 peakStart bcell kidney Input
chr10:39,125,594-39,126,266 peakEnd bcell kidney Input
chr10:39,126,866-39,140,858 noPeaks
"
set.dir <- file.path(Sys.getenv("HOME"), "PeakSegPipeline-test", "demo")
repos.url <- "https://raw.githubusercontent.com/tdhock/input-test-data/master/"
for(bigWig.part in bigWig.part.vec){
  suffix <- ifelse(grepl("MS026601|MS002201", bigWig.part), "/", "_/")
  bigWig.file <- file.path(set.dir, "samples", sub("/", suffix, bigWig.part), "coverage.bigWig")
  bigWig.url <- paste0(repos.url, bigWig.part, ".bigwig")
  download.to(bigWig.url, bigWig.file)
}
labels.file <- file.path(set.dir, "labels", "some_labels.txt")
dir.create(dirname(labels.file), showWarnings=FALSE, recursive=TRUE)
writeLines(label.txt, labels.file)
problems.bed <- file.path(set.dir, "problems.bed")
unlink(problems.bed)
cat("chr10	60000	17974675
chr10	18024675	38818835
chr10	38868835	39154935
chr10	42746000	46426964
chr10	47529169	47792476
chr10	47892476	48055707
chr10	48105707	49095536
chr10	49195536	51137410
chr10	51187410	51398845
chr10	51448845	125869472
chr10	125919472	128616069
", file=problems.bed)

## Whole pipeline.
system(paste("bigWigToBedGraph", bigWig.file, "/dev/stdout|head"))
pipeline(set.dir)
index.html <- file.path(set.dir, "index.html")
test_that("index.html is created", {
  expect_true(file.exists(index.html))
})

## Post-processing to explain the output.
joint.problems.dt <- fread(paste("cat", file.path(
  set.dir, "problems", "*", "jointProblems.bed")))
setnames(joint.problems.dt, c("chrom", "problemStart", "problemEnd"))
peaks.glob <- file.path(
  set.dir, "problems", "*", "jointProblems", "*", "peaks.bed")
##Sys.glob(peaks.glob)
joint.peaks.dt <- fread(paste("cat", peaks.glob))
setnames(
  joint.peaks.dt,
  c("chrom", "peakStart", "peakEnd", "sample.path", "mean"))
joint.peaks.dt[, sample.id := sub(".*/", "", sample.path)]
joint.peaks.dt[, sample.group := sub("/.*", "", sample.path)]
chunks.dt <- fread(paste("cat", file.path(
  set.dir, "problems", "*", "chunks", "*", "chunk.bed")))
setnames(chunks.dt, c("chrom", "chunkStart", "chunkEnd"))
all.problems <- fread(file.path(set.dir, "problems.bed"))
setnames(all.problems, c("chrom", "problemStart", "problemEnd"))
data.start <- min(chunks.dt$chunkStart)
data.end <- max(chunks.dt$chunkEnd)
two.problems <- all.problems[!(
  problemEnd < data.start |
    data.end < problemStart),]
both.problems <- rbind(
  data.table(two.problems, y="separate problems"),
  data.table(joint.problems.dt, y="joint problems"))

limits.list <- list(
  c(34, 35),
  c(35, 36),
  c(36, 37),
  c(37, 38),
  c(38, 39))
limits.dt <- data.table(matrix(unlist(limits.list)*1e6, ncol=2, byrow=TRUE))
setnames(limits.dt, c("plotStart", "plotEnd"))
limits.dt$y <- "plots in un-labeled regions"

labels.bed.vec <- Sys.glob(file.path(
  set.dir, "samples", "*", "*", "labels.bed"))
all.labels.list <- list()
for(labels.bed in labels.bed.vec){
  sample.dir <- dirname(labels.bed)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  sample.group <- basename(group.dir)
  sample.labels <- fread(labels.bed)
  setnames(sample.labels, c("chrom", "labelStart", "labelEnd", "annotation"))
  all.labels.list[[labels.bed]] <- data.table(
    sample.id, sample.group, sample.labels)
}
all.labels <- do.call(rbind, all.labels.list)

input.labels <- all.labels[sample.group=="Input", list(
                             prop.noPeaks=mean(annotation=="noPeaks")
  ), by=.(labelStart, labelEnd)]
setkey(input.labels, labelStart, labelEnd)
input.pred <- joint.peaks.dt[, list(
  n.Input=sum(sample.group=="Input")
  ), by=.(peakStart, peakEnd)]
setkey(input.pred, peakStart, peakEnd)
labeled.input <- foverlaps(input.pred, input.labels, nomatch=0L)
thresh.dt <- labeled.input[, data.table(WeightedROC(
  n.Input, ifelse(prop.noPeaks==0, 1, -1)))]
thresh.best <- thresh.dt[which.min(FP+FN),]
## threshold is smallest n.Input that is classified as non-specific.
setkey(joint.peaks.dt, peakStart, peakEnd)
peaks.with.counts <- input.pred[joint.peaks.dt]
peaks.with.counts[, specificity := ifelse(
                                  n.Input >= thresh.best$threshold, "non-specific", "specific")]

gg <- ggplot()+
  coord_cartesian(xlim=c(data.start, data.end)/1e3, expand=TRUE)+
  geom_segment(aes(
    problemStart/1e3, y,
    xend=problemEnd/1e3, yend=y),
               color="blue",
               data=both.problems)+
  geom_point(aes(
    problemStart/1e3, y),
             color="blue",
             data=both.problems)+
  geom_tallrect(aes(
    xmin=chunkStart/1e3, xmax=chunkEnd/1e3),
                alpha=0.1,
                data=chunks.dt)+
  geom_segment(aes(
    plotStart/1e3, y,
    xend=plotEnd/1e3, yend=y),
               data=limits.dt)+
  geom_point(aes(
    plotStart/1e3, y),
             data=limits.dt)+
  geom_point(aes(
    peakStart/1e3, sample.path, color=specificity),
             data=peaks.with.counts)+
  geom_segment(aes(
    peakStart/1e3, sample.path,
    xend=peakEnd/1e3, yend=sample.path, color=specificity),
               data=peaks.with.counts)+
  ylab("")+
  xlab("position on chr10 (kb = kilo bases)")+
  scale_color_manual(values=c("non-specific"="red", specific="deepskyblue"))
png(file.path(set.dir, "figure-demo-overview.png"), res=100, h=200, w=1000)
print(gg)
dev.off()

coverage.bigWig.vec <- Sys.glob(file.path(
  set.dir, "samples", "*", "*", "coverage.bigWig"))
ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
img.list <- list()
for(limit.i in seq_along(limits.list)){
  limit.vec <- limits.list[[limit.i]]*1e6
  print(limit.vec)
  coverage.list <- list()
  for(coverage.bigWig in coverage.bigWig.vec){
    sample.dir <- dirname(coverage.bigWig)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    cmd <- sprintf(
      "bigWigToBedGraph -chrom=chr10 -start=%d -end=%d %s /dev/stdout",
      limit.vec[1], limit.vec[2],
      coverage.bigWig)
    sample.coverage <- fread(cmd)
    setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
    coverage.list[[coverage.bigWig]] <- data.table(
      sample.id, sample.group, sample.coverage)
  }
  coverage <- do.call(rbind, coverage.list)

  show.peaks <- joint.peaks.dt[
                               !(peakEnd < limit.vec[1] | limit.vec[2] < peakStart),]
  show.labels <- all.labels[
                            !(labelEnd < limit.vec[1] | limit.vec[2] < labelStart),]

  gg <- ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(sample.group + sample.id ~ ., scales="free")+
    geom_tallrect(aes(
      xmin=labelStart/1e3, 
      xmax=labelEnd/1e3,
      fill=annotation), 
                  alpha=0.5,
                  data=show.labels)+
    scale_fill_manual(values=ann.colors)+
    scale_y_continuous(
      "aligned read coverage",
      breaks=function(limits){
        lim <- floor(limits[2])
        if(lim==0){
          Inf
        }else{
          lim
        }
      })+
    scale_x_continuous(paste(
      "position on",
      coverage$chrom[1],
      "(kb = kilo bases)"))+
    geom_step(aes(
      chromStart/1e3, count),
              color="grey50",
              data=coverage)+
    coord_cartesian(xlim=limit.vec/1e3, expand=FALSE)
  if(nrow(show.peaks)){
    gg <- gg+
      geom_point(aes(
        peakStart/1e3, 0),
                 color="deepskyblue",
                 size=3,
                 data=show.peaks)+
      geom_segment(aes(
        peakStart/1e3, 0,
        xend=peakEnd/1e3, yend=0),
                   color="deepskyblue",
                   size=3,
                   data=show.peaks)
  }
  f <- sprintf(
    "%s/figure-demo-%d-%d.png",
    set.dir,
    limit.vec[1]/1e6, limit.vec[2]/1e6)
  png(f, res=100, width=1000, height=600)
  print(gg)
  dev.off()

  img.list[[limit.i]] <- sprintf(
    '<p><img src="%s" /></p>',
    basename(f))
}

index.lines <- readLines(index.html)
new.index.lines <- c('
<h2>Explanation of results below</h2>
<ul>
  <li>
    PeakSegFPOP + PeakSegJoint were trained using four labeled samples
    (two H3K36me3 samples, two Input samples). These labeled samples
    have colored rectangular labels/annotations in the plots below,
    which indicate genomic regions where a scientist has visually
    determined presence or absence of peaks.
  </li>
  <li>
    The black peak predictions represent the PeakSegFPOP model, which
    predicts peaks separately for each sample. These peaks do not
    occur in the exact same positions across samples, but are clustered
    together to define "joint problems" where the PeakSegJoint model
    is fit to all samples.
  </li>
  <li>
    The blue peak predictions represent the PeakSegJoint model, which
    predicts peaks separately for each genomic region, but jointly
    across all samples. These peaks occur in the exact same positions
    across samples, and are reported in the final output files
(<a href="peaks_summary.tsv">peaks_summary.tsv</a>,
   <a href="peaks_matrix.tsv">peaks_matrix.tsv</a>).
  </li>
  <li>
    The models were used to predict peaks in eight samples,
    four of which have no
    labels. The un-labeled samples appear without any colored
    rectangles in the plots below, and can be used for evaluating the
    accuracy of the learned model. For example, consider the genomic region below:</li>
<img src="problems/chr10:18024675-38818835/chunks/chr10:33061897-33974942/figure-predictions.png" />
  <li>On the left, the model predicted a common peak in all bcell and kidney samples
 (even though this peak was not labeled in any of the samples).</li>
  <li>In the middle, the model predicted a common peak in the two kidney samples,
  and no corresponding peak in the two bcell samples
(only one of the two bcell samples was labeled as having noPeaks,
 and only one of the two kidney samples was labeled as having a peakStart and peakEnd).</li>
</ul>
', index.lines, '
<h2>More plots</h2>
<p>The plots below are usually not included in the output,
but we added them here, so you can see peak predictions
at all genomic regions that have coverage in this data set.</p>
<img src="figure-demo-overview.png" />
', unlist(img.list))
writeLines(new.index.lines, index.html)
browseURL(index.html)
