#!/bin/bash
cd ..
rm -rf PeakSegPipeline-release
cp -r PeakSegPipeline PeakSegPipeline-release
rm PeakSegPipeline-release/tests/testthat/test-pipeline-*.R
grep -v Remotes PeakSegPipeline/DESCRIPTION > PeakSegPipeline-release/DESCRIPTION
PKG_TGZ=$(R CMD build PeakSegPipeline-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
