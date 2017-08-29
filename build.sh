#!/bin/bash
cd ..
rm -rf PeakSegPipeline-release
cp -r PeakSegPipeline PeakSegPipeline-release
##grep -v PeakSegPipelineData PeakSegPipeline/DESCRIPTION | grep -v Remotes > PeakSegPipeline-release/DESCRIPTION
##rm PeakSegPipeline-release/tests/testthat/test-PeakSegPipelineData.R
PKG_TGZ=$(R CMD build PeakSegPipeline-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
