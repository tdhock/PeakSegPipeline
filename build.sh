#!/bin/bash
cd ..
rm -rf PeakSegPipeline-release
cp -r PeakSegPipeline PeakSegPipeline-release
grep -v Remotes PeakSegPipeline/DESCRIPTION > PeakSegPipeline-release/DESCRIPTION
rm PeakSegPipeline-release/tests/testthat/*
cp PeakSegPipeline/tests/testthat/test-CRAN*.R PeakSegPipeline-release/tests/testthat
PKG_TGZ=$(R CMD build PeakSegPipeline-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
