#!/usr/bin/env bash

# $File: runSampleCode.sh $
# $Created: 13.09.2018 14:56:13 $
#
# This file is part of GECS, an application to conduct association analysis for rare variants exhaustively in whole-genome or whole-exome data sets.
#
# Copyright George Kanoungi, Dmitriy Drichel (2019) (george.kanoungi@uni-koeln.de, dmitriy@drichel-analytics.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# This bash script apply GECS on the example data set provided in GECS/Docs/DATA in binary plink format:

# Input files for GECS are binary plink formats, i.e. *.bed, *.bim, and *.fam.
# Input parameters and Output file name has to specified in the paremeter file (*.sfile):
#
# In addition to the input and output file names, basic parameters have to be specified:
 
# For SINGLEMARKER = 1 (The case of conducting the standard single-marker test ONLY):
# The file *singlemarker.txt will be generated, which includes the results of association tests.

# For SINGLEMARKER = 0 (The case of conducting GECS with variable binning ONLY):
# The parameter MAFT or NCT have to be specified. If both parameters were specified, then MAFT will be ignored.
# The file *_gecs_nct_*.txt will be generated, which includes the results of the exhaustive scan for association.

# In both cases, two additional files will be generated, namely *.pvals and *.srt_pvals, and include the the results of the correction for multiple testing.
#

# Uncomment one of the following lines and then execute the bash file
# Use for that the command "bash runSampleCode.sh" 


./gecs ./DATA/example_1.param
#./gecs ./DATA/example_2.param
#./gecs ./DATA/example_3.param
#./gecs ./DATA/example_4.param
#./gecs ./DATA/example_5.param