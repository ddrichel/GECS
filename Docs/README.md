
# Genomic Exhaustive Collapsing Scan (GECS)
 ###### _An exhaustive genomic scan for association in genetic data_

### Table of contents
* [Installation](#Instalation)
* [Description](#Description)
* [Components](#Components)
  * [Single marker analysis](#Single_Marker_Analysis)
  * [Variable binning](#Variable_Binning)
* [Getting started](#Getting_started)
  * [Example 1](#Example_1)
  * [Example 2](#Example_2)
  * [Example 3](Example_3)
  * [Example 4](Example_4)
* [Attributions](#Attributions)
  * [Authors](#Authors)
  * [Aknowledgment](#Acknowledgment)
  * [Citation](#Citation)

### Instalation

git clone --recursive https://github.com/ddrichel/GECS
cd GECS
./configure?
make

### Description 
We developed a new approach to conduct association analysis for rare variants exhaustively in whole-genome or whole-exome data sets, by variating bins sizes and MAF tresholds. GECS is an ultra fast program that perform an exhaustive scan for association in case-control genetic data.

##### Prerequisites

GECS is distributed under GPL3 license. Starting from GECS 1.1.1, it supports c++ (?) on linux systems.
This program uses the alglib c++ library.

##### Usage

The usage of gecs is very simple. Only execute the follwoing command.
gecs <~/path/to/file.param>

Keywords in the parameter file [*.param](https://github.com/ddrichel/GECS/tree/master/Docs/DATA/example_1.param) :

**BFILE** _\<string\>_         (prefix of the plink binary file) 

**SINGLEMARKER**	_\<bool\>_		  (whether single-marker analysis should be performed instead of VB (default=0))   

**PERMUTATIONS**	_\<int\>_		   (number of permutations) 

**NCT**		_\<int\>_		           ("rareness" threshold: max. number of carriers per variant)

**MAFT** _\<double\>_          (Minor allele fequency threshold for rare variants) 

**PTHRESHOLD**	_\<double\>_		  (max. nominal p-value for bins to be written to output files)

**ALLBINS**		_\<bool\>_		      (whether locally not-distinct bins should be written to output (useful for plotting, default=0))

**OR**		_\<bool\>_		           (Whether odd ratios will be calculated)

**CORRECTED_P** _\<bool\>_     (Whether calculated p values will be corrected by wilson score interval of CI 95%) 

**OUTPUT**		_\<string\>_ 		    (Prefix of output files)

### Components

GECS is provided by two major features to conduct the analysis for association, namely for signle markers and for all subsequences of contiguous markers. Permutations will be used to make correction for multiple testing. That means if PERMUTATIONS==0, then there is no correction for multiple testing will be done. Moreover, yoe have the possibility to get the corrected p values by wilson score interval for conficence interval of 95%. You have the option to calcutae the odds ratios for all bins by specifying OR==1. (OR=0 is by default)

#### Sigle_Marker_Analysis
If SINGLEMARKER==1, then GECS will conduct only the single-marker test on all variants included in the analysis.
-NOTICE- that in case of SINGLEMARKER==1, MAFT will be set to the maximum and the NCT and ALLBINS parameters will be ignored.
Input and output string parameter are to cpecify only the name of input and output file without any extensions.
###### Resulted files:
1. A text file \[*_singlemarker.txt\], which includes the results of association tests for all included markers.
2. Two files \[*_singlemarker.pvals\] and \[*_singlemarker.srt.pvals\] for the correction of multiple testing based on the null model.
#### Variable_Binning
In the default case (SINGLEMARKER==0) the parameters NCT and MAFT do the same job, which determinig the rareness threshold for the analysis.
-NOTICE- that sypecifying NCT will overwrite the parameter MAFT, otherwise, specifying only MAFT will be enough.
ALLBINS paremeter is set to 0 for the default case, where only distinct bins will be considered. This feature is only usefull in case of scanning small regions for plotting purposes.
###### Resulted files:
1. A text file \[*_gecs_\<nct\>.txt\], which includes the results of association tests for all genomic subsequences.
2. Two files \[*_gecs_\<nct\>.pvals\] and \[*_gecs_\<nct\>.srt.pvals\] for the correction of multiple testing based on the null model.
-NOTICE- the final corrected alpha will be reported with other informations about the analysis in the \[*_gecs_nct_\<nct\>.log\] file.

### Getting_started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

#### Example_1
#### Example_2
#### Example_3
#### Example_4

### Attributions

#### Authors

Dmitriy Drichel -initial work-

#### License

This project is licensed under the MIT License - see the LICENSE.md file for details

#### Acknowledgment
This project is financed by DFG.

#### Citation
Please site: George Kanoungi, Michael Nothnagel and Dmitriy Drichel. An exhaustive genomic scan (2019)


