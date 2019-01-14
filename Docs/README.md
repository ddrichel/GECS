
# Genomic Exhaustive Collapsing Scan (GECS)
 ###### _An exhaustive genomic scan for association in genetic data_

### Table of contents
* [Installation](#Instalation)
* [Description](#Description)
* [Components](#Components)
  * [Single marker analysis](#Single_Marker_Analysis)
  * [Variable binning](#Variable_Binning)
* [Getting started](#Getting_started)
  * [Examples](#Examples)
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
We developed a new approach to conduct association analysis for rare variants exhaustively in whole-genome or whole-exome data sets, by variating bins sizes and MAF tresholds. GECS is an ultra fast program that perform an exhaustive scan for association in case-control genetic data, implemented in c++.

##### Prerequisites

GECS is distributed under GPL3 license. Starting from GECS 1.1.1, it supports c++ (?) on linux systems.
This program uses the alglib c++ library.

##### Usage

The usage of gecs is very simple. Only execute the follwoing command.
gecs <~/path/to/file.param>

Keywords in the parameter file [\*.param](https://github.com/ddrichel/GECS/tree/master/Docs/DATA/example_1.param) :

**BFILE** _\<string\>_         (prefix of the plink binary file) 

**SINGLEMARKER**	_\<bool\>_		  (whether single-marker analysis should be performed instead of VB (default=0))   

**PERMUTATIONS**	_\<int\>_		   (number of permutations for correction of multiple testing, default=0) 

**NCT**		_\<int\>_		           ("rareness" threshold: max. number of carriers per variant)

**MAFT** _\<double\>_          (Minor allele fequency threshold for rare variants) 

**PTHRESHOLD**	_\<double\>_		  (max. nominal p-value for bins to be written to output files, optional)

**ALLBINS**		_\<bool\>_		      (whether locally not-distinct bins should be written to output (useful for plotting, default=0))

**OR**		_\<bool\>_		           (Whether odd ratios will be calculated, optional)

**CORRECTED_P** _\<bool\>_     (Whether calculated p values will be corrected by wilson score interval of CI 95%, optional) 

**OUTPUT**		_\<string\>_ 		    (Prefix of output files)

### Components

GECS provides two major features for conducting association analysis for rare variants, namely for signle markers and for all possible bins (subsequences of contiguous markers) in the genetic data set. In the single marker analysis, all variants will be considered in the analysis, regardless of their frequencies. however, in the variable binning approach we need to specify a threshold of minor allele frequency. Permutations with respect to the case-control labels is applied to make correction for multiple testing. That means if PERMUTATIONS==0, then there is no correction for multiple testing will be done. Moreover, yoe have the possibility to get the corrected p values by wilson score interval for conficence interval of 95%. You have the option to calcutae the odds ratios for all bins by specifying OR==1. (OR=0 is by default)

#### Sigle_Marker_Analysis
If SINGLEMARKER==1, then GECS will conduct only the single-marker test on all variants included in the analysis.
-NOTICE- that in case of SINGLEMARKER==1, MAFT will be set to the maximum and the NCT and ALLBINS parameters will be ignored.
Input and output string parameter are to cpecify only the name of input and output file without any extensions.
###### Resulted files:
1. A text file \[\*_singlemarker.txt\], which includes the results of association tests for all included markers.
2. Two files \[\*_singlemarker.pvals\] and \[\*_singlemarker.srt.pvals\] for the correction of multiple testing based on the null model.
-NOTICE- the final corrected alpha will be reported with other informations about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.

#### Variable_Binning
In the default case (SINGLEMARKER==0) the parameters NCT and MAFT do the same job, which determinig the rareness threshold for the analysis.
-NOTICE- that sypecifying NCT will overwrite the parameter MAFT, otherwise, specifying only MAFT will be enough.
ALLBINS paremeter is set to 0 for the default case, where only distinct bins will be considered. This feature is only usefull in case of scanning small regions for plotting purposes.
Only 
###### Resulted files:
1. A text file \[\*_gecs_\<nct\>.txt\], which includes the results of association tests for all genomic subsequences.
2. Two files \[\*_gecs_\<nct\>.pvals\] and \[\*_gecs_\<nct\>.srt.pvals\] for the correction of multiple testing based on the null model.

-NOTICE- the final corrected alpha will be reported with other informations about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.

### Getting_started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

#### Examples
execute some examples of gecs included in the bash [file](https://github.com/ddrichel/GECS/blob/master/runSampleCode.sh), applied on an example data set. 
In the [Examples](https://github.com/ddrichel/GECS/tree/master/Docs/DATA) you will try different combinations of parameter specifications for a specific goal in the anaylsis.
1) Example 1 : You will perform a single marker analysis and use 999 permutations for correcting for multiple testing. Note no need to specify the minor allele frequency threshold (MAFT) or the threshold of number of cariers (NCT).
2) Example 2: You will perform the variable binning approach instead of the single marker analysis, and here you are obligated to specify the parameter MAFT or the equivalent parameter NCT like in the example 3.
3) Example 4: In this example a conflict in the values of NCT and MAFT will be solved by cosidering NCT always it is specified.
4) Example 5: The SMA approach is applied on the chromosome 22 of a simulated sata set with 1000 sample size.

### Attributions

#### Authors
Dmitriy Drichel -initial work-
George Kanoungi

#### License
This project is licensed under the MIT License - see the LICENSE.md file for details

#### Acknowledgment
This project was supported by the German Research Foundation grant BE 38/28/9-1. The funding organization did not have any influence on the design, conduct or conclusions of the study.

#### Citation
Please site: George Kanoungi, Michael Nothnagel and Dmitriy Drichel. An exhaustive genomic scan (2019)


