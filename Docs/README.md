
# Genomic Exhaustive Collapsing Scan (GECS)
 ###### _An exhaustive genomic scan for association in genomic data_

### Table of contents
* [Installation](#Installation)
* [Description](#Description)
* [Getting started](#GS)
  * [Examples](#Examples)
* [Functionality](#Functionality)
  * [Single-marker analysis](#SMA)
  * [Variable binning](#VB)
  * [Correction for multple testing](#CMT)
* [Attributions](#Attributions)
  * [Authors](#Authors)
  * [License](#License)
  * [Acknowledgment](#Acknowledgment)
  * [Citation](#Citation)

### Installation

git clone https://github.com/ddrichel/GECS.git
cd GECS
make

### Description 
GECS is an implementation of our novel approach to conduct exhaustive region-based association analysis of rare variants in genomic case-control studies. The main idea of the exhaustive scan is to compute test statistics for all contigous subsequences of larger sequences, such as human chromosomes. The idea of the exhaustive one-dimensional scan can be regarded as analogous to Kulldorff's two-dimensional spacial scan (Kulldorff, 1996).

GECS is written in C++ and implements an ultrafast algorithm that enables testing for association of all genomic subsequences using the collapsing test of Li and Leal (2008). The application is scalable suitable for running analyses on large whole-exome and whole-genome studies.

##### Prerequisites

GECS is written in C++ and distributed under the GPL3 license.
This program uses the alglib C++ library for computing the p-values for the Pearson's chi-squared test (http://www.alglib.net/).

<a name="GS"/>

### Getting started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

#### Examples
execute some examples of gecs included in the bash [file](https://github.com/ddrichel/GECS/blob/master/runSampleCode.sh), applied on an example data set. 
In the [Examples](https://github.com/ddrichel/GECS/tree/master/Docs/DATA) you will try different combinations of parameter specifications for a specific goal in the analysis.
1) Example 1 : You will perform a single marker analysis and use 999 permutations for correcting for multiple testing. Note no need to specify the minor allele frequency threshold (MAFT) or the threshold of number of carriers (NCT).
2) Example 2: You will perform the variable binning approach instead of the single marker analysis, and here you are obligated to specify the parameter MAFT or the equivalent parameter NCT like in the example 3.
3) Example 4: In this example a conflict in the values of NCT and MAFT will be solved by considering NCT always it is specified.
4) Example 5: The SMA approach is applied on the chromosome 22 of a simulated data set with 1000 sample size.



##### Usage

The usage of gecs is:

**gecs** <~/path/to/**file.param**>

**gecs --help** will present a brief description of the keywords in the parameter file.

Keywords in the parameter file [\*.param](https://github.com/ddrichel/GECS/tree/master/Docs/DATA/example_1.param) :

| KEY | _\<Input\>_ | Default | Description | Notes |
|----|---|-------|------|------|
|**BFILE**|_\<string\>_|none| prefix of the plink binary file|With this parmater the name of the binary plink files \*.bed, \*.bim, and \*.fam of the data set will be specified.| 
|**SINGLEMARKER**|	_\<bool\>_|0|whether single-marker analysis should be performed instead of VB|If SINGLEMARKER=1, then the standard single marker analysis (SMA) will be performed and all variants will be included in the analysis regardless of their minor allele frequencies. For SINGLEMARKER=0, the variable binning (VB) approach will be preformed, where threshold of the minor allele frequency of variants included in the analysis need to be specified.|
|**NCT**|_\<int\>_|none|"rareness" threshold: max. number of carriers per variant|If SINGLEMARKER=1, then this parameter will be set to the maximum. Otherwise NCT or MAFT has to be specified|
|**MAFT**|_\<double\>_|none|Minor allele frequency threshold for rare variants|If SINGLEMARKER=1, then MAFT will be set on 0.5. Otherwise NCT or MAFT has to be specified| 
|**PERMUTATIONS**|_\<int\>_|999|Number of permutations for correction of multiple testing|In order to control the familiy wise error rate (FWER) at the 5% level, gecs performs in default 999 permutations with respect to the case-control labels. In each permutation, we will obtain the smallest p-value calculated in the analysis, so at the end we will have a list of 999 p values.|
|**CORRECTED_P**|_\<bool\>_|0|Whether calculated p values will be corrected by wilson score interval of CI 95%|Using the results of the FWER approach we can correct the results of the analysis and calculate the upper and lower limits of the confidence interval of 95% using the Wilson score interval.|
|**PTHRESHOLD**|_\<double\>_|1|Max. nominal p-value for bins to be written to output files|By specifying this parameter you can restrict the results on only p-values less than a specific threshold. This feature is usefull in case of having a limited memory size.|
|**ALLBINS**|_\<bool\>_|0|whether locally not-distinct bins should be written to output (useful for plotting)||
|**OR**|_\<bool\>_|0|Whether odd ratios will be calculated||
|**OUTPUT**|_\<string\>_|prefix of the plink binary file of the analyzed data|Prefix of output files||

### Functionality

GECS provides two major features for conducting association analysis for rare variants, namely for single markers and for all possible bins (subsequences of contiguous markers) in the genetic data set. In the single marker analysis, all variants will be considered in the analysis, regardless of their frequencies. However, in the variable binning approach we need to specify a threshold of minor allele frequency. Permutations with respect to the case-control labels is applied to make correction for multiple testing. That means if PERMUTATIONS=0, then there is no correction for multiple testing will be done. Moreover, you have the possibility to get the corrected p-values by Wilson score interval for confidence interval of 95%. You have the option to calculate the odds ratios for all bins by specifying OR 1. (OR=0 is by default)

<a name="SMA"/>

#### Single Marker Analysis 

If SINGLEMARKER=1, then the analysis conducted by GECS will perform only the single-marker test on all variants included in the analysis.

###### *NOTICE* 

In case of SINGLEMARKER=1, MAFT will be set to 0.5 and the NCT and ALLBINS parameters will be ignored.
Input and output string parameter are to specify only the name of input and output file without any extensions.
###### Resulted files:
1. A text file \[\*_singlemarker.txt\], which includes the results of association tests for all included markers.
2. Two files \[\*_singlemarker.pvals\] and \[\*_singlemarker.srt.pvals\] for the correction of multiple testing based on the null model.
###### *NOTICE* 

The final corrected alpha will be reported with other information about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.

<a name="VB"/>

#### Variable Binning
In the default case (SINGLEMARKER=0) the parameters NCT and MAFT do the same job, which is determinig the rareness threshold for the analysis.
###### *NOTICE* 
1) In case of specifying the both parameters NCT and MAFT, the NCT parameter will overwrite the parameter MAFT by default, otherwise, specifying only one parameter will be enough.
2) ALLBINS paremeter is set to 0 for the default case, where only distinct bins will be considered. This feature is only usefull in case of scanning small regions for plotting purposes.
Only 
###### Resulted files:
1. A text file \[\*_gecs_\<nct\>.txt\], which includes the results of association tests for all genomic subsequences.
2. Two files \[\*_gecs_\<nct\>.pvals\] and \[\*_gecs_\<nct\>.srt.pvals\] for the correction of multiple testing based on the null model.

<a name="CMT"/>

#### Correction for multiple testing

The correction for multiple testing is achieved by performing permutation.
the final corrected alpha will be reported with other informations about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.




### Attributions

#### Authors
George Kanoungi, Dmitriy Drichel

#### License
This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details.

#### Acknowledgment
This project was supported by the German Research Foundation grant BE 38/28/9-1. The funding organization did not have any influence on the design, conduct or conclusions of the study.

#### Citation
Please site: George Kanoungi, Michael Nothnagel and Dmitriy Drichel. An exhaustive genomic scan (2019)


