
# Genomic Exhaustive Collapsing Scan (GECS)
 ###### _An exhaustive genomic scan for association in genomic data_

George Kanoungi, Michael Nothnagel, Tim Becker, Dmitriy Drichel. The exhaustive genomic scan approach, with an application to rare-variant association analysis"
https://www.biorxiv.org/content/10.1101/571752v1

### Table of contents
* [Installation](#Installation)
* [Description](#Description)
* [Getting started](#GS)
  * [Examples](#Examples)
* [Functionality](#Functionality)
  * [Single-marker analysis](#SMA)
  * [Exhaustive scan](#ES)
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
In the [Examples](https://github.com/ddrichel/GECS/tree/master/DATA) you will try different combinations of parameter specifications for a specific goal in the analysis.
1) Example 1 : You will perform a single marker analysis and use 999 permutations for correcting for multiple testing. Note no need to specify the minor allele frequency threshold (MAFT) or the threshold of number of carriers (NCT).
2) Example 2: You will perform the variable binning approach instead of the single marker analysis, and here you are obligated to specify the parameter MAFT or the equivalent parameter NCT like in the example 3.
3) Example 4: In this example a conflict in the values of NCT and MAFT will be solved by considering NCT always it is specified.
4) Example 5: The SMA approach is applied on the chromosome 22 of a simulated data set with 1000 sample size.



##### Usage

The usage of gecs is:

**gecs** <~/path/to/**file.param**>

**gecs --help** will present a brief description of the keywords in the parameter file.

Keywords in the parameter file [\*.param](https://github.com/ddrichel/GECS/tree/master/DATA/example_1.param) :

| KEY | _\<Input\>_ | Default | Description | Notes |
|----|---|-------|------|------|
|**BFILE**|_\<string\>_|none| Prefix of the plink binary file|With this parameter, the name of the binary plink files (\*.bed, \*.bim, and \*.fam) will be specified| 
|**SINGLEMARKER**|	_\<bool\>_|0|Whether the single-marker analysis (SMA) should be performed instead of the region-based exhaustive scan|If SINGLEMARKER is set to 1, the standard single marker analysis (SMA) will be performed and all variants will be included in the analysis, regardless of their minor allele frequencies. For SINGLEMARKER 0, the exhaustive scan approach will be preformed, and the threshold of the minor allele frequencies (or number of carriers per variant) needs to be specified.|
|**NCT**|_\<int\>_|none|"Rareness" threshold: max. number of carriers per variant| Either NCT or MAFT has to be specified for the exhaustive scan approach|
|**MAFT**|_\<double\>_|none|Minor allele frequency threshold for rare variants|Ignored if SINGLEMARKER is set to 1, otherwise either NCT or MAFT has to be specified| 
|**PERMUTATIONS**|_\<int\>_|999|Number of permutations for correction of multiple testing|In order to control the familiy wise error rate (FWER) at the 5% level, gecs performs a full analysis of each permutation replicate (random reassigment of case-control labels). The smallest p-value of each replicate is written to an output file. The full list of smallest p-values is used for correction for multiple testing|
|**CORRECTED_P**|_\<bool\>_|0|Whether the permutation-adjusted p-values will be written to output, together with the 95% CI (Wilson score interval) |Useful if the data set can be analyzed in a single run. Adjustment is performed for a single "rareness" threshold |
|**PTHRESHOLD**|_\<double\>_|1|Max. nominal p-value for bins to be written to output|By specifying this parameter you can restrict the results on only p-values less than a specific threshold. This feature is usefull for controlling the size of output files|
|**ALLBINS**|_\<bool\>_|0|Whether locally not-distinct bins should be written to output (useful for plotting and testing)|Use with caution and small input files, otherwise the output can be become extremely large. Therefore, no correction for multiple testing is possible for this process (PERMUTATIONS will be set to 0)|
|**OR**|_\<bool\>_|0|Whether the odd ratios should be calculated||
|**OUTPUT**|_\<string\>_|Prefix of the input file|Prefix of output files||

### Functionality

GECS provides two major modes for conducting association analysis, namely for single markers (SMA) and for all possible bins (subsequences of contiguous markers below a specified "rareness" threshold) in the genetic data set. In the single marker analysis, all variants will be considered in the analysis, regardless of their frequencies. Permutations with respect to the case-control labels is applied to make correction for multiple testing. If PERMUTATIONS is set to 0, no correction for multiple testing will be performed and only nominal p-values will be computed.

<a name="SMA"/>

#### Single Marker Analysis 

If SINGLEMARKER is set to 1, only the single-marker analysis (SMA) will be performed. The SMA can be regarded as a special case of the collapsing test COLL (Li and Leal, 2008) without a maximum "rareness" threshold and with bin size of 1. If SINGLEMARKER is set to 1, MAFT will be effectively set to 0.5, and the NCT and ALLBINS parameters will be ignored. For SMA, the COLL test is identical to Pearson's chi-squared test under the dominant model.


###### Resulting files:
1. A text file \[\*_singlemarker.txt\], which includes the results of association tests for all included markers.
2. Two files \[\*_singlemarker.pvals\] and \[\*_singlemarker.srt.pvals\] for the correction of multiple testing based on the null model.
###### *NOTICE* 

The final corrected alpha will be reported with other information about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.

<a name="ES"/>

#### Exhaustive scan
In the default case (SINGLEMARKER=0) the parameters NCT and MAFT do the same job, which is determinig the rareness threshold for the analysis.
###### *NOTICE* 
1) In case of specifying the both parameters NCT and MAFT, the NCT parameter will overwrite the parameter MAFT by default, otherwise, specifying only one parameter will be enough.
2) ALLBINS paremeter is set to 0 for the default case, where only distinct bins will be considered. This feature is only usefull in case of scanning small regions for plotting purposes.
Only 
###### Resulting files:
1. A text file \[\*_gecs_\<nct\>.txt\], which includes the results of association tests for all genomic subsequences.
2. Two files \[\*_gecs_\<nct\>.pvals\] and \[\*_gecs_\<nct\>.srt.pvals\] for the correction of multiple testing based on the null model.

<a name="CMT"/>

#### Correction for multiple testing

The correction for multiple testing is achieved by performing permutation.
the final corrected alpha will be reported in the \[\*_gecs_nct_\<nct\>.log\] file, along with additional information.




### Attributions

#### Authors
George Kanoungi, Dmitriy Drichel

#### License
This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details.

#### Acknowledgment
This project was supported by the German Research Foundation grant BE 38/28/9-1. The funding organization did not have any influence on the design, conduct or conclusions of the study.

#### Citation
Please site: George Kanoungi, Michael Nothnagel, Tim Becker, Dmitriy Drichel. The exhaustive genomic scan approach, with an application to rare-variant association analysis" (2019)
