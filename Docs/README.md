
# Genomic Exhaustive Collapsing Scan (GECS)
 ###### _An exhaustive genomic scan for association in genetic data_

### Table of contents
* [Installation](#Instalation)
* [Description](#Description)
* [Components](#Components)
  * [Single marker analysis](#SMA)
  * [Variable binning](#VB)
  * [Correction for multple testing](#CMT)
* [Getting started](#GS)
  * [Examples](#Examples)
* [Attributions](#Attributions)
  * [Authors](#Authors)
  * [License](#License)
  * [Aknowledgment](#Acknowledgment)
  * [Citation](#Citation)

### Instalation

git clone --recursive https://github.com/ddrichel/GECS

cd GECS

./configure?

make

### Description 
We developed a new approach to conduct association analysis for rare variants exhaustively in whole-genome or whole-exome data sets, by variating bins sizes and MAF thresholds. GECS is an ultra fast program that perform an exhaustive scan for association in case-control genetic data, implemented in C++.

##### Prerequisites

GECS is distributed under GPL3 license. Starting from GECS 1.1, it supports C++ (?) on linux systems.
This program uses the alglib C++ library.

##### Usage

The usage of gecs is very simple. Only execute the follwoing command.

**gecs** <~/path/to/**file.param**>

**gecs --help** will present a breef description of the important parameters to specify in the parameter file.

Keywords in the parameter file [\*.param](https://github.com/ddrichel/GECS/tree/master/Docs/DATA/example_1.param) :

|Key _\<Input\>_|Default|Description|Notes|
|---------------|-------|-----------|---------------------------------|
|**BFILE**  _\<string\>_|none| prefix of the plink binary file)|With this parmater the name of the binary plink files \*.bed, \*.bim, and \*.fam of the data set will be specified.| 

   

**SINGLEMARKER**	_\<bool\>_		  (whether single-marker analysis should be performed instead of VB (default=0))   

GECS can perform two kind of analysis. The standrad single marker analysis (SMA) and the variable binnig approarch (VB). The default option for GECS is the VB approach, where the parameter SINGLEMARKER=0. If SINGLEMARKER=1, then the single marker analysis will be performed. In  oposite to SMA, where all variants regardless of their minor allele frequencies, in VB approach we need to specify the threshold of variants included in the analysis. In other words defining the rareness of the included variants. GECS use the threshold of number of cariers (NCT) to derive the minor allele frequnecy threshold (MAFT). Therefore we need to specify at least one of the fllowing two equivalent parameters.

**NCT**		_\<int\>_		           ("rareness" threshold: max. number of carriers per variant)

**MAFT** _\<double\>_          (Minor allele fequency threshold for rare variants) 

**PERMUTATIONS**	_\<int\>_		   (Number of permutations for correction of multiple testing, default=999)

In order to control the familiy wise error rate (FWER) at the 5% level, gecs performs in default 999 permutations with respect to the case-control labels. In each permutation, we will obtaine the smallest p-value calculated in the analysis, so at the end we will have a list of 999 p values.

**CORRECTED_P** _\<bool\>_     (Whether calculated p values will be corrected by wilson score interval of CI 95%, default=0) 

Using the results of the FWER approach we can correct the results of the analysis and calculate the upper and lower limits of the confidence interval of 95% usinf the wilson score interval.

**PTHRESHOLD**	_\<double\>_		  (max. nominal p-value for bins to be written to output files, default=1)
By specifying this parameter you can restrict the results on only p-values less than a specific threshold. This feature is usefull in case of having a limited memory size.

**ALLBINS**		_\<bool\>_		      (whether locally not-distinct bins should be written to output (useful for plotting, default=0))

**OR**		_\<bool\>_		           (Whether odd ratios will be calculated, default=0)

**OUTPUT**		_\<string\>_ 		    (Prefix of output files, default=prefix of the plink binary file of the analyzed data)

### Components

GECS provides two major features for conducting association analysis for rare variants, namely for signle markers and for all possible bins (subsequences of contiguous markers) in the genetic data set. In the single marker analysis, all variants will be considered in the analysis, regardless of their frequencies. however, in the variable binning approach we need to specify a threshold of minor allele frequency. Permutations with respect to the case-control labels is applied to make correction for multiple testing. That means if PERMUTATIONS==0, then there is no correction for multiple testing will be done. Moreover, yoe have the possibility to get the corrected p values by wilson score interval for conficence interval of 95%. You have the option to calcutae the odds ratios for all bins by specifying OR 1. (OR=0 is by default)

<a name="SMA"/>

#### Single Marker Analysis 

If SINGLEMARKER==1, then the analysis conducted by GECS will perform only the single-marker test on all variants included in the analysis.

###### *NOTICE* 

In case of SINGLEMARKER==1, MAFT will be set to the maximum and the NCT and ALLBINS parameters will be ignored.
Input and output string parameter are to cpecify only the name of input and output file without any extensions.
###### Resulted files:
1. A text file \[\*_singlemarker.txt\], which includes the results of association tests for all included markers.
2. Two files \[\*_singlemarker.pvals\] and \[\*_singlemarker.srt.pvals\] for the correction of multiple testing based on the null model.
###### *NOTICE* 

The final corrected alpha will be reported with other informations about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.

<a name="VB"/>

#### Variable Binning
In the default case (SINGLEMARKER==0) the parameters NCT and MAFT do the same job, which determinig the rareness threshold for the analysis.
###### *NOTICE* 
1) In case of specifying the both parameters NCT and MAFT, the NCT parameter will overwrite the parameter MAFT by default, otherwise, specifying only one parameter will be enough.
2) ALLBINS paremeter is set to 0 for the default case, where only distinct bins will be considered. This feature is only usefull in case of scanning small regions for plotting purposes.
Only 
###### Resulted files:
1. A text file \[\*_gecs_\<nct\>.txt\], which includes the results of association tests for all genomic subsequences.
2. Two files \[\*_gecs_\<nct\>.pvals\] and \[\*_gecs_\<nct\>.srt.pvals\] for the correction of multiple testing based on the null model.

<a name="CMT"/>

#### Correction for maultiple testing

The correction for multiple testing is achieved by performing permutation 
the final corrected alpha will be reported with other informations about the analysis in the \[\*_gecs_nct_\<nct\>.log\] file.



<a name="GS"/>

### Getting started

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
This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details

#### Acknowledgment
This project was supported by the German Research Foundation grant BE 38/28/9-1. The funding organization did not have any influence on the design, conduct or conclusions of the study.

#### Citation
Please site: George Kanoungi, Michael Nothnagel and Dmitriy Drichel. An exhaustive genomic scan (2019)


