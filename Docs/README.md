
# Genomic Exhaustive Collapsing Scans (GECS)
#### -- Exhaustive scan for association in genmoic data

### Description: 
We developed a new approach to conduct association analysis for rare variants exhaustively in whole-genome or whole-exome data sets, by variating bins sizes and MAF tresholds.

* Citation
Pleadse site George Kanoungi, Michael Nothnagel and Dmitriy Drichel (2018)
- Exhaustive scan in Genome-Wide Association studies -

* Table of contents

Getting strated
Prerequisites
Instalation
Usage
Deployment
Versioning
Authors
License
Acknowledgments
 
* Getting strated

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

* Prerequisites

GECS is distributed under GPL3 license. Starting from GECS 1.1.1, it supports c++ (?) on linux systems.
This program uses the alglib c++ library.

* Instalation

clone from github ddrichel repository
execiute the make command in the cloned directory

* Usage

The usage of gecs is very simple. Only execute the follwoing command.

gecs <~/path/to/file.param>

Keywords in the parameter file [*/file.param] :

BFILE		string		  # prefix of the plink binary file
SINGLEMARKER	bool		  # whether single-marker analysis should be performed instead of VB (default=0)	  
PERMUTATIONS	int		  # number of permutations
NCT		int		  # "rareness" threshold: max. number of carriers per variant
MAFT            double            # Minor allele fequency threshold for rare variants
PTHRESHOLD	double		  # max. nominal p-value for bins to be written to output files
ALLBINS		bool		  # whether locally not-distinct bins should be written to output (useful for plotting, default=0)
OR		bool		  # Whether odd ratios will be calculated or not.
OUTPUT		string 		  # prefix of output files

Interpretation

The default feature of GECS is to conduct the exhaustive collapsing scan for association (SINGLEMARKER==0).
If SINGLEMARKER==1, then GECS will conduct only the single-marker test on all variants included in the analysis.
-NOTICE- that in case of SINGLEMARKER==1, MAFT will be set to the maximum and the NCT and ALLBINS parameters will be ignored.
Input and output string parameter are to cpecify only the name of input and output file without any extensions.

In the default case (SINGLEMARKER==0) the parameters NCT and MAFT do the same job, which determinig the rareness threshold for the analysis.
-NOTICE- that sypecifying NCT will overwrite the parameter MAFT, otherwise, specifying only MAFT will be enough.
ALLBINS paremeter is set to 0 for the default case, where only distinct bins will be considered. This feature is only usefull in case of scanning small regions for plotting purposes.

For all cases, permutations will be used to make correction for multiple testing. That means if PERMUTATIONS==0, then there is no correction for multiple testing will be done.

You have the option to calcutae the odds ratios for all bins by specifying OR==1. (OR=0 is by default)
Resulted files

Depends on the aim of analysis. if SINGLEMARKER==1, then the text file *Singlemarker.txt will be generated, which includes the results of association tests.
If SINGLEMARKER==0, then the text file *VB.txt will be generated, which includes the results of the exhaustive scan for association.
Additionally, two files will be always generated, namely *.pvals and *srt.pvals for the correction of multiple testing.
-NOTICE- the final corrected alpha will be reported with other informations about the analysis in the *.log file.

* Deployment

* Versioning

We use SemVer for versioning. For the versions available, see the tags on this repository.

* Authors

Dmitriy Drichel -initial work-

* License

This project is licensed under the MIT License - see the LICENSE.md file for details

*Acknowledments

