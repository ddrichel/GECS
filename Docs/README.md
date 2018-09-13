* Project name: 
Genomic Exhaustive Collapsing Scan Statistics (GECS)

* Description: 
We developed a new approach to conduct association analysis for rare variants exhaustively in whole-genome or whole-exome data sets, by variating bins sizes and MAF tresholds.

* Table of contents

* Instalation

* Usage

./GECS file.param

Keywords in the parameter file:

BFILE	 	string		  # prefix of the plink binary file
SINGLEMARKER	bool		  # whether single-marker analysis should be performed	  
PERMUTATIONS	int		  # number of permutations
NCT		int		  # "rareness" threshold: max. number of carriers per variant
PTHRESHOLD	double		  # max. nominal p-value for bins to be written to output files
ALLBINS		bool		  # whether locally don-distinct bins should be written to output (useful for plotting, default=0)
OUTPUT		string 		  # prefix of output files


# License



# Update2
