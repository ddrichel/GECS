#!/bin/bash
set -e
# The script is a test of the exhaustive collapsing algorithm in combination with the Variable Threshold (VT) method
# It requires gecs to be recompiled (make -B) with int DEBUG=1 in scan.cpp
# Execute the script in the main GECS directory: DATA/VB_test.sh
# For more information, see the documentation

  
# 100x FT
filelist=""
for i in {1..50}
do
    cat DATA/example_template.param | sed s:XXX:${i}:g > DATA/example_template_NCT${i}.param
    ./gecs DATA/example_template_NCT${i}.param
    sed '1d;2d;' DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPS_gecs_${i}.txt | cut -f9  | sort | uniq > DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPS_gecs_${i}_indlist_uniq.txt
    rm DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPS_gecs_nct_${i}.log
    rm DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPS_gecs_${i}.txt
    filelist=${filelist}" DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPS_gecs_"${i}"_indlist_uniq.txt"
done
cat $filelist > DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBFT_indlist.txt && rm  $filelist

cat DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBFT_indlist.txt  | sort | uniq > DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBFT_indlist_uniq.txt
rm DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBFT_indlist.txt


# 1x VB
./gecs DATA/example_6.param
sed '1d;2d;' DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBVT.txt | cut -f8 | sort | uniq > DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBVT_indlist_uniq.txt
rm DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBVT.txt


diffout=$(diff -q DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBVT_indlist_uniq.txt DATA/1k_sim1_CEU_POLYSNPS_chr22_10k_SNPSVBFT_indlist_uniq.txt)


if [ "$diffout" ]
then
    echo "Output of the two test cases differs!"
else
    echo "Output of the two test cases is identical!"
fi
