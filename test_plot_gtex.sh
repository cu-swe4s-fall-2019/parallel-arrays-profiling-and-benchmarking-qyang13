#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

echo Functional test for plot_gtex
rm *.png

python  plot_gtex.py \
--gene_reads \
GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz \
--sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt \
--gene ACTA2 \
--group_type SMTS \
--out_file ACTA2.png

run test_file_exist test -f ACTA2.png
assert_exit_code 0

rm *.png
