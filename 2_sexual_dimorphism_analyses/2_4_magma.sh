#!/bin/bash

FEMALE_GWAS="Female_Obesity_GWAS_Short.txt"
MALE_GWAS="Male_Obesity_Common_Factor_GWAS_Short.txt"
REF_DATA="g1000_eur"

module load magma/1.08

magma --annotate window=50,50 --snp-loc $FEMALE_GWAS --gene-loc NCBI37.3.gene.loc --out female_annot
magma --annotate window=50,50 --snp-loc $MALE_GWAS   --gene-loc NCBI37.3.gene.loc --out male_annot

magma --bfile $REF_DATA --pval $FEMALE_GWAS N=305573 --gene-annot female_annot.genes.annot --out female_genes
magma --bfile $REF_DATA --pval $MALE_GWAS   N=268683 --gene-annot male_annot.genes.annot   --out male_genes
