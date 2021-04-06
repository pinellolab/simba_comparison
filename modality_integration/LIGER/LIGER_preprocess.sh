#!/bin/bash

#dataset=pbmc
#gunzip -c ../data/$dataset/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz | sort  -k1,1 -k2,2n -k3,3n - > ../data/$dataset/atac_fragments.sort.bed
#sort -k 1,1 -k2,2n -k3,3n ../data/$dataset/annots/hg19_genes.bed > ../data/$dataset/annots/hg19_genes.sort.bed
#sort -k 1,1 -k2,2n -k3,3n ../data/$dataset/annots/hg19_promoters.bed > ../data/$dataset/annots/hg19_promoters.sort.bed
#
#bedmap --ec --delim "\t" --echo --echo-map-id ../data/$dataset/annots/hg19_promoters.sort.bed ../data/$dataset/atac_fragments.sort.bed > ../data/$dataset/atac_promoters_bc.bed
#bedmap --ec --delim "\t" --echo --echo-map-id ../data/$dataset/annots/hg19_genes.sort.bed ../data/$dataset/atac_fragments.sort.bed > ../data/$dataset/atac_genes_bc.bed
#
dataset=mouse-skin
#gunzip -c ../data/$dataset/GSM4156597_skin.late.anagen.atac.fragments.bed.gz > ../data/$dataset/GSM4156597_skin.late.anagen.atac.fragments.bed 
#sort -k1,1 -k2,2n -k3,3n --parallel=20 ../data/$dataset/GSM4156597_skin.late.anagen.atac.fragments.bed > ../data/$dataset/atac_fragments.sort.bed
sort -k 1,1 -k2,2n -k3,3n ../data/$dataset/annots/mm10.genes.bed > ../data/$dataset/annots/mm10_genes.sort.bed
sort -k 1,1 -k2,2n -k3,3n ../data/$dataset/annots/mm10.promoters.bed > ../data/$dataset/annots/mm10_promoters.sort.bed

bedmap --ec --delim "\t" --echo --echo-map-id ../data/$dataset/annots/mm10_genes.sort.bed ../data/$dataset/atac_fragments.sort.bed > ../data/$dataset/atac_promoters_bc.bed
bedmap --ec --delim "\t" --echo --echo-map-id ../data/$dataset/annots/mm10_genes.sort.bed ../data/$dataset/atac_fragments.sort.bed > ../data/$dataset/atac_genes_bc.bed


