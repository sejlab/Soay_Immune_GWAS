# Soay_Immune_GWAS
 
This repository is for code related to the manuscript:

The genetic architecture of helminth-specific immune responses in a wild
population of Soay sheep (Ovis aries).

A. M. Sparks, K. Watt, R. Sinclair, J. G. Pilkington, J. M. Pemberton, 
T. N. McNeilly, D. H. Nussey, S. E. Johnston

In press at PLoS Genetics. For preprint see:  https://www.biorxiv.org/content/10.1101/628271v1

For questions related specifically to these scripts and dataset, please contact Susan
Johnston (Susan.Johnston@ed.ac.uk) or Alexandra Sparks (A.M.Sparks@leeds.ac.uk).

All analysis scripts are provided atin this repository and use the same directory structure as outlined in the file paths below. 

Scripts are run in the following order:

* 1.1_BeastAnModGRM070218.R
* 1.2_Format_Results.R
* 1.3_BEAST_Genetic_Correlations.R
* 1.4_Parse_Genetic_Correlations.R
* brute_GWAS/2_BEAST_GWAS_Prepare.R 
* brute_GWAS/2_BEAST_GWAS_Run.R
* brute_GWAS_log_Lambs/2_BEAST_GWAS_Prepare.R 
* brute_GWAS_log_Lambs/2_BEAST_GWAS_Run.R
* 2.2.2_BEAST_GWAS_Parse_BiomaRt_log_Lambs.R
* 2.2_BEAST_GWAS_Parse_BiomaRt.R
* 2.3_Reparse_and_format_gene_GO_tables.R
* 3.1_Imputation_in_Sig_Regions.R
* 3.1_Imputation_in_Sig_Regions_log_Lambs.R
* alphaimpute/3_Imputed_GWAS_Prepare.R
* alphaimpute/3_Imputed_GWAS_Run.R
* alphaimpute/3_Imputed_GWAS_Run_log_Lambs.R
* 3.2_Parse_Genotypes_Prepare_for_GWAS.R
* 3.2_Parse_Genotypes_Prepare_for_GWAS_log_Lambs.R
* 3.3_Parse_Alphaimpute_GWAS.R
* 3.3_Parse_Alphaimpute_GWAS_log_Lambs.R
* 4_Examine_GWAS_Results_LD_Genes.R 
* 5_Variance_Attributed_to_Regions.R
* 6_BiomaRt_SNP_functions.R
* 6_BLAST_MHC_with_Rambouillet_Genome.R
* 7_Sex_Specific_Effects_v1.R
* 7_SNP_Effects.R

The data files are archived on the Dryad repository and are described as follows:
[ LINK WILL BE UPDATED ]

* data/20181107_BEAST_data_formatted.RData

This is an .RData file containing the data frame `pedigree` (ID, MOTHER, FATHER
information) and data frame `BEASTX` which contains individual immune phenotype
measures.
[N.B. All .RData files can be read into R using the command: load("FILENAME")]


* data/20181107_Genabel_Data.RData

This is an .RData file containing the 50K genotype data in GenABEL format after
quality control (N.B. requires R v3.4.3 and the R library GenABEL v1.8-0).

* data/20181107_SoayGRM.Rdata

This is an .RData file containing the genomic relatedness matrix used for
estimating trait heritability.

* data/Ovis_aries.Oar_v3.1.94.genes.txt

This is a text file containing the gene information for Ovis aries, as obtained
from Ensembl.

* data/20180209_SoayPlates1-83.bed/.bim/.fam

These are the PLINK files containing the 50K genotype data after quality 
control.

* data/20180209_SoayPlates1-83.bed/.bim/.fam

These are the PLINK files containing the 50K genotype data after quality 
control.

* data/20140214_SheepHD_QC1_Polym.bed/.bim/.fam

These are the PLINK files containing the HD genotype data after quality control.





