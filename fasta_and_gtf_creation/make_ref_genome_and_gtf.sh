#!/bin/bash

##
## make_ref_genome_and_gtf.sh
##
## Created by Austin Christofferson on 06/12/2017.
## Copyright 2017 Translational Genomics Research Institute. All rights reserved.
##
## Description - Creates the custom reference genome, gtf, and related tables used in
##               the tgen05 recipe of the jetstream pipeline at TGen.
##
## Usage = Make_gtf_and_related_tables.sh [--help -h]
##
####################################################################
####################################################################

# Set help options

if [ "$*" == "--help" ] || [ "$*" == "-h" ]
then
        grep "^##" $0
	echo
        exit 0
fi

BASEDIR=$(dirname "$0")

#### 1000G reference assembly with ERCC spike-ins, Human ribosomal DNA complete repeting unit, and Cancer associated onco viruses.  

# 1) Download phase2 reference assembly and associated files
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/*

# 2) Unzip the hs37d5.fa.gz reference assembly and edit EBV wild-type genome to be consistent with MT entry
	gunzip hs37d5.fa.gz
	sed 's/NC_007605/NC_007605 gi|82503188|ref|NC_007605.1| Human herpesvirus 4 complete wild type genome/g' hs37d5.fa > hs37d5.anno.fa

# 3) Download the ERCC spike-ins fasta file
	wget https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip
	unzip ERCC92.zip

# 4) Download the Human ribosomal DNA complete repeting unit and edit to be consistent with MT entry 
	wget -O U13369.1.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=U13369.1&amp;rettype=fasta"
	sed 's/U13369.1 Human ribosomal DNA complete repeating unit/U13369 gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit/' U13369.1.fa > U13369.1.anno.fa
	sed -i '/^[[:blank:]]*$/d' U13369.1.anno.fa

# 5) Download Cancer associated onco viruses and edit files to be consistant with MT geneome entry
		
	ACCESSION_LIST="`echo \
				FR872717.1 \
				AF092932.1 \
				NC_001526.2 \
				NC_001357.1 \
				NC_003977.1 \
				NC_004102.1 \
				NC_009823.1 \
				NC_009824.1 \
				NC_009825.1 \
				NC_009826.1 \
				NC_009334.1 \
				NC_006273.2 \
				NC_011800.1 \
				NC_001806.1`"

	# Download
	for line in `echo $ACCESSION_LIST ` 
	do 
		wget -O ${line}.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore;id=$line;rettype=fasta"
	done	

	# Convert to lower case to keep consistency with pipeline reference geneome
	gawk '$0 !~ />/ {$0=tolower($0) } ; { print $0 }' NC_001357.1.fa > NC_001357.1.lower.fa
	gawk '$0 !~ />/ {$0=tolower($0) } ; { print $0 }' NC_001526.2.fa > NC_001526.2.lower.fa

	# Edit 
	sed 's/FR872717.1 Human papillomavirus type 11 complete genome/FR872717 gi|335334258|emb|FR872717.1| Human papillomavirus type 11 complete genome/' FR872717.1.fa >> Cancer_Related_Viruses.fa
	sed 's/AF092932.1 Human papillomavirus type 6, complete genome/AF092932 gi|6002612|gb|AF092932.1| Human papillomavirus type 6, complete genome/' AF092932.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_001526.2 Human papillomavirus type 16, complete genome/NC_001526 gi|310698439|ref|NC_001526.2| Human papillomavirus type 16, complete genome/' NC_001526.2.lower.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_001357.1 Human papillomavirus - 18, complete genome/NC_001357 gi|9626069|ref|NC_001357.1| Human papillomavirus - 18, complete genome/' NC_001357.1.lower.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_003977.1 Hepatitis B virus, complete genome/NC_003977 gi|21326584|ref|NC_003977.1| Hepatitis B virus, complete genome/' NC_003977.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_004102.1 Hepatitis C virus genotype 1, complete genome/NC_004102 gi|22129792|ref|NC_004102.1| Hepatitis C virus genotype 1, complete genome/' NC_004102.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_009823.1 Hepatitis C virus genotype 2, complete genome/NC_009823 gi|157781212|ref|NC_009823.1| Hepatitis C virus genotype 2, complete genome/' NC_009823.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_009824.1 Hepatitis C virus genotype 3, genome/NC_009824 gi|157781216|ref|NC_009824.1| Hepatitis C virus genotype 3, genome/' NC_009824.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_009825.1 Hepatitis C virus genotype 4, genome/NC_009825 gi|157781208|ref|NC_009825.1| Hepatitis C virus genotype 4, genome/' NC_009825.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_009826.1 Hepatitis C virus genotype 5, genome/NC_009826 gi|157781210|ref|NC_009826.1| Hepatitis C virus genotype 5, genome/' NC_009826.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_009334.1 Human herpesvirus 4, complete genome/NC_009334 gi|139424470|ref|NC_009334.1| Human herpesvirus 4, complete genome/' NC_009334.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_006273.2 Human herpesvirus 5 strain Merlin, complete genome/NC_006273 gi|155573622|ref|NC_006273.2| Human herpesvirus 5 strain Merlin, complete genome/' NC_006273.2.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_011800.1 Human T-lymphotropic virus 4, complete genome/NC_011800 gi|219552908|ref|NC_011800.1| Human T-lymphotropic virus 4, complete genome/' NC_011800.1.fa >> Cancer_Related_Viruses.fa
	sed 's/NC_001806.1 Human herpesvirus 1, complete genome/NC_001806 gi|9629378|ref|NC_001806.1| Human herpesvirus 1, complete genome/' NC_001806.1.fa >> Cancer_Related_Viruses.fa

	# remove blank lines 
	sed -i '/^[[:blank:]]*$/d' Cancer_Related_Viruses.fa

# 6) Concat modified files together to make final custom reference sequence
	cat hs37d5.anno.fa U13369.1.anno.fa Cancer_Related_Viruses.fa ERCC92.fa > hs37d5_plusRibo_plusOncoViruses_plusERCC.fa

#### Ensemble Version 74 - Processing 

# 1) Download Ensembl v74 GTF
        wget ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz

# 2) Uncompress GTF
        gunzip Homo_sapiens.GRCh37.74.gtf.gz

# 3) Determine which chromosomes are present in this GTF and remove contigs that are not in the reference
	cut -f1 Homo_sapiens.GRCh37.74.gtf | sort | uniq > Ensembl_v74_Unique_Chromosomes.txt

	# Remove all contigs not in the referenece genome
        # hs37d5_plusRibo_plusOncoViruses_plusERCC.fa

	grep "^>" hs37d5_plusRibo_plusOncoViruses_plusERCC.fa | \
		awk '{print $1}' | \
		cut -d">" -f2 > hs37d5_plusRibo_plusOncoViruses_plusERCC.fa.chrom.contigs.txt

	gawk -F'\t' 'FNR==NR{a[$1]=$0;next} ($1 in a) { print $0 }' hs37d5_plusRibo_plusOncoViruses_plusERCC.fa.chrom.contigs.txt Homo_sapiens.GRCh37.74.gtf > Homo_sapiens.GRCh37.74.gtf.hs37d5.gtf

# 4) Add the ERCC contigs to the end of the GTF
	cat Homo_sapiens.GRCh37.74.gtf.hs37d5.gtf ERCC92.gtf > Homo_sapiens.GRCh37.74.gtf.hs37d5.ERCC92.gtf

# 5) Add in the EGFRvIII gtf right before EGFR-AS1-001
	gawk -F'\t' 'BEGIN{ STOP=0 ; OFS = "\t" } 
			FNR==NR{ a[NR]=$0; next } 
			FNR == 1 { n = asort(a) } ; 
			{ if (STOP == 0 && $0 ~ /EGFR-AS1-001/) 
				{ STOP=1 ; 
				for ( i = 1 ; i <= n ; i++ ) 
					{ print a[i] } ; 
					print $0  } 
			else { print $0 }}' ${BASEDIR}/EGFRvIII.gtf Homo_sapiens.GRCh37.74.gtf.hs37d5.ERCC92.gtf > Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf
 

