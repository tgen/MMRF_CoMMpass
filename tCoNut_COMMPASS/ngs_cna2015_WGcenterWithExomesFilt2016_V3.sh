#!/bin/bash

## tCoNuT Pipeline
##
## MIT License
##
## Copyright (c) 2017 Translational Genomics Research Institute
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
##
## Major Contributors: Jessica Aldrich, Austin Christofferson 
## Minor Contributors: 
##
## Description -
##
##
## Usage = 
##
## Options:

# Usage and Command line options

#{{{

usage="
Program:	tCoNuT_MMRF Pipeline
Version:	v1.0.0
Summary:	TGen Copy Number Tool version specific for MMRF CoMMpass study

Usage:		$(basename "$0") [Options] -b <> -n <> -N <> -t <> -T <> -h <> -m </packages/MCR/9.0/v90>

Options:
	-N	<string>	Exome normal sample name	Required
	-T	<string>	Exome tumor sample name		Required
	-b	<file>		Padded targets file			Required
	-n	<file>		Genome Normal Dat file		Required
	-t	<file>		Genome Tumor Dat file		Required
	-h	<file>		Exome Normal HS Metrics		Required
	-H	<file>		Exome Tumor HS Metrics		Required
	-v	<file>		Exome Haplotype Caller VCF	Required
	-m	<path>		MCR path					Required
	-d	<file>		dbSNP vcf					Required
	-s	<file>		SnpSift jar file			Required
	-c	<file>		CCDS gtf					Required
	-f	<file>		Copy Number targets file	Required
	-o 	<string>	Out file prefix				Required

"

if [ "$*" == "--help" ] || [ "$*" == "-h" ] || [ "$*" == "" ]
then
	echo "${usage}"
	exit 0
fi

while getopts ":b:n:N:t:T:h:m:d:s:c:f:v:H:o:" opt
do
	case $opt in
		N)
			NORMALSAMPLE=$OPTARG ;;
		T)
			TUMORSAMPLE=$OPTARG ;;
		b)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				BEDFILE=$OPTARG
			fi ;;
		n)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				NORMALDAT=$OPTARG
			fi ;;
		h)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				NORMALHS=$OPTARG
			fi ;;
		t)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				TUMORDAT=$OPTARG
			fi ;;
		H)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				TUMORHS=$OPTARG
			fi ;;
		v)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				VCF=$OPTARG
			fi ;;
		m)
			if [ ! -d $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				MCRPATH=$OPTARG
			fi ;;
		d)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				DBSNP_VCF=$OPTARG
			fi ;;
		s)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				SNPSIFT=$OPTARG
			fi ;;
		c)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				CCDSLIST=$OPTARG
			fi ;;
		f)
			if [ ! -f $OPTARG ]
			then
				echo
				echo $OPTARG is not a valid file path.
				echo Please check the path to your input list.
				exit 1
			else
				TARGETSFILE=$OPTARG
			fi ;;
		o)
			OFILE=$OPTARG ;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1 ;;
		\?)
			echo "Invalid option: -$OPTARG" >&2 ;;

	esac
done

if [ -z $NORMALSAMPLE ]
then
	echo "The -N option is required."
	exit 1
elif [ -z $TUMORSAMPLE ]
then
	echo "The -T option is required."
	exit 1
elif [ -z $BEDFILE ]
then
	echo "The -b option is required."
	exit 1
elif [ -z $NORMALDAT ]
then
	echo "The -n option is required."
	exit 1
elif [ -z $NORMALHS ]
then
	echo "The -h option is required."
	exit 1
elif [ -z $TUMORDAT ]
then
	echo "The -t option is required."
	exit 1
elif [ -z $TUMORHS ]
then
	echo "The -H option is required."
	exit 1
elif [ -z $VCF ]
then
	echo "The -v option is required."
	exit 1
elif [ -z $MCRPATH ]
then
	echo "The -m option is required."
	exit 1
elif [ -z $DBSNP_VCF ]
then
	echo "The -d option is required."
	exit 1
elif [ -z $SNPSIFT ]
then
	echo "The -s option is required."
	exit 1
elif [ -z $CCDSLIST ]
then
	echo "The -c option is required."
	exit 1
elif [ -z $TARGETSFILE ]
then
	echo "The -f option is required."
	exit 1
elif [ -z $OFILE ]
then
	echo "The -o option is required."
	exit 1
fi

# Confirm all required tools are available
type bedtools >/dev/null 2>&1 || { echo >&2 "Require \"\bedtools\"  but it's not in the PATH.  Aborting."; exit -1; }
type perl >/dev/null 2>&1 || { echo >&2 "Require \"\perl\"  but it's not in the PATH.  Aborting."; exit -1; }
type java >/dev/null 2>&1 || { echo >&2 "Require \"\java\"  but it's not in the PATH.  Aborting."; exit -1; }
type Rscript >/dev/null 2>&1 || { echo >&2 "Require \"\Rscript\"  but it's not in the PATH.  Aborting."; exit -1; }


#}}}

#module load VCFtools/0.1.10		Not sure if this is even used.

#module load BEDTools/2.26.0
#module load perl/5.14.2
#module load R/3.0.0

#export PERL5LIB=$PERL5LIB:/home/jaldrich/perl5/lib/perl5


SCRIPTDIR=$(dirname $(readlink -f "$0"))
ZTABLE=${SCRIPTDIR}/ztable.txt
GENDER_GRAPH_R=${SCRIPTDIR}/Gender_Graph.R

FAILED=0

# Possible fix for matlab runtime
TD=`pwd`
mkdir TempDir
TMPDIR=${TD}/TempDir
TMP=${TD}/TempDir
TEMPDIR=${TD}/TempDir

export MCR_CACHE_ROOT=${TD}/TempDir

# Filter vcf to Targets

cat ${VCF} | java -jar ${SNPSIFT} intervals ${BEDFILE} > ${VCF}.temp.vcf

if [ $? -ne 0 ]
then
		FAILED=1
fi

mv ${VCF} ${VCF}.original
mv ${VCF}.temp.vcf ${VCF}

# Filter out pos35, waterfall type 2 from the vcf

cat ${VCF} | java -jar ${SNPSIFT} intervals ${SCRIPTDIR}/CNA_type2_and_pos35_filter_X_fixed_intervals_to_keep_X_Y.bed > ${VCF}.temp.vcf

if [ $? -ne 0 ]
then
		FAILED=1
fi

mv ${VCF} ${VCF}.filtered.original
mv ${VCF}.temp.vcf ${VCF}

# Gender correction

mkdir genderTest

#Haplotyper Caller outputs the genotype columns in alphabetical-numeric order.	Need to determine if the order is Tumor-Normal (expected) or Normal-Tumor (it happens)
echo "The following line contains the Haplotyper Caller VCF header:"
grep -v "##" ${VCF} | grep "#"
echo Extracting the first genotype column to determine if it is tumor or normal
FIRST_GENOTYPE_COLUMN=`grep -v "##" ${VCF} | grep "#" | cut -f10`
echo The first genotype column header is: ${FIRST_GENOTYPE_COLUMN}
FRACTION_LETTER=`echo ${FIRST_GENOTYPE_COLUMN} | cut -d_ -f6 | cut -c1`
echo The fraction letter in the first genotype column is: ${FRACTION_LETTER}

if [ "${FRACTION_LETTER}" == "T" ]
then 
		echo Found expected genotype order - Proceeding with filtering
		cat ${VCF} | java -jar ${SNPSIFT} filter \
		"( GEN[1].DP > 50 ) & isVariant( GEN[1] ) \
		& ((REF = 'A') | (REF = 'C') | (REF = 'G') | (REF = 'T')) \
		& ((ALT = 'A') | (ALT = 'C') | (ALT = 'G') | (ALT = 'T')) \
		& (exists ID) & ( ID =~ 'rs' ) \
		& ( CHROM = 'X' )" | \
		java -jar ${SNPSIFT} annotate -info GMAF ${DBSNP_VCF} - | \
		java -jar ${SNPSIFT} filter "( GMAF >= 0.05 )" | \
		java -jar ${SNPSIFT} extractFields -e "." - \
		CHROM POS ID REF ALT FILTER GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP | \
		awk 'NR>1' | \
		awk 'BEGIN{OFS = "\t" ; print "Chr", "Pos", "Normal_Ratio", "Normal_DP"}{OFS = "\t" ; print $1, $2, $8/($7+$8), $9}' > GenderDataTable.txt
		
elif [ "${FRACTION_LETTER}" == "C" ]
then 
		echo Found the wrong genotype order - Reordering genotype columns and then Proceeding with filtering
		awk -F'\t' '{OFS="\t" ; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$10}' ${VCF} | java -jar ${SNPSIFT} filter \
		"( GEN[1].DP > 50 ) & isVariant( GEN[1] ) \
		& ((REF = 'A') | (REF = 'C') | (REF = 'G') | (REF = 'T')) \
		& ((ALT = 'A') | (ALT = 'C') | (ALT = 'G') | (ALT = 'T')) \
		& (exists ID) & ( ID =~ 'rs' ) \
		& ( CHROM = 'X' )" | \
		java -jar ${SNPSIFT} annotate -info GMAF ${DBSNP_VCF} - | \
		java -jar ${SNPSIFT} filter "( GMAF >= 0.05 )" | \
		java -jar ${SNPSIFT} extractFields -e "." - \
		CHROM POS ID REF ALT FILTER GEN[1].AD[0] GEN[1].AD[1] GEN[1].DP | \
		awk 'NR>1' | \
		awk 'BEGIN{OFS = "\t" ; print "Chr", "Pos", "Normal_Ratio", "Normal_DP"}{OFS = "\t" ; print $1, $2, $8/($7+$8), $9}' > GenderDataTable.txt
else	
		#This should not happen
		echo WARNING!!! ERROR - PLEASE LOOK WE DID NOT FIND A T or C, what the hell!
		exit 1
fi		

#Calculate Mean and SD
MEAN=`awk 'NR>1' GenderDataTable.txt | awk '{sum+=$3} END {print sum/NR}'`
SD=`awk 'NR>1' GenderDataTable.txt | awk '{sum+=$3; array[NR]=$3} END {for(x=1;x<=NR;x++){sumsq+=((array[x]-(sum/NR))**2);}print sqrt(sumsq/NR)}'`
ALLELES=`awk 'NR>1' GenderDataTable.txt | wc -l | cut -f1`

if [[ $MEAN > 0.8 ]]
then
	MALE=yes
	awk -F'\t' '{ if ($1 == 23 || $1 == 24) { OFS = "\t" ; print $1,$2,$3*2 } else { print $0 }}' ${NORMALDAT} > ${NORMALDAT}.male
	mv ${NORMALDAT} ${NORMALDAT}.original
	mv ${NORMALDAT}.male ${NORMALDAT}
fi

NORMX=`cat ${NORMALHS} | grep -v ^# | grep -v ^BAIT | cut -f 22`
echo $NORMX
TUMORX=`cat ${TUMORHS} | grep -v ^# | grep -v ^BAIT | cut -f 22`
echo $TUMORX

## parse filtered to targets HC VCF for BAF

${SCRIPTDIR}/parseHC_VCF_BAF.pl ${VCF}.filtered.original ${NORMALSAMPLE} ${TUMORSAMPLE}

if [ $? -ne 0 ]
then
		FAILED=1
fi

# parse filtered to pos35 and type2 HC VCF for Hets

${SCRIPTDIR}/parseHC_VCF_Hets.pl ${VCF} ${NORMALSAMPLE} ${TUMORSAMPLE}

if [ $? -ne 0 ]
then
		FAILED=1
fi

##Copy Number Analysis (MATLAB)
HETFILE=merged.vcf.txt
assayID="Genome"
#HETFILE=0
smWin=10				 #	 <<<< THIS CAN BE ADJUSTED - smoothing window size >>>>
fcThresh=0.2		   #   <<<< THIS CAN BE ADJUSTED - fold-change threshold - plot >>>>
res=10					 #	 <<<< THIS CAN BE ADJUSTED - min resolution >>>>
#maxGap=1200
maxGap=100000000		#	<<<< THIS CAN BE ADJUSTED - max distance between consecutive postions >>>>

hetDepthN=${NORMX}		#	<<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDepthT=${TUMORX}		#	<<<< THIS CAN BE ADJUSTED - min depth of diploid het position >>>>
hetDev=0.025			 #	 <<<< THIS CAN BE ADJUSTED - allowable deviation from ref allele frequency of 0.5

readDepth=100
#readDepth=$( echo "$hetDepthN * 3" | bc )			  #   <<<< THIS CAN BE ADJUSTED - min number of counts >>>
echo $readDepth
echo "time ${SCRIPTDIR}/pegasusCNA_MMRF/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE}_filt ${HETFILE} ${smWin} ${fcThresh} ${assayID} ${res} ${readDepth} ${maxGap} ${hetDepthN} ${hetDepthT} ${hetDev} ${TARGETSFILE}"

time ${SCRIPTDIR}/pegasusCNA_MMRF/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE}_filt ${HETFILE} ${smWin} ${fcThresh} ${assayID} ${res} ${readDepth} ${maxGap} ${hetDepthN} ${hetDepthT} ${hetDev} ${TARGETSFILE}

if [ $? -ne 0 ]
then
		FAILED=1
fi

Rscript --vanilla ${SCRIPTDIR}/validateCNAHets_v2.R ${OFILE}_filt.hets.tsv 

if [ $? -ne 0 ]
then
		FAILED=1
fi

if [ -f CNA_unimodal_pass ]
then
	echo
	echo ${OFILE}_filt.hets.tsv fits a unimodal distribution.
	echo Continuing CNA script.

	awk -F'\t' 'NR%5==1' ${OFILE}_filt.cna.tsv > temp.cna.tsv
	mv ${OFILE}_filt.cna.tsv ${OFILE}_filt.cna.tsv.original
	mv temp.cna.tsv ${OFILE}_filt.cna.tsv

	## CBS segmentation
	Rscript --vanilla ${SCRIPTDIR}/runDNAcopyV4_filt2016.R ${OFILE}_filt.cna.tsv ${OFILE}_filt.seg

	if [ $? -ne 0 ]
	then
			FAILED=1
	fi
else
	# keep hets.tsv original
	mv ${OFILE}_filt.hets.tsv ${OFILE}_filt.hets.tsv.original

	awk -F'\t' '{ OFS = "\t" ; print $2,$3-10000,$3+10000 }' CNA_positions_to_exclude.txt > CNA_intervals_to_exclude.bed

	bedtools merge -i CNA_intervals_to_exclude.bed > CNA_intervals_to_exclude_merged.bed

	awk -F'\t' 'NR > 1 { OFS = "\t" ; print $1,$2,$2+1,$3,$4,$5,$6,$7,$8 }' merged.vcf.txt > merged.vcf.bed

	mv merged.vcf.txt merged.vcf.txt.original

	bedtools intersect -v -a merged.vcf.bed -b CNA_intervals_to_exclude_merged.bed > temp
	
	awk -F'\t' 'BEGIN { OFS = "\t" ; print "Chromosome","Position","NormalReadDepth","NormalRefAllele","NormalAltAllele","TumorReadDepth","TumorRefAllele","TumorAltAllele" } { OFS = "\t" ; print $1,$2,$4,$5,$6,$7,$8,$9 }' temp > merged.vcf.txt

	rm temp

	time ${SCRIPTDIR}/pegasusCNA_MMRF/run_ngsCNA.sh ${MCRPATH} ${NORMALDAT} ${TUMORDAT} ${OFILE}_filt ${HETFILE} ${smWin} ${fcThresh} ${assayID} ${res} ${readDepth} ${maxGap} ${hetDepthN} ${hetDepthT} ${hetDev} ${TARGETSFILE}

	if [ $? -ne 0 ]
	then
		FAILED=1
	fi

	awk -F'\t' 'NR%5==1' ${OFILE}_filt.cna.tsv > temp.cna.tsv
	mv ${OFILE}_filt.cna.tsv ${OFILE}_filt.cna.tsv.filt.original
	mv temp.cna.tsv ${OFILE}_filt.cna.tsv
	
	## CBS segmentation
	Rscript --vanilla ${SCRIPTDIR}/runDNAcopyV4_filt2016.R ${OFILE}_filt.cna.tsv ${OFILE}_filt.seg

	if [ $? -ne 0 ]
	then
		FAILED=1
	fi
fi	

##plotting
Rscript --vanilla ${SCRIPTDIR}/plotCGH.R ${OFILE}_filt.cna.tsv ${OFILE}_filt.amp.tsv ${OFILE}_filt.del.tsv ${OFILE}_filt

if [ $? -ne 0 ]
then
	FAILED=1
fi

if [ -f ${OFILE}_filt.hets.tsv ]
then
	Rscript --vanilla ${SCRIPTDIR}/plotCGHwithHets.R ${OFILE}_filt.cna.tsv ${OFILE}_filt.amp.tsv ${OFILE}_filt.del.tsv ${OFILE}_filt.hets.tsv ${OFILE}_filt_withhets

	if [ $? -ne 0 ]
	then
			FAILED=1
	fi
fi

Rscript --vanilla ${SCRIPTDIR}/plotBAF.R baf.txt ${OFILE}_filt.baf

if [ $? -ne 0 ]
then
		FAILED=1
fi

Rscript --vanilla ${SCRIPTDIR}/runDNAcopyBAF.R baf.txt ${OFILE}_filt.baf

if [ $? -ne 0 ]
then
		FAILED=1
fi

##Annotate and convert SEG file to gVCF 
DUPTHRESH=0.58	   #   <<<< THIS CAN BE ADJUSTED - Amplification Threshold - log2 fold-change >>>>
DELTHRESH=-0.99    #   <<<< THIS CAN BE ADJUSTED - Deletion Threshold - log2 fold-change >>>>

${SCRIPTDIR}/annotSeg.pl ${CCDSLIST} ${OFILE}_filt.cna.seg ${DUPTHRESH} ${DELTHRESH}

if [ $? -ne 0 ]
then
		FAILED=1
fi

${SCRIPTDIR}/validateCNAVariantsVCF.pl ${OFILE}_filt.cna.seg.vcf baf.txt ${ZTABLE}

if [ $? -ne 0 ]
then
		FAILED=1
fi

# Plot Linear CNA and BAF

# Removed error capture on Matlab image creation. The image may not be created depending on sample quality and depth.
${SCRIPTDIR}/plotLinearCNA/run_plotLinearCNAandBAF.sh ${MCRPATH} ${OFILE}_filt.cna.tsv baf.txt ${OFILE}_filt.cnaBAF.png

if [ $? -ne 0 ]
then
		FAILED=1
fi

# Plot Linear CNA and Abs BAF

${SCRIPTDIR}/plotLinearCNAabsBAF/run_plotLinearCNAandAbsBAF.sh ${MCRPATH} ${OFILE}_filt.cna.tsv baf.txt ${OFILE}_filt.cnaAbsBAF.png

# Removed error capture on Matlab image creation. The image may not be created depending on sample quality and depth.
if [ $? -ne 0 ]
then
		FAILED=1
fi

exit $FAILED

