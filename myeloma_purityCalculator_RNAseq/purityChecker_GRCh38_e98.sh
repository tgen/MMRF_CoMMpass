#######################################################################
#######################################################################
## RNAseq based Myeloma Purity Calculator
##
## Copyright (c) 2021 Translational Genomics Research Institute
##
## This software may be modified and distributed under the terms
## of the MIT license.  See the LICENSE file for details.
##
## Major Contributors: Jonathan J. Keats
## Minor Contributors:
#######################################################################
#######################################################################



##############################
##
## Require applications, these must be in your $PATH at runtime
##
##############################


# Samtools is required (http://www.htslib.org/)
# Needed for idxstats step
#SAMTOOLS=samtools

# Featurecounts is required (http://subread.sourceforge.net/)
# It is used because it is an extremely fast counting application
#FEATURECOUNTS=/home/jkeats/downloads/subread-1.4.6-source/bin/featureCounts

# Add R to your path at TGEN
#module load R/3.1.1
# R and the ggplot2 package are required
# Tested with R version 3.1.1 (2014-07-10) -- "Sock it to Me"
# ggplot2 version 0.9.3.1
# These must be availalbe in your path
# To Test enter the following commands
# R			# Starts R
# library(ggplot2)	# Loads ggplot2 library
# q()			# Exits R

usage()
{
  # `cat << EOF` This means that cat should stop reading when EOF is detected
  cat << EOF

  Usage: $0 [ -bph options]
  -h      Display Help
  -b      Inpute RNAseq BAM file from STAR alignment
  -p      Prepare Resource Files (Yes/No)
  -g      Input GTF (only needed to prepare resource files)
  -l      List of Genes not expressed in B-cells (only needed to prepare resource files)
  -i      Immunoglobulin Loci in GTF format (only needed to prepare resource files)
  -d      Cleanup Temp files (Yes/No) - Set to "no" to debug

EOF
# EOF is found above and hence cat command stops reading. This is equivalent to echo but much neater when printing out.
  exit 2
}

while getopts 'b:p:g:l:i:d:?h' flag
do
    case ${flag} in
        b) bam=${OPTARG};;
        p) BUILD_FILES=${OPTARG};;
        g) GTF=${OPTARG};;
        l) NON_BCELL_CONTAMINATION_GENELIST=${OPTARG};;
        i) IG_LOCI_REGIONS=${OPTARG};;
        d) REMOVE_TEMP_FILES=${OPTARG};;
        h|?) usage;;
    esac
done

echo The bam name is: $bam
echo The prepare status is: $BUILD_FILES
echo The GTF is: $GTF
echo The non B-cell Gene list is: $NON_BCELL_CONTAMINATION_GENELIST
echo The Immunoglobulin loci file is: $IG_LOCI_REGIONS
echo The File Cleanup setting is: $REMOVE_TEMP_FILES

##############################
##
## Required Variables
##
##############################


# Table for verbose tabulated results
RESULTS_TABLE=/MMRF/purityChecker.txt

# Set the build directory, this makes it more dynamic
# Capture the current script path to dynamically link to required graphing script
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
BUILD_DIR=`dirname ${SCRIPT_PATH}`
echo BUILD DIRECTORY = ${BUILD_DIR}
#BUILD_DIR=/home/jkeats/toolkit_jkeats/myelomaPurityChecker/

# Set the path to the immunoglobulin loci GTF for the respective genome version
#IG_LOCI_REGIONS=${BUILD_DIR}Immunoglobulin_GRCh37_Loci.gtf
#IG_LOCI_REGIONS=${BUILD_DIR}Immunoglobulin_GRCh38_Loci.gtf

# The geometric mean is calculated from the RPKM of the non-B cell contamination genes
CALCULATE_GEOMEAN_R=${BUILD_DIR}calculate_geomean.R

# Two graphs will be automatically generated, this code is provided but the path needs to be updated
IGH_GRAPH_R=${BUILD_DIR}igh_graph.R

# This is only used in the first run to build all the needed files
#GTF=/home/tgenref/pecan/ensembl_v74/Homo_sapiens.GRCh37.74.gtf.hs37d5.EGFRvIII.gtf

#NON_BCELL_CONTAMINATION_GENELIST=${BUILD_DIR}Non_Bcell_Contamination_GeneList_e98.txt

# These are created in the first run using the code provided below
IGH_CONSTANT_GENES=${BUILD_DIR}IgH_Constant_Genes.txt
IGH_VARIABLE_GENES=${BUILD_DIR}IgH_Variable_Genes.txt
IGK_CONSTANT_GENES=${BUILD_DIR}IgK_Constant_Genes.txt
IGK_VARIABLE_GENES=${BUILD_DIR}IgK_Variable_Genes.txt
IGL_CONSTANT_GENES=${BUILD_DIR}IgL_Constant_Genes.txt
IGL_VARIABLE_GENES=${BUILD_DIR}IgL_Variable_Genes.txt
IG_GTF=${BUILD_DIR}Immunoglobulin_RegionsToCount.gtf

## CODE TO GENERATE NEEDED FILES FROM INPUT GTF FILE ##

#BUILD_FILES=No
#BUILD_FILES=Yes

if [ ${BUILD_FILES} == "Yes" ]
then
	echo
	echo
	echo "------------------------------------------------"
	echo Building required files
	CURRENT_DIR=`pwd`
	cd ${BUILD_DIR}
	#grep "IG_" ${GTF} | grep -w "exon" | awk '{gsub(/"/,"\t",$0); print;}' | cut -f1-8,10,12,16,22 | grep -v "pseudogene" > Immunoglobulin_Positions.txt
	grep "IG_" ${GTF} | grep -w "exon" | awk '{gsub(/"/,"\t",$0); print;}' | cut -f1-8,10,14,20,24 | grep -v "pseudogene" > Immunoglobulin_Positions.txt
	grep -w "IG_C_gene" Immunoglobulin_Positions.txt | awk '{OFS="\t" ; if($1=="chr14") print $11}' | sort | uniq > IgH_Constant_Genes.txt
	grep -w "IG_V_gene" Immunoglobulin_Positions.txt | awk '{OFS="\t" ; if($1=="chr14") print $11}' | sort | uniq > IgH_Variable_Genes.txt
	grep -w "IG_C_gene" Immunoglobulin_Positions.txt | awk '{OFS="\t" ; if($1=="chr2") print $11}' | sort | uniq > IgK_Constant_Genes.txt
	grep -w "IG_V_gene" Immunoglobulin_Positions.txt | awk '{OFS="\t" ; if($1=="chr2") print $11}' | sort | uniq > IgK_Variable_Genes.txt
	grep -w "IG_C_gene" Immunoglobulin_Positions.txt | awk '{OFS="\t" ; if($1=="chr22") print $11}' | sort | uniq > IgL_Constant_Genes.txt
	grep -w "IG_V_gene" Immunoglobulin_Positions.txt | awk '{OFS="\t" ; if($1=="chr22") print $11}' | sort | uniq > IgL_Variable_Genes.txt
	grep "IG_" ${GTF} | grep -w "exon" | grep -v "pseudogene" | grep "IG_C_gene" | awk '{OFS="\t" ; if($1=="chr14") print $0}' > Immunoglobulin_Regions.gtf
	grep "IG_" ${GTF} | grep -w "exon" | grep -v "pseudogene" | grep "IG_C_gene" | awk '{OFS="\t" ; if($1=="chr2") print $0}' >> Immunoglobulin_Regions.gtf
	grep "IG_" ${GTF} | grep -w "exon" | grep -v "pseudogene" | grep "IG_C_gene" | awk '{OFS="\t" ; if($1=="chr22") print $0}' >> Immunoglobulin_Regions.gtf
	grep "IG_" ${GTF} | grep -w "exon" | grep -v "pseudogene" | grep "IG_V_gene" | awk '{OFS="\t" ; if($1=="chr14") print $0}' >> Immunoglobulin_Regions.gtf
	grep "IG_" ${GTF} | grep -w "exon" | grep -v "pseudogene" | grep "IG_V_gene" | awk '{OFS="\t" ; if($1=="chr2") print $0}' >> Immunoglobulin_Regions.gtf
	grep "IG_" ${GTF} | grep -w "exon" | grep -v "pseudogene" | grep "IG_V_gene" | awk '{OFS="\t" ; if($1=="chr22") print $0}' >> Immunoglobulin_Regions.gtf
	# Create GTF gene query that should only find the specific gene in question
	for line in `cat ${NON_BCELL_CONTAMINATION_GENELIST}`
	do
	  echo gene_name \"${line}\" >> gene_list
	done
	grep -w "exon" ${GTF}  | grep -v "pseudogene" | grep -f ${NON_BCELL_CONTAMINATION_GENELIST} > Non_Bcell_Contamination_GeneList.gtf
	## Merge the provided Loci list with the calculated list of Ig elements, update the kappa elements from the inversed segment "D" so they will merge automatically
	cat ${IG_LOCI_REGIONS} Immunoglobulin_Regions.gtf Non_Bcell_Contamination_GeneList.gtf | sed 's/IGKV1D-12/IGKV1-12/g' | sed 's/IGKV1D-16/IGKV1-16/g' | sed 's/IGKV1D-17/IGKV1-17/g' | sed 's/IGKV1D-33/IGKV1-33/g' | sed 's/IGKV1D-37/IGKV1-37/g' | sed 's/IGKV1D-39/IGKV1-39/g' | sed 's/IGKV1D-8/IGKV1-8/g' | sed 's/IGKV2D-24/IGKV2-24/g' | sed 's/IGKV2D-28/IGKV2-28/g' | sed 's/IGKV2D-30/IGKV2-30/g' | sed 's/IGKV2D-40/IGKV2-40/g' | sed 's/IGKV3D-11/IGKV3-11/g' | sed 's/IGKV3D-15/IGKV3-15/g' | sed 's/IGKV3D-20/IGKV3-20/g' | sed 's/IGKV3D-7/IGKV3-7/g' | sed 's/IGKV6D-21/IGKV6-21/g' > Immunoglobulin_RegionsToCount.gtf
	rm Immunoglobulin_Regions.gtf
	rm Non_Bcell_Contamination_GeneList.gtf
	cd ${CURRENT_DIR}
	echo
	echo Build complete, the script will now exit
	echo Please reset the build variable to run the program
	echo "------------------------------------------------"
	echo
	echo
	exit 1
fi

## Developement (THIS MUST BE YES FOR DAY-2-DAY Usage, set to No in development to see temp files
#REMOVE_TEMP_FILES=No
#REMOVE_TEMP_FILES=Yes


####################################
##
##   Parameterized Code
##
####################################

for INPUT_BAM in `ls *.starAligned.final.bam`
do

# Get the filename base
FILENAME=`basename ${INPUT_BAM} ".bam"`

# Detect if a *.bam.bai file exists, if not determine if you need a symbolic link to an existing *.bai or if you need to create one as samtools idxstats and bedtools multicov requires and index
if [ -e ${INPUT_BAM}.bai ]
then
	echo Found an existing index
else
	#No Index found so check for *.bai file
	if [ -e ${FILENAME}.bai ]
	then
		#Make symbolic link
		echo Found and existing picard index, creating samtools style symbolic link
		ln -s ${FILENAME}.bai ${INPUT_BAM}.bai
	else
		#No Index exists in folder, make a new one with samtools
		echo No index found, creating an index with samtools
		samtools index ${INPUT_BAM}
	fi
fi

# Get total read count
TOTAL_READ_COUNT=`${SAMTOOLS} idxstats ${INPUT_BAM} | grep -v "*" | awk '{sum+=$3} END {print sum}'`
echo "Total Reads = ${TOTAL_READ_COUNT}"

# Generate a table for the counts by gene
### FEATURECOUNTS PARAMETERS
##   -g The attribute to summarize exon counts to gene level counts, default = gene_id
##   -O Allow reads to be assigned to multiple features (ie. counted multiple times)
##   -s Strand of assay ( 0 = unstranded , 1 = stranded, 2 = reversely stranded )
##   -Q Minimum mapping quality for counting (TOPHAT/STAR outputs 1,2,3, 255 with 255 being unique alignments)
##   -T Number of threads
##   -C Do Not Count chimeric reads
##   -p Fragments will be counted instead of reads (if bam is coordinate sorted this will be slow
##   -a Name of the annotation file

echo Starting Feature Counts
date '+%m/%d/%y %H:%M:%S'
${FEATURECOUNTS} -g gene_name -O -s 0 -Q 10 -T 8 -C -a ${IG_GTF} -o temp_featureCounts_Counts.txt ${INPUT_BAM}
date '+%m/%d/%y %H:%M:%S'
echo Finished Coverage Analysis

## Calculate the total number of good mapped reads
ASSIGNED=`grep -w "Assigned" temp_featureCounts_Counts.txt.summary | cut -f2`
NO_FEATURE=`grep -w "Unassigned_NoFeatures" temp_featureCounts_Counts.txt.summary | cut -f2`
echo "Assigned is ${ASSIGNED}"
echo "No Feature ${NO_FEATURE}"
FEATURECOUNT_TOTAL=`echo ${ASSIGNED}+${NO_FEATURE} | bc`

echo "The total read count from featurecounts is ${FEATURECOUNT_TOTAL}"

# Clean-up Results, output columns (Gene, ReadCount, GeneSize, ReadsPerBp, RPKM)
awk 'NR>2' temp_featureCounts_Counts.txt | awk -v var1="$FEATURECOUNT_TOTAL" '{OFMT = "%.12f" ; OFS = "\t" ; if( $7 + 0 != 0) print $1, $7, $6, $7/$6, (1000000000*$7)/(var1*$6) ; else print $1, $7, $6, "0", "0"}' > featureCounts_Counts.txt

# Subset to the IgH Constant Genes
grep -w -f ${IGH_CONSTANT_GENES} featureCounts_Counts.txt > IgH_Constant_Counts.txt
# Subset to the IgH Variable Genes
grep -w -f ${IGH_VARIABLE_GENES} featureCounts_Counts.txt > IgH_Variable_Counts.txt
# Subset to the IgK Constant Genes
grep -w -f ${IGK_CONSTANT_GENES} featureCounts_Counts.txt > IgK_Constant_Counts.txt
# Subset to the IgK Variable Genes
grep -w -f ${IGK_VARIABLE_GENES} featureCounts_Counts.txt > IgK_Variable_Counts.txt
# Subset to the IgL Constant Genes
grep -w -f ${IGL_CONSTANT_GENES} featureCounts_Counts.txt > IgL_Constant_Counts.txt
# Subset to the IgL Variable Genes
grep -w -f ${IGL_VARIABLE_GENES} featureCounts_Counts.txt > IgL_Variable_Counts.txt
# Subset to the Non-Bcell Contamination Genes
grep -w -f ${NON_BCELL_CONTAMINATION_GENELIST} featureCounts_Counts.txt > NonBcell_Contamination_Counts.txt

#Calculate Column totals
TOTAL_IGHC_READS=`awk '{sum+=$2} END {print sum}' IgH_Constant_Counts.txt`
echo "Total IgH Constant Reads = ${TOTAL_IGHC_READS}"

TOTAL_IGHV_READS=`awk '{sum+=$2} END {print sum}' IgH_Variable_Counts.txt`
echo "Total IgH Variable Reads = ${TOTAL_IGHV_READS}"

TOTAL_IGKC_READS=`awk '{sum+=$2} END {print sum}' IgK_Constant_Counts.txt`
echo "Total IgK Constant Reads = ${TOTAL_IGKC_READS}"

TOTAL_IGKV_READS=`awk '{sum+=$2} END {print sum}' IgK_Variable_Counts.txt`
echo "Total IgK Variable Reads = ${TOTAL_IGKV_READS}"

TOTAL_IGLC_READS=`awk '{sum+=$2} END {print sum}' IgL_Constant_Counts.txt`
echo "Total IgL Constant Reads = ${TOTAL_IGLC_READS}"

TOTAL_IGLV_READS=`awk '{sum+=$2} END {print sum}' IgL_Variable_Counts.txt`
echo "Total IgL Variable Reads = ${TOTAL_IGLV_READS}"

# Get the Locus totals
TOTAL_IGH=`grep HEAVY_Locus featureCounts_Counts.txt | cut -f2`
echo "Total IgH Locus Reads = ${TOTAL_IGH}"

TOTAL_IGK=`grep KAPPA_Locus featureCounts_Counts.txt | cut -f2`
echo "Total IgK Locus Reads = ${TOTAL_IGK}"

TOTAL_IGL=`grep LAMBDA_Locus featureCounts_Counts.txt | cut -f2`
echo "Total IgL Locus Reads = ${TOTAL_IGL}"

TOTAL_IG=`echo ${TOTAL_IGH}+${TOTAL_IGK}+${TOTAL_IGL} | bc`
PERCENT_IG=`echo " scale=5 ; ${TOTAL_IG}/${FEATURECOUNT_TOTAL}" | bc`
TOTAL_LIGHT_CHAIN=`echo ${TOTAL_IGK}+${TOTAL_IGL} | bc`
TOTAL_LIGHT_VARIABLE=`echo ${TOTAL_IGKV_READS}+${TOTAL_IGLV_READS} | bc`
TOTAL_LIGHT_CONSTANT=`echo ${TOTAL_IGKC_READS}+${TOTAL_IGLC_READS} | bc`

PERCENT_KAPPA=`echo " scale=5 ; ${TOTAL_IGK}/${TOTAL_LIGHT_CHAIN}" | bc`
PERCENT_LAMBDA=`echo " scale=5 ; ${TOTAL_IGL}/${TOTAL_LIGHT_CHAIN}" | bc`

# Calculate geomean of the contamination gene RPKM
echo Launching R
R <${CALCULATE_GEOMEAN_R} --no-save
echo Finished R
# This will export a file with a single column/row with the geometric mean
GEOMEAN=`cut -f1 geomean.txt | head -n1`

# Create Title Line
echo "Percent Ig = ${PERCENT_IG} ; Kappa/(K+L) = ${PERCENT_KAPPA} ; Lambda/(K+L) = ${PERCENT_LAMBDA} ; Non B Contamination = ${GEOMEAN}" > title.txt

# Create tables for graphing
#  OFMT = "%.12f"  - This prevents awk from printing in scientific notation
echo -e CommonName"\t"Count"\t"Percentage"\t"TotalFrequency"\t"Locus"\t"ElementSize > Header.txt
awk -v var1="$TOTAL_IGHV_READS" -v var2="$FEATURECOUNT_TOTAL" '{ OFMT = "%.12f" ; OFS = "\t" ; if( $2 + 0 != 0) print $1, $2, $2/var1, $2/var2, "IGHV", $3 ; else print $1, $2, "0", "0", "IGHV", $3}' IgH_Variable_Counts.txt > temp_IgHV_Calculation.txt
awk -v var1="$TOTAL_IGHC_READS" -v var2="$FEATURECOUNT_TOTAL" '{ OFMT = "%.12f" ; OFS = "\t" ; if( $2 + 0 != 0) print $1, $2, $2/var1, $2/var2, "IGHC", $3 ; else print $1, $2, "0", "0", "IGHC", $3}' IgH_Constant_Counts.txt > temp_IgHC_Calculation.txt

awk -v var1="$TOTAL_LIGHT_VARIABLE" -v var2="$FEATURECOUNT_TOTAL" '{ OFMT = "%.12f" ; OFS = "\t" ; if( $2 + 0 != 0) print $1, $2, $2/var1, $2/var2, "IGKV", $3 ; else print $1, $2, "0", "0", "IGKV", $3}' IgK_Variable_Counts.txt > temp_IgKV_Calculation.txt
awk -v var1="$TOTAL_LIGHT_CONSTANT" -v var2="$FEATURECOUNT_TOTAL" '{ OFMT = "%.12f" ; OFS = "\t" ; if( $2 + 0 != 0) print $1, $2, $2/var1, $2/var2, "IgLC", $3 ; else print $1, $2, "0", "0", "IgLC", $3}' IgK_Constant_Counts.txt > temp_IgKC_Calculation.txt

awk -v var1="$TOTAL_LIGHT_VARIABLE" -v var2="$FEATURECOUNT_TOTAL" '{ OFMT = "%.12f" ; OFS = "\t" ; if( $2 + 0 != 0) print $1, $2, $2/var1, $2/var2, "IGLV", $3 ; else print $1, $2, "0", "0", "IGLV", $3}' IgL_Variable_Counts.txt > temp_IgLV_Calculation.txt
awk -v var1="$TOTAL_LIGHT_CONSTANT" -v var2="$FEATURECOUNT_TOTAL" '{ OFMT = "%.12f" ; OFS = "\t" ; if( $2 + 0 != 0) print $1, $2, $2/var1, $2/var2, "IgLC", $3 ; else print $1, $2, "0", "0", "IgLC", $3}' IgL_Constant_Counts.txt > temp_IgLC_Calculation.txt

#Create Graphing tables
cat Header.txt temp_IgHC_Calculation.txt temp_IgHV_Calculation.txt > Graph_IgH.txt
cat Graph_IgH.txt

cat Header.txt temp_IgKC_Calculation.txt temp_IgKV_Calculation.txt temp_IgLC_Calculation.txt temp_IgLV_Calculation.txt > Graph_IgL.txt
cat Graph_IgL.txt

#Graph out the IgH distribution
echo Launching R
R <${IGH_GRAPH_R} --no-save
echo Finished R


# Record results to table
#Determine if table exists or not
if [ -e ${RESULTS_TABLE} ]
then
	echo Results Table Exists
else
	echo Results Table Does Not Exist, It will be created
	echo -e Sample"\t"PrimaryIgHC"\t"PrimaryIgHC_Freq"\t"SecondaryIgHC"\t"DeltaIgHC"\t"PrimaryIgHV"\t"PrimaryIgHV_Freq"\t"SecondaryIgHV"\t"DeltaIgHV"\t"PrimaryIgLC"\t"PrimaryIgLC_Freq"\t"SecondaryIgLC"\t"DeltaIgLC"\t"PrimaryIgLV"\t"PrimaryIgLV_Freq"\t"SecondaryIgLV"\t"DeltaIgLV"\t"TOTAL_IGHC_READS"\t"TOTAL_IGHV_READS"\t"TOTAL_IGKC_READS"\t"TOTAL_IGKV_READS"\t"TOTAL_IGLC_READS"\t"TOTAL_IGLV_READS"\t"TOTAL_IGH"\t"TOTAL_IGK"\t"TOTAL_IGL"\t"TOTAL_IG"\t"PERCENT_IG"\t"TOTAL_LIGHT_CHAIN"\t"TOTAL_LIGHT_VARIABLE"\t"TOTAL_LIGHT_CONSTANT"\t"PERCENT_KAPPA"\t"PERCENT_LAMBDA"\t"Top1"\t"Top2"\t"Mean_Top_Delta"\t"NonB_Contamination > ${RESULTS_TABLE}

fi

#Calculate Results and Add to table
PrimaryIgHC=`sort -nrk3 temp_IgHC_Calculation.txt | head -n1 | cut -f1`
PrimaryIgHC_Freq=`sort -nrk3 temp_IgHC_Calculation.txt | head -n1 | cut -f3`
SecondaryIgHC=`sort -nrk3 temp_IgHC_Calculation.txt | head -n2 | tail -n1 | cut -f1`
SecondaryIgHC_Freq=`sort -nrk3 temp_IgHC_Calculation.txt | head -n2 | tail -n1 | cut -f3`
DeltaIgHC=`echo " scale=5 ; ${PrimaryIgHC_Freq}-${SecondaryIgHC_Freq}" | bc`

PrimaryIgHV=`sort -nrk3 temp_IgHV_Calculation.txt | head -n1 | cut -f1`
PrimaryIgHV_Freq=`sort -nrk3 temp_IgHV_Calculation.txt | head -n1 | cut -f3`
SecondaryIgHV=`sort -nrk3 temp_IgHV_Calculation.txt | head -n2 | tail -n1 | cut -f1`
SecondaryIgHV_Freq=`sort -nrk3 temp_IgHV_Calculation.txt | head -n2 | tail -n1 | cut -f3`
DeltaIgHV=`echo " scale=5 ; ${PrimaryIgHV_Freq}-${SecondaryIgHV_Freq}" | bc`

cat temp_IgKC_Calculation.txt temp_IgLC_Calculation.txt > forTable_IgLC.txt
cat temp_IgKV_Calculation.txt temp_IgLV_Calculation.txt > forTable_IgLV.txt

PrimaryIgLC=`sort -nrk3 forTable_IgLC.txt | head -n1 | cut -f1`
PrimaryIgLC_Freq=`sort -nrk3 forTable_IgLC.txt | head -n1 | cut -f3`
SecondaryIgLC=`sort -nrk3 forTable_IgLC.txt | head -n2 | tail -n1 | cut -f1`
SecondaryIgLC_Freq=`sort -nrk3 forTable_IgLC.txt | head -n2 | tail -n1 | cut -f3`
DeltaIgLC=`echo " scale=5 ; ${PrimaryIgLC_Freq}-${SecondaryIgLC_Freq}" | bc`

PrimaryIgLV=`sort -nrk3 forTable_IgLV.txt | head -n1 | cut -f1`
PrimaryIgLV_Freq=`sort -nrk3 forTable_IgLV.txt | head -n1 | cut -f3`
SecondaryIgLV=`sort -nrk3 forTable_IgLV.txt | head -n2 | tail -n1 | cut -f1`
SecondaryIgLV_Freq=`sort -nrk3 forTable_IgLV.txt | head -n2 | tail -n1 | cut -f3`
DeltaIgLV=`echo " scale=5 ; ${PrimaryIgLV_Freq}-${SecondaryIgLV_Freq}" | bc`

echo -e $PrimaryIgHC_Freq"\t"$DeltaIgHC > Top.txt
echo -e $PrimaryIgHV_Freq"\t"$DeltaIgHV >> Top.txt
echo -e $PrimaryIgLC_Freq"\t"$DeltaIgLC >> Top.txt
echo -e $PrimaryIgLV_Freq"\t"$DeltaIgLV >> Top.txt

Top1=`sort -nrk1 Top.txt | head -n1 | cut -f1`
Top2=`sort -nrk1 Top.txt | head -n2 | tail -n1 | cut -f1`
Top1_Delta=`sort -nrk1 Top.txt | head -n1 | cut -f2`
Top2_Delta=`sort -nrk1 Top.txt | head -n2 | tail -n1 | cut -f2`
Mean_Top_Delta=`echo " scale=5 ; (${Top1_Delta}+${Top2_Delta})/2" | bc`

# Write data to summary table
echo -e ${FILENAME}"\t"$PrimaryIgHC"\t"$PrimaryIgHC_Freq"\t"$SecondaryIgHC"\t"$DeltaIgHC"\t"$PrimaryIgHV"\t"$PrimaryIgHV_Freq"\t"$SecondaryIgHV"\t"$DeltaIgHV"\t"$PrimaryIgLC"\t"$PrimaryIgLC_Freq"\t"$SecondaryIgLC"\t"$DeltaIgLC"\t"$PrimaryIgLV"\t"$PrimaryIgLV_Freq"\t"$SecondaryIgLV"\t"$DeltaIgLV"\t"$TOTAL_IGHC_READS"\t"$TOTAL_IGHV_READS"\t"$TOTAL_IGKC_READS"\t"$TOTAL_IGKV_READS"\t"$TOTAL_IGLC_READS"\t"$TOTAL_IGLV_READS"\t"$TOTAL_IGH"\t"$TOTAL_IGK"\t"$TOTAL_IGL"\t"$TOTAL_IG"\t"$PERCENT_IG"\t"$TOTAL_LIGHT_CHAIN"\t"$TOTAL_LIGHT_VARIABLE"\t"$TOTAL_LIGHT_CONSTANT"\t"$PERCENT_KAPPA"\t"$PERCENT_LAMBDA"\t"$Top1"\t"$Top2"\t"$Mean_Top_Delta"\t"$GEOMEAN >> ${RESULTS_TABLE}

# Write data to local file
echo -e Sample"\t"PrimaryIgHC"\t"PrimaryIgHC_Freq"\t"SecondaryIgHC"\t"DeltaIgHC"\t"PrimaryIgHV"\t"PrimaryIgHV_Freq"\t"SecondaryIgHV"\t"DeltaIgHV"\t"PrimaryIgLC"\t"PrimaryIgLC_Freq"\t"SecondaryIgLC"\t"DeltaIgLC"\t"PrimaryIgLV"\t"PrimaryIgLV_Freq"\t"SecondaryIgLV"\t"DeltaIgLV"\t"TOTAL_IGHC_READS"\t"TOTAL_IGHV_READS"\t"TOTAL_IGKC_READS"\t"TOTAL_IGKV_READS"\t"TOTAL_IGLC_READS"\t"TOTAL_IGLV_READS"\t"TOTAL_IGH"\t"TOTAL_IGK"\t"TOTAL_IGL"\t"TOTAL_IG"\t"PERCENT_IG"\t"TOTAL_LIGHT_CHAIN"\t"TOTAL_LIGHT_VARIABLE"\t"TOTAL_LIGHT_CONSTANT"\t"PERCENT_KAPPA"\t"PERCENT_LAMBDA"\t"Top1"\t"Top2"\t"Mean_Top_Delta"\t"NonB_Contamination > ${FILENAME}.purityResults.txt
echo -e ${FILENAME}"\t"$PrimaryIgHC"\t"$PrimaryIgHC_Freq"\t"$SecondaryIgHC"\t"$DeltaIgHC"\t"$PrimaryIgHV"\t"$PrimaryIgHV_Freq"\t"$SecondaryIgHV"\t"$DeltaIgHV"\t"$PrimaryIgLC"\t"$PrimaryIgLC_Freq"\t"$SecondaryIgLC"\t"$DeltaIgLC"\t"$PrimaryIgLV"\t"$PrimaryIgLV_Freq"\t"$SecondaryIgLV"\t"$DeltaIgLV"\t"$TOTAL_IGHC_READS"\t"$TOTAL_IGHV_READS"\t"$TOTAL_IGKC_READS"\t"$TOTAL_IGKV_READS"\t"$TOTAL_IGLC_READS"\t"$TOTAL_IGLV_READS"\t"$TOTAL_IGH"\t"$TOTAL_IGK"\t"$TOTAL_IGL"\t"$TOTAL_IG"\t"$PERCENT_IG"\t"$TOTAL_LIGHT_CHAIN"\t"$TOTAL_LIGHT_VARIABLE"\t"$TOTAL_LIGHT_CONSTANT"\t"$PERCENT_KAPPA"\t"$PERCENT_LAMBDA"\t"$Top1"\t"$Top2"\t"$Mean_Top_Delta"\t"$GEOMEAN >> ${FILENAME}.purityResults.txt

# Rename output files and cleanup all the temp files

if [ $REMOVE_TEMP_FILES == "Yes" ]
then
	rm forTable_IgLC.txt
	rm forTable_IgLV.txt
	rm Graph_IgL.txt
	rm Graph_IgH.txt
	rm IgL_Constant_Counts.txt
	rm temp_IgLC_Calculation.txt
	rm IgL_Variable_Counts.txt
	rm temp_IgLV_Calculation.txt
	rm IgK_Constant_Counts.txt
	rm temp_IgKC_Calculation.txt
	rm IgK_Variable_Counts.txt
	rm temp_IgKV_Calculation.txt
	rm IgH_Constant_Counts.txt
	rm temp_IgHC_Calculation.txt
	rm IgH_Variable_Counts.txt
	rm temp_IgHV_Calculation.txt
	rm title.txt
	rm temp_featureCounts_Counts.txt
	rm featureCounts_Counts.txt
	rm Rplots.pdf
	rm Header.txt
	rm temp_featureCounts_Counts.txt.summary
	rm Top.txt
	rm NonBcell_Contamination_Counts.txt
	rm geomean.txt
fi
mv IGH.png ${FILENAME}_IgH.png
mv IGL.png ${FILENAME}_IgL.png

# Record results to table

done
