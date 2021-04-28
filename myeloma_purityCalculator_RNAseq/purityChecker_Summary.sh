#!/usr/bin/env bash

## Usage help message

usage()
{
  # `cat << EOF` This means that cat should stop reading when EOF is detected
  cat << EOF

  Usage: $0 [ options ]

  -h      Display Help

  Required Options
  -s      Study name [Test] - Basename for output files and graphs

EOF
# EOF is found above and hence cat command stops reading. This is equivalent to echo but much neater when printing out.
  exit 2
}

## Assign default values
STUDY=Test

## Capture and assign variables from inputs
while getopts 's:?h' flag
do
    case ${flag} in
        s) STUDY=${OPTARG};;
        h|?) usage;;
    esac
done

## Send variable information to log
echo The assigned study name is: $STUDY

### From a jetstream project directory generate a summary matrix file

# Create a template file with header lines
echo -e Sample"\t"PrimaryIgHC"\t"PrimaryIgHC_Freq"\t"SecondaryIgHC"\t"DeltaIgHC"\t"PrimaryIgHV"\t"PrimaryIgHV_Freq"\t"SecondaryIgHV"\t"DeltaIgHV"\t"PrimaryIgLC"\t"PrimaryIgLC_Freq"\t"SecondaryIgLC"\t"DeltaIgLC"\t"PrimaryIgLV"\t"PrimaryIgLV_Freq"\t"SecondaryIgLV"\t"DeltaIgLV"\t"TOTAL_IGHC_READS"\t"TOTAL_IGHV_READS"\t"TOTAL_IGKC_READS"\t"TOTAL_IGKV_READS"\t"TOTAL_IGLC_READS"\t"TOTAL_IGLV_READS"\t"TOTAL_IGH"\t"TOTAL_IGK"\t"TOTAL_IGL"\t"TOTAL_IG"\t"PERCENT_IG"\t"TOTAL_LIGHT_CHAIN"\t"TOTAL_LIGHT_VARIABLE"\t"TOTAL_LIGHT_CONSTANT"\t"PERCENT_KAPPA"\t"PERCENT_LAMBDA"\t"Top1"\t"Top2"\t"Mean_Top_Delta"\t"NonB_Contamination > ${STUDY}.puritySummary.txt

# Now find the results files and add their result lines to the summary file
find . -name "*.purityResults.txt" -exec awk 'NR>1' {} >> ${STUDY}.puritySummary.txt \;

# Now plot results
