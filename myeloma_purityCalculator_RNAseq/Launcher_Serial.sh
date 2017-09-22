#!/bin/sh

#######################################################################
#######################################################################
## RNAseq based Myeloma Purity Calculator
##
## Copyright (c) 2017 Translational Genomics Research Institute
##
## This software may be modified and distributed under the terms
## of the MIT license.  See the LICENSE file for details.
##
## Major Contributors: Jonathan J. Keats 
## Minor Contributors:
#######################################################################
#######################################################################


PURITY_CHECK_SCRIPT=/home/jkeats/toolkit_jkeats/myelomaPurityChecker/purityChecker.sh
START_DIR=`pwd`

# Add for loop to enter project folders based on a list
for line in `cat $1`
do

echo
echo
echo "--------------------------------"
echo Processing patient ${line} folder

# Enter the project folder
cd ${line}

# Capture Primary Directory
PRIMARY_DIR=`pwd`

if [ -d "TSMRU" ]
then
    # The TSMRU Directory Exists
    echo Found a TSMRU folder
    # Enter the TSMRU folder and determine how many specimens exist
    cd TSMRU

    FOLDER_COUNT=`ls -d */ | rev | cut -c 2- | rev | wc -l`

    echo "Found ${FOLDER_COUNT} folders in the TSMRU folder"

    # Create a list of the TSMRU folders detected
    ls -d */ | rev | cut -c 2- | rev > Folder_List
    # Total number of exome folder
    TOTAL_EXOME_FOLDERS=`wc -l Folder_List | awk '{print $1}'`
    echo Found a total of ${TOTAL_EXOME_FOLDERS} TSMRU folders

    # Loop through each folder and launch each analysis job

	for rna in `cat Folder_List`
	do
	    echo Starting purity check on ${rna}
	    INT_DIR_A=`pwd`
	    cd ${rna}/${rna}.starDir
	    # Submit job to qsub
	    ${PURITY_CHECK_SCRIPT}
	    cd ${INT_DIR_A}
	    echo
	done

    cd ${PRIMARY_DIR}
    touch Purity_Check_Complete

else

    echo WARNING - NO TSMRU FOLDER, SKIPPING

fi

cd ${START_DIR}
echo

done

echo
echo
echo "--------------------------------"
echo Finished Starting All Jobs
echo
echo
echo
echo
