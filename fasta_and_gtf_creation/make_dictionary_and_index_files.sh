#!/bin/bash

##
## make_dictionary_and_index_files.sh
##
## Created by Austin Christofferson on 06/12/2017.
## Copyright 2017 Translational Genomics Research Institute. All rights reserved.
##
## Description - Creates the dictionary and index files for the custom reference 
##		 genome and gtf used in the tgen05 recipe of the jetstream pipeline at TGen.
##
## Usage = Make_gtf_and_related_tables.sh [--help -h] [-p /path/to/picard-tools -s /path/to/samtools]
##
## Options:
##      -p      Full path to the picard tools directory containing CreateSequenceDictionary.jar
##      -s      Full path to the samtools directory containing samtools
##
####################################################################
####################################################################

# Set help options

if [ "$*" == "--help" ] || [ "$*" == "-h" ] || [ "$*" == "" ]
then
        grep "^##" $0
        echo
        exit 0
fi

# Set command line options

while getopts ":p:s:" opt
do
        case $opt in
                p)
                        PICARD_PATH=$OPTARG ;;
                s)
                        SAMTOOLS_PATH=$OPTARG ;;
                :)
                        echo "Option -$OPTARG requires an argument." >&2
                        exit 1 ;;
                \?)
                        echo "Invalid option: -$OPTARG" >&2 ;;

        esac
done

# 7) Create dict and fai files using picard tools
        java -jar ${PICARD_PATH}/CreateSequenceDictionary.jar R=hs37d5_plusRibo_plusOncoViruses_plusERCC.fa O=hs37d5_plusRibo_plusOncoViruses_plusERCC.dict

        ${SAMTOOLS_PATH}/samtools faidx hs37d5_plusRibo_plusOncoViruses_plusERCC.fa


