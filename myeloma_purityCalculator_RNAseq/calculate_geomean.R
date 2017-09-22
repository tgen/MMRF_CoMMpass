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

# Load the psych library to calculate a geomean
#library(psych)
if (!require('psych')) install.packages('psych'); library('psych')


# Read in the non-B cell contamination table
NONB <- read.table("NonBcell_Contamination_Counts.txt", header=F)

# Calcualte the geomean and print output to new file with no more than 3 decimal places
gm <- geometric.mean(NONB$V5+1,na.rm=TRUE)
rounded <- round(gm, digits=3)
write.table(rounded, file="geomean.txt", sep="\t", row.names=F, col.names=F)

# Exit
q()


