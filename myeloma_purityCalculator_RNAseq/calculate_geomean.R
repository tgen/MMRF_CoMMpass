# Load the psych library to calculate a geomean
library(psych)

# Read in the non-B cell contamination table
NONB <- read.table("NonBcell_Contamination_Counts.txt", header=F)

# Calcualte the geomean and print output to new file with no more than 3 decimal places
gm <- geometric.mean(NONB$V5+1,na.rm=TRUE)
rounded <- round(gm, digits=3)
write.table(rounded, file="geomean.txt", sep="\t", row.names=F, col.names=F)

# Exit
q()


