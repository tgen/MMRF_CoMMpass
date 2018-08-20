#Start

library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

hist <- read.table("GenderDataTable.txt", header=T)

ggplot(hist, aes(x=Normal_Ratio)) + geom_histogram(aes(y=..density..), binwidth = 0.01, fill = "magenta") + geom_density(alpha=.2) + scale_y_continuous(limits=c(0,10))


