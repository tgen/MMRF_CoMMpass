library(tidyverse)

# Read in datafile
purity <- read_tsv("CLARION.puritySummary.txt")
purity <- read_tsv("purityChecker.txt")


# Plot result summary
ggplot(purity, aes(Sample, Top1, ymin = Top2, ymax = Top1, fill = Mean_Top_Delta)) +
  geom_crossbar(stat='identity', linetype='blank', width=0.75) +
  geom_hline(aes(yintercept=0.9), colour="blue", linetype="dashed") +
  scale_x_discrete(name="Samples Tested") +
  scale_y_continuous(name = "Top Two Isoform Range",
                     limits=c(0,1),
                     breaks=seq(0, 1, 0.1)) +
  scale_fill_gradient2("Mean\nInter\nIsoform\nDelta",
                       low="blue", mid = "yellow", high="red",
                       limits=c(0.0,1.0),
                       midpoint = 0.5) + 
  theme(axis.text.x=element_blank(),
        axis.title = element_text(size=16, vjust=0.4))

# Plot XY results
ggplot(purity, aes(Top1, Top2, color=NonB_Contamination)) + 
  geom_point() + 
  scale_color_gradientn(colors = rainbow(6))

ggplot(purity, aes(Mean_Top_Delta, NonB_Contamination)) + 
  geom_point()

ggplot(purity, aes(PrimaryIgHC_Freq, PrimaryIgHV_Freq, color=NonB_Contamination)) + 
  geom_point() + 
  scale_color_gradientn(colors = rainbow(6))

ggplot(purity, aes(PrimaryIgLC_Freq, PrimaryIgLV_Freq, color=NonB_Contamination)) + 
  geom_point() + 
  scale_color_gradientn(colors = rainbow(6))

# histogram
ggplot(purity, aes(NonB_Contamination)) + 
  geom_histogram(binwidth = 0.2)

mean(purity$NonB_Contamination)
sd(purity$NonB_Contamination)
median(purity$NonB_Contamination)
quantile(purity$NonB_Contamination,0.9)
quantile(purity$NonB_Contamination,0.95)
quantile(purity$NonB_Contamination,0.99)
