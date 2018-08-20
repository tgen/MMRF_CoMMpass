

args = commandArgs(TRUE)
hetsTSV = args[1]

dat <- read.table(hetsTSV, header = TRUE)

all_modes <- function(x) {
  mode <- NULL
  for ( CN in 2:(length(x)-1)) {
    if (( x[CN] > x[CN-1]) & (x[CN] > x[CN+1])) {
      mode <- c(mode,CN)
    } 
    if (( x[CN] < x[CN-1]) & (x[CN] < x[CN+1])) {
      mode <- c(mode,CN)
    } 
  }
  if (length(mode) == 0) {
    mode = 'Fail'
  }
  return(mode)
}

modes_values <- all_modes(density(dat$Fold.Change, adjust = 1.4)$y)
density(dat$Fold.Change, adjust = 1.4)$y[modes_values]
density(dat$Fold.Change, adjust = 1.4)$x[modes_values]
density(dat$Fold.Change, adjust = 1.4)$x[modes_values[2]]

if (length(density(dat$Fold.Change, adjust = 1.4)$x[modes_values]) == 1) {
  STATUS <- data.frame()
  write.table(STATUS, file='CNA_unimodal_fail', col.names = FALSE)
  MAINPEAK <- density(dat$Fold.Change, adjust = 1.4)$x[modes_values[1]]
  CUTOFF <- MAINPEAK+.3
  CUTOFF2 <- MAINPEAK-.3
  POSITIONS <- dat[ dat$Fold.Change > CUTOFF | dat$Fold.Change < CUTOFF2 ,]
  write.table(POSITIONS, file='CNA_positions_to_exclude.txt', col.names = FALSE, sep = "\t")
} else {
  STATUS <- data.frame()
  write.table(STATUS, file='CNA_unimodal_fail', col.names = FALSE)
  densum <- sum(density(dat$Fold.Change, adjust = 1.4)$y)
  modepos <- density(dat$Fold.Change, adjust = 1.4)$y[modes_values[1]]
  percent <- modepos/densum
  if (percent > 0.00176) {
    MAINPEAK <- density(dat$Fold.Change, adjust = 1.4)$x[modes_values[1]]
    CUTOFF <- MAINPEAK+.3
    CUTOFF2 <- MAINPEAK-.3
    POSITIONS <- dat[ dat$Fold.Change > CUTOFF | dat$Fold.Change < CUTOFF2 ,]
    write.table(POSITIONS, file='CNA_positions_to_exclude.txt', col.names = FALSE, sep = "\t")
  } else {
    INCVAR=1
    PASS=0
    for (peak in density(dat$Fold.Change, adjust = 1.4)$y[modes_values]) {
      print(peak)
      percent <- peak/densum
      print(percent)
      if (percent > 0.00176 & INCVAR < length(density(dat$Fold.Change, adjust = 1.4)$x[modes_values])) {
        MAINPEAK <- density(dat$Fold.Change, adjust = 1.4)$x[modes_values[INCVAR]]
        CUTOFF <- MAINPEAK+.3
        CUTOFF2 <- MAINPEAK-.3
        POSITIONS <- dat[ dat$Fold.Change > CUTOFF | dat$Fold.Change < CUTOFF2 ,]
        write.table(POSITIONS, file='CNA_positions_to_exclude.txt', col.names = FALSE, sep = "\t")
        break
      } else if (percent > 0.00176 & INCVAR == length(density(dat$Fold.Change, adjust = 1.4)$x[modes_values])) {
        MAINPEAK <- density(dat$Fold.Change, adjust = 1.4)$x[modes_values[INCVAR]]
        CUTOFF <- MAINPEAK+.3
        CUTOFF2 <- MAINPEAK-.3
        POSITIONS <- dat[ dat$Fold.Change > CUTOFF | dat$Fold.Change < CUTOFF2 ,]
        write.table(POSITIONS, file='CNA_positions_to_exclude.txt', col.names = FALSE, sep = "\t")
        break
      } else {
        INCVAR=INCVAR + 1
      }
    }
  }
}

png(filename = paste(hetsTSV, ".density.png", sep=""), width = 648, height = 368, units = "mm", res = 300)
plot(density(dat$Fold.Change, adjust = 1.4)) + abline(v=MAINPEAK+0.3) + abline(v=MAINPEAK-0.3) + abline(v=MAINPEAK)
dev.off()


