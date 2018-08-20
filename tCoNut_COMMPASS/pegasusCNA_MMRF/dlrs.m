function dlrsStat = dlrs(log2Ratios)
%
% derivative log-ratio spread

xdiffSTD = std(diff(log2Ratios));
dlrsStat = xdiffSTD/sqrt(2);
