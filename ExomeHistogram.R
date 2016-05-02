library(plyr);
library(ggplot2);
library(data.table);

# USER INPUT
setwd("/data/users/nate/GenomePlot");
NUMBEROFBINS = 1000; # e.g. 1000 bins corresponds to approx. 3 billion / 1000 = 3 million bp per bin.

GenericStarHeader = read.csv("GenericStarHeader.sam", sep="\t");
ChromosomeLengths = as.numeric(gsub("LN:", "", GenericStarHeader$SO.coordinate))[1:24];
CumulativeLengths = c(0, cumsum(ChromosomeLengths)[1:23]);
GenomeSize = sum(ChromosomeLengths);
BinLengths = GenomeSize / NUMBEROFBINS;
breaks = c(seq(1, GenomeSize, BinLengths), GenomeSize); 

ExonCountsOriginal = read.csv("ACC2ExonCounts.csv");
ExonCountsOriginal$chr = gsub("X", "23", ExonCountsOriginal$chr);
ExonCountsOriginal$chr = gsub("Y", "24", ExonCountsOriginal$chr);
ExonCountsOriginal$chr = as.integer(ExonCountsOriginal$chr);
ExonCountsOriginal$start = unlist(mclapply(1:24, function(x) ExonCountsOriginal$start[ExonCountsOriginal$chr==x] + CumulativeLengths[x], mc.cores=12));
ExonCountsOriginal$end = unlist(mclapply(1:24, function(x) ExonCountsOriginal$end[ExonCountsOriginal$chr==x] + CumulativeLengths[x], mc.cores=12));

BinID = rep(0, nrow(ExonCountsOriginal));
for (i in 1:(length(breaks)-1)) {
  BinMin = breaks[i];
  BinMax = breaks[i+1];
  BinID[ExonCountsOriginal$start >= BinMin & ExonCountsOriginal$start < BinMax] <- i;
} 

# If a bin contains 0 for all samples.
EmptyBins = data.frame(BinID=setdiff(1:(length(breaks)-1), unique(BinID)));
ZeroCounts = as.data.frame(matrix(0, nrow=nrow(EmptyBins), ncol=(ncol(ExonCountsOriginal)-5)));
colnames(ZeroCounts) = colnames(ExonCountsOriginal)[6:82];
EmptyBins = cbind(EmptyBins, ZeroCounts);

RegionCounts = aggregate(ExonCountsOriginal[ , 6:ncol(ExonCountsOriginal)], by=list(BinID), FUN=sum);
RegionCounts = rename(RegionCounts, c("Group.1"="BinID"));
RegionCounts = rbind(RegionCounts, EmptyBins);
RegionCounts = RegionCounts[order(RegionCounts$BinID), ];
RegionCounts$start = breaks[1:(length(breaks)-1)];
RegionCounts$end = breaks[2:length(breaks)];

RegionCountsFlat = melt(RegionCounts, id=c("BinID", "start", "end"));
RegionCountsFlat = rename(RegionCountsFlat, c("variable"="FileName", "value"="count"));
#RegionCountsFlat$count[RegionCountsFlat$count>10000] <- 10000; # Set a maximum count

# Plot read depth/count.
AnnotateDF = data.frame(Chrom=c(1:22, "X", "Y"));
AnnotateDF$ChromStart= CumulativeLengths;
AnnotateDF$ChromEnd = c(CumulativeLengths[2:length(CumulativeLengths)], GenomeSize);
HighlightDF = AnnotateDF[seq(2, nrow(AnnotateDF), by=2), ];
TextDF = data.frame(Chrom=AnnotateDF$Chrom, xpos=rowMeans(AnnotateDF[, c("ChromStart", "ChromEnd")]));

g = ggplot(RegionCountsFlat, aes(x=start, y=count)) + geom_bar(fill="red", stat="identity");
g = g + facet_grid(FileName ~ ., scales="free_y") + labs(x="");
g = g + geom_rect(data=HighlightDF, inherit.aes=F, aes(xmin=ChromStart, xmax=ChromEnd, ymin=0, ymax=Inf, group=Chrom), fill='yellow', alpha=0.1);
g = g + geom_text(data=TextDF, inherit.aes=F, aes(x=xpos, y=Inf, label=Chrom), size=1, alpha=0.3, colour="blue", vjust="top", hjust="left");
g = g + theme_bw() + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), legend.position="none", panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"));
ggsave("ExomeHist.png", plot=g, width=10, height=49);


