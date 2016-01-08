library(VennDiagram)
# Count_1 = 53
# Count_2 = 2
# Count_3 = 8032
# Count_4 = 4487
# Count_5 = 2
# Count_12 = 3
# Count_13 = 2249
# Count_14 = 309
# Count_15 = 2
# Count_23 = 3
# Count_24 = 3
# Count_25 = 0
# Count_34 = 6733
# Count_35 = 2
# Count_45 = 2
# Count_123 = 10
# Count_124 = 10
# Count_125 = 0
# Count_134 = 7722
# Count_135 = 2
# Count_145 = 2
# Count_234 = 8
# Count_235 = 0
# Count_245 = 0
# Count_345 = 6
# Count_1234 = 2
# Count_1235 = 0
# Count_1245 = 0
# Count_1345 = 1
# Count_2345 = 0
# Count_12345 = 0
Count_1 = 53
Count_2 = 2
Count_3 = 8032
Count_4 = 4487
Count_5 = 2
Count_12 = 2
Count_13 = 39
Count_14 = 36
Count_15 = 1
Count_23 = 2
Count_24 = 2
Count_25 = 0
Count_34 = 671
Count_35 = 1
Count_45 = 1
Count_123 = 2
Count_124 = 2
Count_125 = 0
Count_134 = 36
Count_135 = 0
Count_145 = 0
Count_234 = 2
Count_235 = 0
Count_245 = 0
Count_345 = 1
Count_1234 = 2
Count_1235 = 0
Count_1245 = 0
Count_1345 = 0
Count_2345 = 0
Count_12345 = 0
# generate pdf
draw.quintuple.venn(area1=Count_1, area2=Count_2, area3=Count_3, area4=Count_4, area5 = Count_5, n12=Count_12, n13=Count_13, n14=Count_14, n15=Count_15, n23=Count_23, n24=Count_24, n25=Count_25, n34=Count_34, n35=Count_35, n45=Count_45, n123=Count_123, n124=Count_124, n125=Count_125, n134=Count_134, n135=Count_135, n145=Count_145, n234=Count_234, n235=Count_235, n245=Count_245, n345=Count_345, n1234=Count_1234, n1235=Count_1235, n1245=Count_1245, n1345=Count_1345, n2345=Count_2345, n12345=Count_12345, category = c('SeqCNV','CoNIFER','CNVnator','CNVer','XHMM'), lwd=1, lty="blank", col=c('red','green','blue','orange','orchid3'), cat.cex = 1.2, margin = 0.1, fill=c('red','green','blue','orange','orchid3'), cat.col=c('red','green','blue','orange','orchid3'), reverse = FALSE, scaled=TRUE, ind = TRUE, cex = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.2, 1, 1.2, 1, 1.2, 1, 1.2, 1, 1.2, 1, 1.2, 0.9, 1.2, 0.9, 1.2, 0.9, 1.2, 0.9, 1.2, 0.9, 1.2, 1.2, 1.2, 1.2, 1.2, 1.8))
venn.plot <- draw.quintuple.venn(area1=Count_1, area2=Count_2, area3=Count_3, area4=Count_4, area5 = Count_5, n12=Count_12, n13=Count_13, n14=Count_14, n15=Count_15, n23=Count_23, n24=Count_24, n25=Count_25, n34=Count_34, n35=Count_35, n45=Count_45, n123=Count_123, n124=Count_124, n125=Count_125, n134=Count_134, n135=Count_135, n145=Count_145, n234=Count_234, n235=Count_235, n245=Count_245, n345=Count_345, n1234=Count_1234, n1235=Count_1235, n1245=Count_1245, n1345=Count_1345, n2345=Count_2345, n12345=Count_12345, category = c('SeqCNV','CoNIFER','CNVnator','CNVer','XHMM'), lwd=1, lty="blank", col=c('red','green','blue','orange','orchid3'), cat.cex = 1.2, margin = 0.1, fill=c('red','green','blue','orange','orchid3'), cat.col=c('red','green','blue','orange','orchid3'), reverse = FALSE, scaled=TRUE, ind = TRUE, cex = c(1.8, 1.8, 1.8, 1.8, 1.8, 1.2, 1, 1.2, 1, 1.2, 1, 1.2, 1, 1.2, 1, 1.2, 0.9, 1.2, 0.9, 1.2, 0.9, 1.2, 0.9, 1.2, 0.9, 1.2, 1.2, 1.2, 1.2, 1.2, 1.8))
# Writing to file
tiff(filename = "venn.tiff", compression = "none", antialias = "none", width=800, height=800);
grid.draw(venn.plot);
dev.off();
