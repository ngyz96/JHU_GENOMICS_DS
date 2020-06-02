library(AnnotationHub)

ahub = AnnotationHub()
ahub = subset(ahub, species == 'Homo sapiens')
qhs = query(ahub, "CpG Islands")
cpg = qhs[[1]]
#q1 26641
autosomes <- paste0('chr', 1:22)
table(cpg@seqnames %in% autosomes)
#q2 1031
table(cpg@seqnames == 'chr4')
#q3 41135164
qhs = query(ahub, c('H3K4me3', 'H1', "EpigenomeRoadMap", "E003"))
data <- qhs[[2]]
data = subset(data, seqnames %in% autosomes)
sum(width(data))
#q4 4.770728
qhs = query(ahub, c('H3K27me3', 'H1', "EpigenomeRoadMap", "E003"))
qhs
data2 <- qhs[[2]]
data2 = subset(data2, seqnames %in% autosomes)
mean(data2$signalValue)
#q5 10289096
sum(width(intersect(data,data2)))
#6 0.5383644
bivalent <- intersect(data, data2)
cpg = subset(cpg, seqnames %in% autosomes)
length(subsetByOverlaps(bivalent, cpg)) / length(bivalent)
#7 0.241688
sum(width(intersect(bivalent, cpg))) / sum(width(cpg))
#8 9774777
cpg10k <- resize(cpg, width=20000, fix='center')
sum(width(intersect(bivalent,cpg10k)))
#9 0.006768009
sum(width(cpg)) / 3000000000
#10 176.1415
inout <- matrix(0, ncol=2,nrow=2)
colnames(inout) = c('in', 'out')
rownames(inout) = c('in', 'out')
inout[1,1] <- sum(width(intersect(cpg, bivalent, ignore.strand=T)))
inout[1,2] <- sum(width(setdiff(cpg, bivalent, ignore.strand=T)))
inout[2,1] <- sum(width(setdiff(bivalent, cpg, ignore.strand=T)))
inout[2,2] <- 3e9 - sum(inout)
oddsR <- inout[1,1] * inout [2,2] / (inout[1,2] * inout[2,1])
