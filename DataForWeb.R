library(stringr)
library(ggplot2)

master <- read.csv('../../website/data/master_table.tsv', sep='\t')
dists <- read.csv('results/dist2man.csv')
dists <- dists[master$Unique.ID,]
master$dist2human <- dists$name

rdb.seqs.nogap <- read.FASTA("results/rbd.ape.fasta", type="AA")
dat <- read.phyDat("results/rbd.ape.fasta", type="AA", format = "fasta")
dist.mat <- dist.ml(rdb.seqs.nogap, "WAG")
treeNJ <- NJ(dist.mat)
dist.wide <- as.data.frame(as.matrix(dist.mat))
master$dist2human <- dist.wide[master$Unique.ID,653]
write.table(master, file='../../website/data/master_table.tsv', sep = '\t', quote = F, row.names = F)
