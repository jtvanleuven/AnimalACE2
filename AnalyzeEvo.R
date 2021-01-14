library(stringr)
library(ape)
library(phytools)
library(reshape2)
library(dplyr)

seqs <- read.table("results/jalview_blc_aligment.txt", stringsAsFactors = F)
names <- seqs$V1[which(str_detect(seqs$V1,">"))]
seqs <- seqs$V1[which(!str_detect(seqs$V1,">"))]
seqs <- seqs[!nchar(seqs) < length(names)]
#align_length <- length(seqs)
seqs.tab <- data.frame(matrix(, nrow=length(seqs), ncol=length(names)), stringsAsFactors = F)   
names.s <- str_replace_all(names,">", "")
names.s <- str_split(names.s,"/",simplify = T)[,1]
names(seqs.tab) <- names.s
for(i in 1:length(seqs)){
  seqs.tab[i,] <- unlist(str_split(seqs[i],""))
}
seqs.tab.t <- data.frame(t(seqs.tab), stringsAsFactors = F)
row.names(seqs.tab.t) <- names
man_index <- which(str_detect(names,"sapiens"))
row.names(seqs.tab.t)[man_index]
seqs.tab.t[man_index,]
man_numbs <- which(str_detect(seqs.tab.t[man_index,], "[A-Z]"))
seqs.tab.t.short <- seqs.tab.t[,man_numbs]

#conservation score from jalview
#This is an automatically calculated quantitative alignment annotation which measures the number of conserved physico-chemical properties conserved for each column of the alignment. Its calculation is based on the one used in the AMAS method of multiple sequence alignment analysis :
#Livingstone C.D. and Barton G.J. (1993), Protein Sequence Alignments: A Strategy for the Hierarchical Analysis of Residue Conservation. CABIOS Vol. 9 No. 6 (745-756)).
conserv <- read.csv("results/jalview_conservation.csv", header = F, row.names = 1, stringsAsFactors = F)
conserv.short <- conserv[,man_numbs]
names(conserv.short) <- as.character(1:length(man_numbs))
plot.conserv <- data.frame(score=t(conserv.short), type='All residues')
add <- data.frame(score=t(conserv.short[rbd$sites]), type='S-binding residues')
plot.conserv <- rbind(plot.conserv, add)

ggplot(plot.conserv, aes(Conservation, fill=type, colour=type)) +
  geom_density(alpha=0.2) +
  #geom_histogram(bins = 10, alpha=0.5, position='stack') +
  theme_classic() +
  labs(x='Per-site variation', y='') +
  scale_fill_nejm() +
  scale_color_nejm()

ggplot(plot.conserv, aes(x=Conservation, y=type, fill=type, colour=type)) +
  geom_violin(alpha=0.9, adjust = 1.2) +
  theme_classic() +
  labs(x='Per-site variation', y='') +
  scale_fill_nejm() +
  scale_color_nejm() +
  theme(legend.position = 'none')


#read in phylogenetic info
#google sheets file that has been curated
tax_info <- read.csv("results/cdhit_clusters_nodups - 100620_drops_removed.csv")
tax_info <- read.csv("results/cdhit_clusters_nodups-jeremy.csv")
matchup <- data.frame(name=row.names(seqs.tab.t), seqs = NA)
for(i in 1:length(matchup$name)){
  tmp <- paste(as.character(as.matrix(seqs.tab.t)[i,]),collapse="")
  matchup$seqs[i] <- str_replace_all(tmp, "-","")
}
table(matchup$seqs %in% tax_info$seq)  ##make sure using sequence as index will work
tax_info_order <- tax_info[match(matchup$seqs, tax_info$seq), ]
#write.csv(tax_info_order, "results/tax_info_order.csv", quote = F, row.names = F)

#join tax info with homo sapien-matched index table
seqs.tab.t.short$kingdom <- c(NA, tax_info_order$kingdom)
seqs.tab.t.short$phylum <- c(NA, tax_info_order$phylum)
seqs.tab.t.short$class <- c(NA, tax_info_order$class)
seqs.tab.t.short$order <- c(NA, tax_info_order$order)
seqs.tab.t.short$family <- c(NA, tax_info_order$family)
seqs.tab.t.short$genus <- c(NA, tax_info_order$genus)
seqs.tab.t.short$species <- c(NA, tax_info_order$sci)
#write.csv(seqs.tab.t.short, "results/google_sheet_homo_sapien_aas.csv", quote = F)

rbd <- data.frame(sites = c(19, 24, 27, 28, 30, 31, 33, 34, 35, 37, 38, 41, 42, 45, 79, 82, 83, 321, 322, 323, 324, 325, 326, 327, 329, 330, 353, 354, 355, 356, 357, 383, 386, 387, 389, 393, 555), homo_aa=NA)
jalview <- read.csv("results/google_sheet_homo_sapien_aas.csv")
jalview.plot <- jalview[2:712,2:806]
jalview.plot <- jalview.plot[,rbd$sites]
get <- which(tax_info_order$name %in% tax_info$name)
jalview.plot.s <- jalview.plot[get,]
names(jalview.plot.s) <- rbd$sites
tax_info_order <- na.omit(tax_info[match(matchup$seqs, tax_info$seq), ])
jalview.plot.s <- cbind(tax_info_order[,c("sci", "name","header","kingdom","phylum","class","order","family","genus", "seq")], jalview.plot.s)
jalview.plot.long <- melt(jalview.plot.s, id.vars = c("name","header","kingdom","phylum","class","order","family","genus", "sci", "seq"))
names(jalview.plot.long) <- c("name","header","kingdom","phylum","class","order","family","genus","sci",'seq',"site","aa")
jalview.plot.long[jalview.plot.long$aa %in% c("-","X"),]$aa <- NA
jalview.plot.long <- jalview.plot.long %>%
  arrange(kingdom, phylum, class, order, family, genus)

rdb.seqs <- c(apply(jalview.plot.s[,as.character(rbd$sites)], 1, paste, collapse=""))
rdb.seqs <- str_replace_all(rdb.seqs, "X", "-")
##some of these have bad gaps. just get rid of them
gap.seqs <- which(str_detect(rdb.seqs, "---"))
rdb.seqs.nogap <- rdb.seqs[-gap.seqs]
jalview.plot.s <- jalview.plot.s[-gap.seqs,]
write.csv(jalview.plot.s, "results/jalview.plot.s.csv", quote = F, row.names = F)

#write.FASTA(as.AAbin(as.list(rdb.seqs.nogap)), file="results/rbd.ape.fasta", header = tax_info_order[-gap.seqs,]$header)
rdb.seqs.nogap <- read.FASTA("results/rbd.ape.fasta", type="AA")
dat <- read.phyDat("results/rbd.ape.fasta", type="AA", format = "fasta")
#modelTest(dat, model=c("JTT", "LG", "WAG"))
# Model   df     logLik      AIC          AICw AICc AICcw      BIC
# 1      JTT 1349 -10460.253 23618.51 5.502489e-283  NaN   NaN 25791.64
# 2    JTT+I 1350 -10366.036 23432.07 1.676318e-242  NaN   NaN 25606.81
# 3    JTT+G 1350  -9817.332 22334.66  3.337175e-04  NaN   NaN 24509.40
# 4  JTT+G+I 1351  -9814.464 22330.93  2.162118e-03  NaN   NaN 24507.28
# 5       LG 1349 -10540.781 23779.56 5.856980e-318  NaN   NaN 25952.69
# 6     LG+I 1350 -10443.842 23587.68 2.713058e-276  NaN   NaN 25762.42
# 7     LG+G 1350  -9841.945 22383.89  6.826035e-15  NaN   NaN 24558.63
# 8   LG+G+I 1351  -9840.162 22382.32  1.493519e-14  NaN   NaN 24558.67
# 9      WAG 1349 -10393.601 23485.20 4.866972e-254  NaN   NaN 25658.33
# 10   WAG+I 1350 -10306.867 23313.73 8.338060e-217  NaN   NaN 25488.47
# 11   WAG+G 1350  -9811.057 22322.11  1.772759e-01  NaN   NaN 24496.85
# 12 WAG+G+I 1351  -9808.525 22319.05  8.202283e-01  NaN   NaN 24495.40
##sort of stupid. alignment is only 37 aa long.
dist.mat <- dist.ml(rdb.seqs.nogap, "WAG")
treeNJ <- NJ(dist.mat)
dist.wide <- as.data.frame(as.matrix(dist.mat))
names(dist.wide) <- tax_info_order[-gap.seqs,]$header
dist.man <- data.frame(name=tax_info_order[-gap.seqs,]$header, val=dist.wide[,"man"], stringsAsFactors = F)
dist.man <- cbind(dist.man, tax_info_order[-gap.seqs,c('sci', 'common', 'kingdom', 'phylum','class','order','family','genus')])
dist.man <- dist.man %>%
  arrange(val,kingdom, phylum, class, order, family, genus)
#write.csv(dist.man, file="results/dist2man.csv", quote = F, row.names = F)



##north american animal species
#from https://www.mammaldiversity.org/
animals <- read.csv("raw_data/mdd.csv")
animals$sciName <- str_replace(animals$sciName, "_", " ")
na.animals <- animals[str_detect(animals$countryDistribution, 'States|Canada'),]
na.animals <- na.animals[!is.na(na.animals$sciName),]
row.names(na.animals) <- na.animals$sciName
nrow(na.animals)

animals.s <- animals[,c('sciName',"biogeographicRealm", 'countryDistribution', 'iucnStatus', 'extinct', 'domestic')]
jalview.plot.s.animaldist <- left_join(jalview.plot.s, animals.s, by=c('sci' = 'sciName'))
write.table(jalview.plot.s.animaldist, "results/jalview_animalDist.tab", quote = F, row.names = F, sep="\t")

