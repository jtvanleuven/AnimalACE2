library(stringr)

seqs <- read.table("results/jalview_blc_aligment.txt", stringsAsFactors = F)
names <- seqs$V1[which(str_detect(seqs$V1,">"))]
seqs <- seqs$V1[which(!str_detect(seqs$V1,">"))]
seqs <- seqs[!nchar(seqs) < length(names)]
aling_length <- length()
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
names(seqs.tab.t.short) <- as.character(1:length(man_numbs))
seqs.tab.t.short <- rbind(conserv.short,seqs.tab.t.short)


#read in phylogenetic info
#google sheets file that has been curated
tax_info <- read.csv("results/cdhit_clusters_nodups - 100620_drops_removed.csv")
matchup <- data.frame(name=row.names(seqs.tab.t), seqs = NA)
for(i in 1:length(matchup$name)){
  tmp <- paste(as.character(as.matrix(seqs.tab.t)[i,]),collapse="")
  matchup$seqs[i] <- str_replace_all(tmp, "-","")
}
table(matchup$seqs %in% tax_info$seq)  ##make sure using sequence as index will work
tax_info_order <- tax_info[match(matchup$seqs, tax_info$seq), ]

#join tax info with homo sapien-matched index table
seqs.tab.t.short$kingdom <- c(NA, tax_info_order$kingdom)
seqs.tab.t.short$phylum <- c(NA, tax_info_order$phylum)
seqs.tab.t.short$class <- c(NA, tax_info_order$class)
seqs.tab.t.short$order <- c(NA, tax_info_order$order)
seqs.tab.t.short$family <- c(NA, tax_info_order$family)
seqs.tab.t.short$genus <- c(NA, tax_info_order$genus)
seqs.tab.t.short$species <- c(NA, tax_info_order$sci)
write.csv(seqs.tab.t.short, "results/google_sheet_homo_sapien_aas.csv", quote = F)

