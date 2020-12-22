library(seqinr)
library(stringr)
library(ggplot2)
library(Biostrings)
library(cowplot)
library(ape)
library(genbankr)
library(ggseqlogo)
library(dplyr)
library(reshape2)
library(gggenes)
library(vegan)
library(ggridges)

s_seq_aln <- readDNAStringSet("raw_data/COVID19_S.fst")
#exported alignment at AA from seaview
#still contains a few indels. maybe check lit to see if correct?
s_seq_aln_AA <- readAAStringSet("raw_data/COVID19_S_AA.fst")
sAAmat <- as.matrix(s_seq_aln_AA)
sAAtab <- table(sAAmat)
sAAvars <- vector(mode='list', length = ncol(sAAmat))
sAAvarcnt <- vector(mode='list',length = ncol(sAAmat))
for(i in 1:length(sAAvars)){
  sAAvars[[i]] <- table(sAAmat[,i])
  sAAvarcnt[[i]] <- length(table(sAAmat[,i]))
}
sAAvarcnt.df <- data.frame(sAAvarcnt)
colnames(sAAvarcnt.df) <- 1:length(sAAvarcnt)
plot(t(sAAvarcnt.df))
which(sAAvarcnt.df > 1)

#splot <- as.character(s_seq_aln_AA)
#splot <- as.character(AA)
#ggplot() + geom_logo(splot, seq_type="aa")

#library(pegas)
#s_seq_aln_stings <- as.DNAbin(s_seq_aln)
#nuc.div(s_seq_aln_stings)
#allelicrichness(s_seq_aln_stings)




###load ACE2 data

#SARS-CoV interaction residues: 30-41, 82-84, 353-357
ace <- readxl::read_excel("raw_data/41421_2020_147_MOESM2_ESM.xlsx", skip = 1)  ##data from doi:10.1038/s41421-020-0147-1
ace.all <- data.frame(ace)
ace <- ace[,c("POS", "REF", "ALT", "Function", "Transcript")]
ace$Mut <- str_split(ace$Transcript,"\\:p\\.", simplify = T)[,2]
ace$Mut <- str_split(ace$Mut,"\\/", simplify = T)[,1]
ace$AApos <- str_extract(ace$Mut,"\\d+")
ace <- ace[!ace$Function == "synonymous_variant",]
ace$freq <- rowMeans(ace.all[,9:20], na.rm=T)
ace.pos.count <- as.matrix(table(ace$AApos))
ace.pos.count <- ace.pos.count[order(as.numeric(row.names(ace.pos.count))),]


##plots
ace.pdb <- read.table("raw_data/ace2_interface_residues.pdb")      
s.pdb <- read.table("raw_data/S-RBD_interface_residues.pdb")

ace.bind <- data.frame(unique(ace.pdb$V6))
#ace.gene <- data.frame("molecule"="chrX","gene"="ace2", "start"=15494520, "end"=15602158, "direction"=-1)
ace.gene <- data.frame(molecule="chrX",gene="ACE2", start=1, end=805)
ace.gene[2:6,] <- ace.gene[1,]
ace.gene$bind <- c("motif1","motif2","motif3","motif4","motif5","motif6")
ace.gene$from <- c(26, 79, 324, 330, 353, 393)
ace.gene$to <- c(45, 83, 325, 331, 357, 394)
rbd <- c(19, 24, 27, 28, 30, 31, 33, 34, 35, 37, 38, 41, 42, 45, 79, 82, 83, 321, 322, 323, 324, 325, 326, 327, 329, 330, 353, 354, 355, 356, 357, 383, 386, 387, 389, 393, 555)

s.bind <- data.frame(unique(s.pdb$V6))
s.gene <- data.frame(molecule="COVID19",gene="S", start=1, end=1274)
s.gene[2:4,] <- s.gene[1,]
s.gene$bind <- c("motif1","motif2","motif3","motif4")
s.gene$from <- c(403, 446, 473, 486)
s.gene$to <- c(404, 456, 477, 505)

plot <- rbind(ace.gene, s.gene)
a <- as.numeric(row.names(plot))
s <- which(sAAvarcnt.df > 1)
vlines <- data.frame(x=c(rbd, s), molecule=c(rep("chrX", length(rbd)),rep("COVID19", length(s))))
                     
ggplot(plot, aes(xmin=start, xmax=end, y=molecule, fill=gene)) +
  geom_gene_arrow() +
  scale_fill_brewer(palette="Set3") +
  geom_vline(data=vlines, aes(xintercept = x, color=molecule)) +
  facet_wrap(~molecule, scale="free", ncol=1) +
  geom_subgene_arrow(data=plot, aes(xmin=start, xmax=end, y=molecule, xsubmin=from, xsubmax=to), color="black", fill="black") +
  #geom_segment(xintercept = as.numeric(row.names(ace.plot)), color="red") +
  theme_genes() 
  #theme(axis.title.y = element_blank())



#############################try to make a conservation colored pdb file
##check out list that paul made

data(ggseqlogo_sample)
ggplot() + geom_logo( seqs_aa$AKT1 ) + theme_logo()
names <- read.csv("raw_data/EnsemblSpeciesTable.csv")
#aln <- read.nexus.data("ENSGT00940000158077_gene_tree.nex")
fa <- read.FASTA("ENSGT00940000158077_gene_tree.fa")
for (i in 1:length(fa)){
  string <- paste(as.character(fa[i]),collapse = "")
  fa[i] <- str_replace_all(as.character(string),"-","")
}
tre <- read.tree("ENSGT00940000158077_gene_tree.newick")
aln <- read.alignment("ENSGT00940000158077_gene_tree.msf", format = "msf")

human_num <- which(str_detect(aln$nam, "Hsap"))
human_seq <- aln$seq[[human_num]][1]
seq.dat <- data.frame(str_split(unlist(aln$seq),"",simplify = T), stringsAsFactors = F)
seq.dat.s <- seq.dat[,which(!as.character(seq.dat[human_num,])=="-")]
names(seq.dat.s) <- 1:805
human_seq_str <- str_remove_all(human_seq, "-")
shannon <- matrix(ncol=1, nrow=805)
simpson <- matrix(ncol=1, nrow=805)
aacnt <- matrix(ncol=1, nrow=805)
for(i in 1:nrow(shannon)){
  col <- seq.dat.s[,i]
  col <- col[!col=="-"]
  shannon[i] <- diversity(as.numeric(factor(col)), index = "shannon")
  simpson[i] <- diversity(as.numeric(factor(col)), index = "simpson")
  aacnt[i] <- length(table(col))
}

hace2 <- na.omit(read.table("hace2.pdb", sep = "", comment.char = "", fill=T, header=F, stringsAsFactors = F))  ###watch out becuase this drop the last column of the pdb file
#gpdstring <- read.table("G_penta.pdb", sep = "\t", fill=T, header=F, stringsAsFactors = F)
hace2string <- as.matrix(readLines("hace2.pdb"))
ace.pdb <- read.table("ace2_interface_residues.pdb")
interface_res <- unique(ace.pdb$V6)
###urg, pdb formats are terrible because of the weird spacing. Need to fix
hace2shannon <- hace2string
hace2simpson <- hace2string
hace2div <- hace2string
for (i in 2:nrow(hace2string)-1){
  site <- hace2[i-1,]$V6
  #gpdb2[i,]$V11
  hace2shannon[i,] <- str_replace(hace2shannon[i,],"(\\w+\\s+\\w+\\s+\\w+\\s+\\w+\\s+\\w+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+)\\S+(.*)", paste("\\1", round(shannon[site],digits=3), "\\2",sep='')) 
  hace2simpson[i,] <- str_replace(hace2simpson[i,],"(\\w+\\s+\\w+\\s+\\w+\\s+\\w+\\s+\\w+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+)\\S+(.*)", paste("\\1", round(simpson[site],digits=3), "\\2",sep='')) 
  #print(round(shannon[site],digits=3))
  hace2div[i,] <- str_replace(hace2shannon[i,],"(\\w+\\s+\\w+\\s+\\w+\\s+\\w+\\s+\\w+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+)\\S+(.*)", paste("\\1", aacnt[site], "\\2",sep=''))

}

#write.table(hace2shannon, file="hace2_shannonbeta.pdb", quote = F, row.names = F, col.names = F)
#write.table(hace2div, file="hace2_aacntsbeta.pdb", quote = F, row.names = F, col.names = F)
#write.table(hace2simpson, file="hace2_simpsonbeta.pdb", quote = F, row.names = F, col.names = F)

plot.dat <- data.frame(shannon)
ggplot(data=plot.dat, aes(x=shannon)) +
  geom_histogram() +
  theme_classic() +
  xlab("Amino acid diversity (Shannon Index)") +
  ylab("")

plot.dat2 <- data.frame(simpson)
ggplot(data=plot.data, aes(x=simpson)) +
  geom_histogram() +
  theme_classic() +
  xlab("Amino acid diversity (Simpson Index)") +
  ylab("")

seq.dat.s[115,interface_res]
simpson[interface_res]

seq.dat.heat <- seq.dat.s
row.names(seq.dat.heat) <- aln$nam
seq.dat.heat$name <- row.names(seq.dat.heat) 
seq.dat.heat.l <- melt(seq.dat.heat, id.vars = "name")
seq.dat.heat.l$name <- str_replace(seq.dat.heat.l$name, "\\/.*", "_")
seq.dat.heat.l[which(seq.dat.heat.l$value=="-"),]$value <- NA
ggplot(seq.dat.heat.l, aes(variable, name, fill= value)) + 
  geom_tile() +
  scale_y_discrete(limits=tre$tip.label) +
  scale_fill_discrete(na.value = 'white')

simpson.plot <- data.frame(simpson)
simpson.plot$position <- 1:805
ggplot(simpson.plot, aes(x=position, y=simpson)) +
  geom_vline(xintercept = interface_res, color='red') +
  geom_line()

  
aacnt.plot <- data.frame(aacnt)
aacnt.plot$position <- 1:805
ggplot(aacnt.plot, aes(x=position, y=aacnt)) +
  geom_vline(xintercept = interface_res, color='red') +
  geom_bar(stat = "identity") +
  theme_classic()
    

#run rate4site
#alignment and tree from ENSEMBL  
#rate4site.exe -s ENSGT00940000158077_gene_tree.fa -t ENSGT00940000158077_gene_tree.newick -o ENSGT00940000158077_gene_tree.rate4site
###didn't run. need to troubleshoot



##manually curated list of sequences then aligned with muscle
##extracted human amino acid index
##plot gglogo
rbd <- data.frame(sites = c(19, 24, 27, 28, 30, 31, 33, 34, 35, 37, 38, 41, 42, 45, 79, 82, 83, 321, 322, 323, 324, 325, 326, 327, 329, 330, 353, 354, 355, 356, 357, 383, 386, 387, 389, 393, 555), homo_aa=NA)
jalview <- read.csv("results/google_sheet_homo_sapien_aas.csv")
jalview.plot <- jalview[2:712,2:806]
jalview.plot <- jalview.plot[,rbd$sites]
rdb.seqs <- c(apply(jalview.plot[,1:37], 1, paste, collapse=""))
rbd$homo_aa <- as.character(jalview.plot[which(str_detect(jalview$X, "sapien")),])
rbd$label <- paste(rbd$homo_aa, rbd$sites, sep="")
##update with reduced sequence set
rdb.seqs.nogap <- read.FASTA("results/rbd.ape.fasta", type="AA")
rdb.seqs2 <- vector(length = length(rdb.seqs.nogap), mode = "character")
for(i in 1:length(rdb.seqs.nogap)){
  rdb.seqs2[i] <- paste(as.character(rdb.seqs.nogap)[[i]],collapse="")
}
names(rdb.seqs2) <- as.character(1:length(rdb.seq2))

ggplot() + 
  geom_logo(rdb.seqs2, method = "prob") + 
  scale_x_discrete(limits=as.character(rbd$label)) +
  theme_logo()+
  theme(axis.text.x = element_text(angle=90), axis.text.y=element_blank(), legend.position = "none") +
  ylab('')


#jalview and tax_info_order have same tax
#read in phylogenetic info
jalview.plot.s <- read.csv("results/jalview.plot.s.csv", stringsAsFactors = F, colClasses = c("character"))
names(jalview.plot.s)[11:47] <- as.character(rbd$sites)
jalview.plot.long <- melt(jalview.plot.s, id.vars = c("name","header","kingdom","phylum","class","order","family","genus", "sci", "seq"))
names(jalview.plot.long) <- c("name","header","kingdom","phylum","class","order","family","genus","sci",'seq',"site","aa")
jalview.plot.long[jalview.plot.long$aa=="-",]$aa <- 'NA'
jalview.plot.long[jalview.plot.long$aa=="X",]$aa <- 'NA'
jalview.plot.long <- jalview.plot.long %>%
  arrange(kingdom, phylum, class, order, family, genus)

aa_palette <- data.frame(aa=c("D", "E", "C", "M", "K", "R", "S", "T", "F", "Y", "N", "Q", "G", "L", "V", "I", "A", "W", "H", "P", "NA"),
                         col=c("#E60A0A", "#E60A0A", "#E6E600", "#E6E600", "#145AFF", "#145AFF", "#FA9600", "#FA9600", "#3232AA", "#3232AA", "#00DCDC", "#00DCDC", "#EBEBEB", "#0F820F", "#0F820F", "#0F820F", "#C8C8C8", "#B45AB4", "#8282D2", "#DC9682", "black"))
lookup <- aa_palette$col
names(lookup) <- aa_palette$aa

ggplot(jalview.plot.long, aes(x=as.factor(site), y=header, fill=aa)) +
  geom_tile() +
  scale_y_discrete(limits=unique(jalview.plot.long$header)) +
  scale_fill_manual(values =  unname(lookup[names(table(na.omit(jalview.plot.long$aa)))]), na.value="white")


#####################################evolutionary distance of rbd
dist.man <- read.csv("results/dist2man.csv")
tmp <- table(dist.man$order)
cut <- names(tmp[!tmp<5])
n <- length(cut)
rainbow <- rev(inlmisc::GetTolColors(n, scheme = "smooth rainbow"))
pie(rep(1,n), col = rainbow)
top_orders <- dist.man[dist.man$order %in% cut,]
order.list <- unique(top_orders$order)
top_orders$order <- factor(top_orders$order, levels=order.list)
order.cols <- data.frame(order=order.list, col=rainbow)

ggplot(top_orders, aes(x=val, y=order, fill=order, color=order)) +
  geom_density_ridges() +
  scale_y_discrete(limits=rev(order.list), labels=paste(rev(order.list),' (', as.character(tmp[rev(order.list)]), ')',sep="")) +
  scale_fill_manual(values = rainbow) +
  scale_color_manual(values = rainbow) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x="Corrected AA distance", y="")
#ggsave(filename = "plots/aa_dist_histos.pdf", width = 5.5, height= 7, units = "in")


##reorder aa aligment plot
order.means <- dist.man %>%
  group_by(order) %>%
  summarise(order_mean = mean(val)) %>%
  arrange(order_mean)
#mean may not be the best. order it by max similarity

jalview.plot.long$order <- as.factor(jalview.plot.long$order)
align.p <- ggplot(jalview.plot.long, aes(x=site, y=header, fill=aa)) +
  geom_tile() +
  scale_y_discrete(limits=rev(unique(dist.man$name))) +
  scale_fill_manual(values =  unname(lookup[names(table(na.omit(jalview.plot.long$aa)))])) +
  theme_classic() +
  theme(axis.text.x =  element_text(angle=-90), axis.text.y = element_text(size=5))

row.names(order.cols) <- order.cols$order
#tmp$col <- order.cols[as.character(tmp$order),]$col
ggplot(order.cols, aes(x=1, y=order, fill=order)) +
  geom_tile() +
  scale_y_discrete(limits=rev(order.list)) +
  scale_fill_manual(values = order.cols[names(table(order.cols$order)),]$col) +
  theme(legend.position = "none")

tmp <- jalview.plot.long[jalview.plot.long$site==19,]
tmp$order <- as.character(tmp$order)
col <- tmp$col
order.p <- ggplot(tmp, aes(x=site, y=header, fill=order)) +
  geom_tile() +
  scale_y_discrete(limits=rev(unique(dist.man$name))) +
  scale_fill_manual(values=order.cols[names(table(tmp$order)),]$col, na.value="white") +
  theme_classic() +
  theme(axis.text.y = element_text(size=5), legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Order")


plot_grid(order.p, 
          align.p + 
            theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()),
          nrow = 1, 
          rel_widths = c(1,4),
          align = "h", axis="bt"
          )
ggsave(filename = "plots/alignment_worders.pdf", width = 10, height= 30, units = "in")


##redo plot with just NA animals
jalview.plot.s.animaldist <- read.delim('results/jalview_animalDist.tab', header = T, sep = "\t", stringsAsFactors = F, colClasses = c("character"))
northamerica <- which(str_detect(jalview.plot.s.animaldist$countryDistribution, "United States|Canada|Mexico"))
man <- which(jalview.plot.s.animaldist$header=='man')
plot.w <- jalview.plot.s.animaldist[c(northamerica,man),]
plot.l <- melt(plot.w, id.vars = c("name","header","kingdom","phylum","class","order","family","genus", "sci", "seq", 'biogeographicRealm', 'countryDistribution', 'iucnStatus', 'extinct', 'domestic'))
names(plot.l) <- c("name","header","kingdom","phylum","class","order","family","genus","sci",'seq', 'biogeographicRealm', 'countryDistribution', 'iucnStatus', 'extinct', 'domestic',"site","aa")
plot.l$site <- str_replace(plot.l$site, "X", "")
plot.l$site <- as.numeric(as.character(plot.l$site))
dist.short <- dist.man[dist.man$name %in% plot.l$header,]
align.na.p <- ggplot(plot.l, aes(x=factor(site), y=header, fill=aa)) +
  geom_tile() +
  scale_y_discrete(limits=rev(unique(dist.short$name))) +
  scale_fill_manual(values =  unname(lookup[names(table(na.omit(plot.l$aa)))])) +
  theme_classic() +
  theme(axis.text.x =  element_text(angle=-90), axis.text.y = element_text(size=5))

tmp <- plot.l[plot.l$site==19,]
tmp$order <- as.character(tmp$order)
order.na.p <- ggplot(tmp, aes(x=factor(site), y=header, fill=order)) +
  geom_tile() +
  scale_y_discrete(limits=rev(unique(dist.short$name))) +
  scale_fill_manual(values=order.cols[names(table(tmp$order)),]$col, na.value="white") +
  theme_classic() +
  theme(axis.text.y = element_text(size=5), legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  xlab("Order")

plot_grid(order.na.p, 
          align.na.p + 
            theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank()),
          nrow = 1, 
          rel_widths = c(1,3.5),
          align = "h", axis="bt"
)
ggsave(filename = "plots/alignment_northamerica.pdf", width = 8, height= 8, units = "in")
