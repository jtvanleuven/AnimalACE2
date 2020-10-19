library(seqinr)
library(stringr)
library(ggplot2)
library(Biostrings)
library(cowplot)
library(ape)
library(genbankr)
library(ggseqlogo)
library(tidyr)
library(vegan)
library(data.table)
library(rentrez)
library(taxize)

###ACE2 sequences
##get data from a variety of sources and combine them. 
##438 from blast search with human ACE2 followed by 90% length filter
##239 from Ensembl ACE2 orthologs:
##318 from NCBI ACE2 orthologs: https://www.ncbi.nlm.nih.gov/gene/59272/ortholog/?scope=89593
##410 from Damas et al. 2020 PNAS paper: https://doi.org/10.1073/pnas.2010146117
##63 additional bat accesion numbers from search done by Jeremy Ellis

##seqs from blast search with 90% length coverage cutoff
cov90 <- read.alignment("raw_data/ACE2_animals_90percoverage.fst", format="fasta")
cov90names <- data.frame(readLines("raw_data/ACE2_animals_90percoverage.names"), stringsAsFactors = F)
cov90names[,2] <- str_remove(str_split(cov90names[,1], "\\[", simplify = T)[,2],"\\]")
cov90names$acc_num <- str_split(cov90names[,1]," ", simplify = T)[,1]
cov90names$acc_num <- str_replace(cov90names$acc_num, ">", "")
cov90names <- cov90names[,c("V2","acc_num")]
cov90_fa <- read.fasta("raw_data/ACE2_animals__90percoverage.fasta", as.string = T)
#common names
cov90_common <- read.csv("raw_data/blast90_common_names.csv", sep=",")

##seqs from ensemble ortholog search
ensembl_names <- read.csv("raw_data/EnsemblSpeciesTable.csv")
ensembl_fa <- read.fasta("raw_data/ENSGT00940000158077_gene_tree.fa", as.string = T)
for (i in 1:length(ensembl_fa)){
  string <- paste(as.character(ensembl_fa[i]),collapse = "")
  ensembl_fa[i] <- str_replace_all(as.character(string),"-","")
}
ensembl_seq_tab <- data.frame(abb=str_sub(str_extract(names(ensembl_fa), "_\\w\\w\\w\\w\\/"), 2, -2),
                              seq=toupper(unlist(ensembl_fa)))
str_sub(str_extract(names(ensembl_fa), "_\\w\\w\\w\\w\\/"), 2, -2)
ensembl_names$abbr <- paste(str_extract(str_split(ensembl_names$Scientific.name, " ", simplify = T)[,1], "\\w"),
                            str_extract(str_split(ensembl_names$Scientific.name, " ", simplify = T)[,2], "\\w\\w\\w"), sep = "")
abbr <- str_split(str_split(names(ensembl_fa), "\\/", simplify = T)[,1],"_", simplify=T)[,2]
write.fasta(sequences = as.list(ensembl_seq_tab$seq), names = as.list(row.names(ensembl_seq_tab)), file.out = "raw_data/ENSGT00940000158077.fa")
ensembl_seq_tab$name <- "NA"
run <- F
if(run){
  for(i in 1:nrow(ensembl_seq_tab)){
    term <- str_split(row.names(ensembl_seq_tab)[i],"_",simplify = T)[,1]
    tmp <- entrez_search(db="gene", term=term)
    if(tmp$count > 0){
      tmp_sum <- entrez_summary(db="gene", id=tmp$ids)
      ensembl_seq_tab[i,]$name <- tmp_sum$organism$scientificname
    }
  }
  write.csv(ensembl_seq_tab, "raw_data/EnsemblSeqTab.csv", quote = F)
}
#ensembl_names$seq <- ensembl_fa[ensembl_names$Taxon.ID]
#some problems with the names. some species have same three-letter abbreviation. Hard to parse out in could not easily fix 


##pulled ACE2 orthologs from NCBI refseq
ncbi_fa <- read.fasta("raw_data/ACE2_refseq_protein.fasta", as.string = T)
ncbi_names <- fread("raw_data/ACE2_refseq_protein.fasta", header = F, stringsAsFactors = F)
ncbi_names <- ncbi_names[str_detect(ncbi_names$V1, ">"),] 
ncbi_names$V1 <- str_replace(ncbi_names$V1, ">", "")
ncbi_names$scientific <- str_split(ncbi_names$V1,"\\[", simplify = T)[,2]
ncbi_names$acc_num <- str_split(ncbi_names$V1," ", simplify = T)[,1]
ncbi_names$scientific <- str_replace(ncbi_names$scientific,"\\]","")
ncbi_names <- ncbi_names[,c("scientific","acc_num")]
ncbi_names$seq <- "NA"
for(i in 1:length(ncbi_fa)){
  ncbi_names[i,]$seq <-  toupper(as.character(unlist(ncbi_fa[i])))
}
#write.csv(ncbi_names, file="raw_data/ACE2_refseq_table.csv", quote = F, row.names = F)


###we searched for additional ACE2 sequence from bats. Here is what we found
bats <- c("ABU54053.1","ALJ94034.1","ADN93471.1","QJF77811.1","ADN93475.1","ALJ94035.1","QJF77803.1","QJF77801.1","QJF77802.1","QJF77791.1","BAF50705.1","ADN93477.1","QJF77796.1","QJF77795.1","QJF77789.1","QJF77815.1","QJF77840.1","QJF77809.1","QJF77806.1","QJF77804.1","QJF77810.1","QJF77841.1","QJF77816.1","QJF77807.1","QJF77822.1","QJF77829.1","QJF77823.1","QJF77842.1","QJF77821.1","QJF77824.1","QJF77825.1","QJF77830.1","QJF77831.1","QJF77834.1","QJF77827.1","QJF77835.1","QJF77808.1","QJF77794.1","QJF77836.1","QJF77792.1","QJF77826.1","QJF77828.1","QJF77793.1","QJF77813.1","ACT66266.1","QJF77832.1","QJF77797.1","QJF77833.1","QJF77819.1","QJF77805.1","QJF77814.1","QJF77820.1","QJF77839.1","QJF77837.1","QJF77798.1","QJF77799.1","QJF77790.1","QJF77812.1","QJF77838.1"," QJF77800.1","QJF77817.1","ADN93476.1","ADN93473.1")
run <- F
if(run){   ###slow step
  bats_seqs <- data.frame(matrix(ncol=3,nrow=length(bats)))
  names(bats_seqs) <- c("name","locus_tag","seq")
  for(i in 1:length(bats)){
    fetch <- entrez_fetch(db="protein", id = bats[i], rettype = "fasta")
    bats_seqs[i,]$name <- str_split(str_split(fetch,"\\[", simplify = T)[,2],"\\]",simplify=T)[,1]
    bats_seqs[i,]$locus_tag <- bats[i]
    bats_seqs[i,]$seq <- str_replace_all(str_split(fetch,"\\]", simplify = T)[,2],"\n","")
  }
  write.csv(bats_seqs,file="raw_data/batsTable.csv",quote = F, row.names = F)
  write.fasta(sequences = as.list(bats_seqs$seq), names=paste(bats_seqs$locus_tag, bats_seqs$name, sep=" "), file.out = "raw_data/bats.fa", nbchar=60)
}
bats_seqs <- read.csv("raw_data/batsTable.csv", stringsAsFactors = F)

##need to extract sequences from accession numbers.
damas <- read.csv("raw_data/pnas.2010146117.sd01.csv", stringsAsFactors = F, header = T)
damas$seq <- NA
run <- F
if(run){  ##slow step
  for(i in 1:nrow(damas)){
    id <- str_split(damas[i,]$Data.availability," ",simplify = T)[,1]
    if(str_detect(damas[i,]$Data.availability, "protein")){
      fetch <- entrez_fetch(db="protein", id = id, rettype = "fasta")
      damas[i,]$seq <- str_replace_all(str_split(fetch,"\\]", simplify = T)[,2],"\n","")
    }
    if(str_detect(damas[i,]$Data.availability, "gene")){
      fetch <- entrez_fetch(db="nucleotide", id = id, rettype = "fasta")
      seq <- str_replace_all(str_split(fetch,"cds", simplify = T)[,2],"\n","")
      seq <- str_replace_all(seq, "N", "")
      damas[i,]$seq <- as.character(translate(DNAString(seq)))
    }
  }
}
#write.csv(damas, "raw_data/damas_table.csv", quote = F, row.names = F)
damas <- read.csv("raw_data/damas_table.csv")
damas_short <- damas[!is.na(damas$seq),]
#write.fasta(sequences=as.list(damas_short$seq), names = damas_short$Data.availability, as.string = T, file.out="raw_data/damas.fa")



##what is the overlap??
species <- list(damas=damas_short$Species, bats=bats_seqs$name, refseq=ncbi_names$scientific, blast=cov90names$V2, ensembl=ensembl_names$Scientific.name)
library(UpSetR)
upset(fromList(species), order.by = "freq")


##Try another way: what is the overlap if we use the sequences??
species_seqs <- list(damas=damas_short$seq, bats=bats_seqs$seq, refseq=ncbi_names$seq, blast=unname(toupper(unlist(cov90_fa))), ensembl=ensembl_seq_tab$seq)
upset(fromList(species_seqs), order.by = "freq")
upset_tab <- fromList(species_seqs)
allseq<- data.frame(seq = unname(unlist(species_seqs)),
                    length=NA, name=names(unlist(species_seqs)))
for(i in 1:length(allseq$seq)){
  allseq[i,]$length <- length(str_split(allseq[i,]$seq,"", simplify = T))
}

#cat ACE2_animals__90percoverage.fasta ACE2_refseq_protein.fasta bats.fa ENSGT00940000158077.fa damas.fa > all.fa
#./cd-hit -i ~/GoogleDrive/projects/SarsCov2/analysis/AnimalACE2/raw_data/all.fa -o all_cdhit.fa -c 1 -d 30

cdhit.aa <- read.fasta("raw_data/all_cdhit.fa", as.string = T)
cdhit <- read.csv("raw_data/all_cdhit.fa.clstr", sep="\t", fill = T, row.names=NULL, header = F)
###get names of cdhit groups.
names_all <- read.csv("raw_data/all_names.csv", sep=",")
names_all$acc_num <- str_extract(names_all$acc_num, "[A-Z0-9_]+")
names_all$str<-"NA"
for(i in 1:nrow(names_all)){
  names_all[i,]$str <- paste(names_all[i,]$scientific, names_all[i,]$acc_num, collapse = ",")
}
for(i in 1:nrow(cdhit)){
  if(str_detect(cdhit[i,]$V1, "Cluster")){
    id <- str_extract(cdhit[i,]$V1,">Cluster\\s\\d+")
  }else{
    cdhit[i,]$V1 <- id
  }
}
cdhit <- cdhit[str_detect(cdhit$V2, "\\w"),]
cdhit$V1 <- str_replace(cdhit$V1,">","")
cdhit$V3 <- str_split(cdhit$V2, ",", simplify = T)[,1]
cdhit$V3 <- str_replace(cdhit$V3, "aa", "")
cdhit$V2 <- str_extract(cdhit$V2, ">[A-Z0-9_]+")
cdhit$V2 <- str_replace(cdhit$V2, ">", "")
cdhit$sci <- "NA"
cdhit$common <- "NA"
cdhit$kingdom <- "NA"
cdhit$phylum <- "NA"
cdhit$class <- "NA"
cdhit$order <- "NA"
cdhit$family <- "NA"
cdhit$genus <- "NA"
run <- F
if(run){
  for(i in 1:nrow(cdhit)){
    cdhit[i,]$sci <- unique(names_all[str_detect(names_all$str, cdhit[i,]$V2),]$scientific)
    if(!is.na(cdhit[i,]$sci)){
      tax_class <- classification(cdhit[i,]$sci, db = 'ncbi', return_id = F)
      if(!is.na(tax_class[1])){
        if("kingdom" %in% tax_class[[1]]$rank){
          cdhit[i,]$kingdom <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="kingdom")]
        }
        if("phylum" %in% tax_class[[1]]$rank){
          cdhit[i,]$phylum <- tax_class[[i]]$name[which(tax_class[[1]]$rank=="phylum")]
        }
        if("class" %in% tax_class[[1]]$rank){
          cdhit[i,]$class <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="class")]
        }
        if("family" %in% tax_class[[1]]$rank){
          cdhit[i,]$family <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="family")]
        }
        if("genus" %in% tax_class[[1]]$rank){
          cdhit[i,]$genus <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="genus")]
        }
      }
      x <- sci2comm(sci=cdhit[i,]$sci, db='itis',simplify = T)
      cdhit[i,]$common <- x[[1]][1]
    }
  }
  cdhit$header <- cdhit$common
  cdhit[is.na(cdhit$header),]$header <- cdhit[is.na(cdhit$header),]$sci
  cdhit[is.na(cdhit$header),]$header <- cdhit[is.na(cdhit$header),]$V2
  write.csv(cdhit, "results/cdhit_clusters.csv", quote = F, row.names = F)
}
cdhit <- read.csv("results/cdhit_clusters.csv")

###merge clusters with sequence info
cdhit.orig <- read.csv("raw_data/all_cdhit.fa.clstr", sep="\t", fill = T, row.names=NULL, header = F)
for(i in 1:nrow(cdhit.orig)){
  if(str_detect(cdhit.orig[i,]$V1, "Cluster")){
    id <- str_extract(cdhit.orig[i,]$V1,">Cluster\\s\\d+")
  }else{
    cdhit.orig[i,]$V1 <- id
  }
}
#cdhit.orig[which(!str_detect(cdhit.orig$V2, "aa")),]$V2 <- cdhit.orig[which(!str_detect(cdhit.orig$V2, "aa"))+1,]$V2
cdhit.orig <- cdhit.orig[str_detect(cdhit.orig$V2,"\\.\\.\\. \\*"),]
cdhit.orig$V1 <- cdhit.orig$V1 <- str_replace(cdhit.orig$V1,">","")
cdhit.orig$V2 <- str_split(cdhit.orig$V2, ", ", simplify = T)[,2]
cdhit.orig$V2 <- str_replace(cdhit.orig$V2, "\\.\\.\\. \\*", "")
#cdhit.orig$V2 <- str_extract(cdhit.orig$V2, ">[A-Z0-9_]+")
cdhit.orig$V2 <- str_replace(cdhit.orig$V2, ">", "")
names(cdhit.orig) <- c("cluster", "name")
for (i in 1:length(cdhit.aa)){
  string <- toupper(paste(as.character(cdhit.aa[i]),collapse = ""))
  cdhit.aa[i] <- str_replace_all(as.character(string),"-","")
}
cdhit.aa.tab <- data.frame(name=names(cdhit.aa), seq=unlist(cdhit.aa), stringsAsFactors = F)
#cdhit.aa.tab$name <- str_extract(cdhit.aa.tab$name, "[A-Z0-9_]+")
cdhit.aa.tab <- merge(cdhit.orig, cdhit.aa.tab, by='name')

names(cdhit) <- c("cluster", "name", "length", "sci", "common", "kingdom", "phylum", "class", "family", "genus", "header", "order")
cdhit_short <- cdhit.aa.tab
cdhit_short$sci <- "NA"
cdhit_short$common <- "NA"
cdhit_short$kingdom <- "NA"
cdhit_short$phylum <- "NA"
cdhit_short$class <- "NA"
cdhit_short$order <- "NA"
cdhit_short$family <- "NA"
cdhit_short$genus <- "NA"
cdhit_short$header <- "NA"
cdhit_short$length <- "NA"
for(i in 1:nrow(cdhit_short)){
  cdhit_short[i,]$sci <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$sci[1]
  cdhit_short[i,]$common <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$common[1]
  cdhit_short[i,]$kingdom <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$kingdom[1]
  cdhit_short[i,]$phylum <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$phylum[1]
  cdhit_short[i,]$class <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$class[1]
  cdhit_short[i,]$order <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$order[1]
  cdhit_short[i,]$family <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$family[1]
  cdhit_short[i,]$genus <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$genus[1]
  cdhit_short[i,]$header <- cdhit[which(cdhit$cluster == cdhit_short[i,]$cluster),]$header[1]
  cdhit_short[i,]$length <- nchar(cdhit_short[i,]$seq)
}
cdhit[which(str_detect(cdhit$sci, "Homo")),]  #human seq is 805aa long
cdhit_short <- cdhit_short[which(cdhit_short$length > 805*.9 & cdhit_short$length < 805*1.1),]
#number of unique species
length(unique(cdhit_short$header))

write.fasta(sequences=as.list(cdhit_short$seq), names = str_replace(cdhit_short$header," ","_"), as.string = T, file.out = "results/combined_filtered.faa")
write.csv(cdhit_short, "results/cdhit_clusters_nodups.csv", quote = F, row.names = F)
rodents <- which(cdhit_short$order == "Rodentia")
write.fasta(sequences=as.list(cdhit_short[rodents,]$seq), names = str_replace(cdhit_short[rodents,]$header, " ", "_"), as.string = T, file.out = "results/combined_filtered_rodents.faa")


###opened cdhit_short in google docs. hand editing list by 
#####1. Fixing ensembl gene names
#####2. remove duplicate sequences
#####3. incorporate ACE2 f##authors of Damas et al 2020 shared sequences with me


damas_fa <- read.fasta("../../data/damas/one_PROTEIN_perSPS.410species.4APR2020.fasta", as.string = T)
#just overwrite the sequences that I pulled from ncbi for damas data
for (i in 1:length(damas_fa)){
  damas_fa[i]  <- toupper(str_replace_all(as.character(damas_fa[i]),"-",""))
}

add_damas <- damas[is.na(damas$seq),]
#add_damas <- damas
for(i in 1:nrow(add_damas)){
  species <- str_replace_all(add_damas$Species[i]," ","_")
  get <- which(str_detect(names(damas_fa), species))
  if(length(get) > 1){
    print(paste("multiple hits for ",species,sep = ""))
  }else if(length(get) == 1){
    add_damas$seq[i] <- as.character(damas_fa[get])
  }
}

add_damas$sci <- "NA"
add_damas$kingdom <- "NA"
add_damas$phylum <- "NA"
add_damas$class <- "NA"
add_damas$order <- "NA"
add_damas$family <- "NA"
add_damas$genus <- "NA"
run <- F
if(run){
  for(i in 1:nrow(add_damas)){
    if(!is.na(add_damas[i,]$Species)){
      tax_class <- classification(add_damas[i,]$Species, db = 'itis', return_id = F)
      if(!is.na(tax_class[1])){
        if("kingdom" %in% tax_class[[1]]$rank){
          add_damas[i,]$kingdom <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="kingdom")]
        }
        if("phylum" %in% tax_class[[1]]$rank){
          add_damas[i,]$phylum <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="phylum")]
        }
        if("class" %in% tax_class[[1]]$rank){
          add_damas[i,]$class <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="class")]
        }
        if("family" %in% tax_class[[1]]$rank){
          add_damas[i,]$family <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="family")]
        }
        if("genus" %in% tax_class[[1]]$rank){
          add_damas[i,]$genus <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="genus")]
        }
      }
    }
  }
  add_damas$order <- add_damas$Order
  add_damas$header <- add_damas$Common.name
  add_damas[is.na(add_damas$header),]$header <- add_damas[is.na(add_damas$header),]$sci
  add_damas[is.na(add_damas$header),]$header <-add_damas[is.na(add_damas$header),]$V2
  add_damas_short <- add_damas[,c("Data.availability", "seq", "Species", "Common.name", "kingdom", "phylum", "class", "order", "family", "genus", "header")]
  write.csv(add_damas_short, "results/add_damas.csv", quote = F, row.names = F)
}



#I couldn't get the R import function to work with google sheets, so it must be exported as .csv occassionaly
googlesheet <- read.csv("results/cdhit_clusters_nodups - cdhit_clusters_nodups.csv", stringsAsFactors = F)
length(unique(googlesheet$sci))

##after looking up names by had, I went back and found taxonomy info
run <- F  ###only needs run once
if(run){
  for(i in 1:nrow(googlesheet)){
    if(is.na(googlesheet[i,]$genus)){
      tax_class <- classification(googlesheet[i,]$sci, db = 'ncbi', return_id = F)
      if(!is.na(tax_class[1])){
        if("kingdom" %in% tax_class[[1]]$rank){
          googlesheet[i,]$kingdom <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="kingdom")]
        }
        if("phylum" %in% tax_class[[1]]$rank){
          googlesheet[i,]$phylum <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="phylum")]
        }
        if("class" %in% tax_class[[1]]$rank){
          googlesheet[i,]$class <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="class")]
        }
        if("family" %in% tax_class[[1]]$rank){
          googlesheet[i,]$family <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="family")]
        }
        if("genus" %in% tax_class[[1]]$rank){
          googlesheet[i,]$genus <- tax_class[[1]]$name[which(tax_class[[1]]$rank=="genus")]
        }
      }
      x <- sci2comm(sci=googlesheet[i,]$sci, db='itis',simplify = T)
      googlesheet[i,]$common <- x[[1]][1]
    }
  }
}





library(ggsci)
write.csv(googlesheet, "results/googlesheet_update.csv", quote = F, row.names = F)
plot.tab <- data.frame(table(googlesheet$class), stringsAsFactors = F)
plot.tab <- plot.tab[order(plot.tab$Freq),]
ggplot(plot.tab, aes(x=reorder(Var1,Freq), y=Freq, fill=Var1)) +
  geom_bar(stat="identity", show.legend = F) +
  coord_flip() +
  theme_classic() +
  scale_fill_aaas() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  


align <- read.fasta("results/combine_filtered_muscle.fa", as.string = T)
cdhit_short_tab <- data.frame(name=names(align), seq=unlist(align))
write.csv(cdhit_short_tab, file="results/combine_filtered_muscle_tab.csv", quote = F)




