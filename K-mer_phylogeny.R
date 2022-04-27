library(seqinr)
library(Biostrings)
library(ape)
library(textmineR)
library(ggtree)
library(ggplot2)
library(tibble)
library(phytools)
library(ggstar)
#refernce method (https://bioinformaticshome.com/bioinformatics_tutorials/R/phylogeny_estimation.html, http://yulab-smu.top/treedata-book/chapter13.html)
#Path to seprated fasta files directory
setwd("usr/dir/")

#In the folder containing only the viral sequences
nfiles <- length(dir())
seqdat <- vector("list", nfiles)

for(i in 1:nfiles){
seqdat[[i]] <- read.fasta(file=dir()[i])
}

label <- sapply(1:nfiles, function(k) unlist(strsplit(dir()[k], "[@]"))[1])
names(seqdat) <- label

seqdat_join <- lapply(seqdat, function(k) paste(toupper(unlist(k)),collapse="") )


#
#enumerate all 5-mers
dna <- c("A","C","G","T")
kmer5 <- expand.grid(dna, dna, dna, dna, dna)
kmer5 <- apply(kmer5, 1, function(k) paste(k, collapse=""))

#function for counting all possible kmers (k=5) given a single dna string
kmercount <- function(data){
  sapply(1:length(kmer5), function(k)
  length(unlist(gregexpr2(kmer5[k], data)))
  )}

#vector of counts for 
#all possible kmers (k=5) for all viral sequences
kmer_features <- lapply(seqdat_join, function(k) kmercount(k))


#Collect k-mer counts into a data frame
M <- do.call(rbind, kmer_features)


#The correct input for CalcJSDivergence is the (unnormalised) count vector
JSdist <- CalcJSDivergence(M)

write.table(JSdist, file = "../JSdist-ML-Matrix.tab", row.names = FALSE, dec = ".", sep = "\t", quote = FALSE)

#plot.phylo(bionj(JSdist), type="unrooted", cex=0.8, rotate.tree=95)

#save tree file
as.phylo(bionj(JSdist))->tr
gsub(".fasta", "", tr$tip.label)->tr$tip.label
write.csv(ttree,"../Nodes.csv", row.names = TRUE)
#add tree info 
#add tree info
tp <- read.delim("../treeinfo.csv", sep = ",")

#plot final
ggtree(tr, branch.length='none', layout='circular') %<+% tp +  geom_tiplab(aes(label=genera, color=type, angle=angle, hjust=1, offset=-10, size=0.2), check_overlap = TRUE, size=3) +  xlim(-150, NA) + geom_tippoint(aes(shape=type, color=type), alpha=0.25)

ggsave("../finall-domain-phylogeny-curated.pdf", width = 12, height = 12)
