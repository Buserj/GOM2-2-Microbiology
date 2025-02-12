## Calculate Nearest Taxon Index, Net Relatedness Index, Faith's Phylogenetic Diversity,
## beta-Nearest Taxon Index, and Raup-Crick Bray Abundance
# bNTI code sourced from Stegen_etal_ISME_2013

# Input data are a species x sample matrix and phylogenetic tree

library(picante)
library(vegan)
library(ape)
library(iCAMP)
library(dplyr)
library(rstatix)

setwd("your/wd")

# can be set to "shallow" or "deep"
# to reflect data groups
depth <- "deep"
beta.reps = 999

##################################
# load data and square away tree #
##################################
tree <- read.tree(paste0("tree_",depth,".nwk"))
otu <- read.csv(paste0("otu_",depth,".csv"),row.names = 1)
otu <- as.data.frame(t(otu))

# Check needs to be TRUE
identical(rownames(otu), tree$tip.label)
rownames(otu) <- tree$tip.label
match.phylo.data <- picante::match.phylo.data(tree, otu)

#########################################
# phylogenetic analyses: bMNTD and bNTI #
#########################################
# calculate emperical (actual) bMNTD
beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.data$data), 
                                         cophenetic(match.phylo.data$phy),
                                         abundance.weighted = TRUE))
write.csv(beta.mntd.weighted, paste0("weighted_bMNTD_",depth,".csv"), quote = FALSE)

# both are checks: both should be TRUE
identical(colnames(match.phylo.data$data),colnames(beta.mntd.weighted))
identical(colnames(match.phylo.data$data),rownames(beta.mntd.weighted))

# calculate randomized bMNTD
rand.weighted.bMNTD.comp = array(c(-999), dim=c(ncol(match.phylo.data$data),
                                                ncol(match.phylo.data$data),beta.reps))
dim(rand.weighted.bMNTD.comp)

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = 
    as.matrix(comdistnt(t(match.phylo.data$data),
                        taxaShuffle(cophenetic(match.phylo.data$phy)),
                        abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),
                       nrow=ncol(match.phylo.data$data),ncol=ncol(match.phylo.data$data))
dim(weighted.bNTI)

for (columns in 1:(ncol(match.phylo.data$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.data$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
}
rownames(weighted.bNTI) = colnames(match.phylo.data$data)
colnames(weighted.bNTI) = colnames(match.phylo.data$data)
write.csv(weighted.bNTI,paste0("weighted_bNTI_",depth,".csv"),quote=F)

###############################################
# phylogenetic analyses: Raup-Crick Abundance #
###############################################
rc.index <- RC.pc(otu, rand = 999)
rc.bray <- as.data.frame(rc.index$index)

#convert the pairwise matrix into df with cols
rc.bray.df <- rc.bray %>%
  cor_gather(drop.na=TRUE)

write.csv(rc.bray.df,paste0("RCbray_",depth,".csv"),quote=F)

##########################################
# phylogenetic analyses: alpha distances #
##########################################
# Faiths #
otu.pd <- pd(t(otu),tree,include.root=FALSE)

write.csv(otu.pd,paste0("faith_",depth,".csv"))

# Net relatedness index #
cophenDist <- cophenetic.phylo(tree)
otu.ses.pd <- ses.mpd(t(otu), cophenDist)
NRI <- as.matrix(-1 * ((otu.ses.pd$mpd.obs - otu.ses.pd$mpd.rand.mean) / otu.ses.pd$mpd.rand.sd))
rownames(NRI) <- rownames(otu.ses.pd)
colnames(NRI) <- "NRI"

write.csv(NRI,paste0("NRI_",depth,".csv"),quote=F)

# Nearest taxon index #
otu.ses.mntd <- ses.mntd(t(otu), cophenDist)
NTI <- as.matrix(-1 * ((otu.ses.mntd[,2] - otu.ses.mntd[,3]) / otu.ses.mntd[,4]))
rownames(NTI) <- rownames(otu.ses.mntd)
colnames(NTI) <- "NTI"

write.csv(NTI,paste0("NTI_",depth,".csv"),quote=F)
