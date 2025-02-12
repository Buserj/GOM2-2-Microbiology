## R Script to filter taxa contaminants using stringent negative controls
# Data are QIIME2 outputs in *.qza format,
# sourcetracker2 output in *.tsv format, and sample metadata

setwd("/set/wd")

library(qiime2R)
library(phyloseq)
library(decontam)
library(ape)
library(forcats)
library(rstatix)
library(ggplot2)

colors2 <- c("#D4A5A5", "#4F4A4A", "#C9CBA3", "#A1A46D", "#6E7B8B","#F1C6A5" , "#A1D6A1", "#B4A7D6")
high_contrast_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#8E44AD", "#3498DB", "#2ECC71", "#9B59B6", "#F39C12", "#E67E22", "#1ABC9C", "#2C3E50",
  "#D35400", "#16A085", "#F4D03F", "#D5DBDB", "#C0392B", "#2980B9", "#27AE60", "#8E44AD",
  "#FF6347", "#FFD700", "#D3D3D3", "#1F618D", "#7D3C98")


#########################################
# make phyloseq object and run decontam #
#########################################
# "pool" refers to pseudopooling during DADA2
physeq<-qza_to_phyloseq(features="./pseudopool/pooltable.qza",tax="./pseudopool/pooltax.qza",
                        tree= "./pseudopool/r.pooltree.qza",metadata = "metadata.tsv")
# load in a metadata df if separate from phyloseq
metadata = read.csv("GOM2-2Metadata.csv")
physeq = subset_samples(physeq, sample_names(physeq) != 'lane1-s036-index--CTGCGTAG-CGTGAGTG-29CS5D_S36')
physeq = subset_samples(physeq, sample_names(physeq) != 'lane1-s033-index--CTGCGTAG-TAGCGAGT-26X4C_S33')
physeq = subset_samples(physeq, sample_names(physeq) != 'lane1-s017-index--ACGCTACT-CGTGAGTG-3H1CXX_S17')
physeq = subset_samples(physeq, sample_names(physeq) != 'lane1-s045-index--CGAGCGAC-TAGCGAGT-2H5CXX_S45')

# Contam ASV removal using 'decontam'
# find contams based on prevelance
sample_data(physeq)$is.neg <- sample_data(physeq)$sampletype == "TRUE"
contamdf.prev <- isContaminant(physeq, method="prevalence", neg="is.neg")
contamdf.prev <- tibble::rownames_to_column(contamdf.prev, var="taxa")
table(contamdf.prev$contaminant)

# consider the truth of the matter
contamrow <- which(contamdf.prev$contaminant == "TRUE")
vecp <- contamdf.prev$taxa[contamdf.prev$contaminant == "TRUE"]
tax_table <- (as.data.frame(tax_table(physeq)))
tax_table <- tibble::rownames_to_column(tax_table, var="taxa")
prev_contam <- tax_table[tax_table$taxa %in% vecp,]
# see the Family names
sort(unique(prev_contam$Family), decreasing=FALSE)
# plot the contams
#plot_frequency(physeq, taxa_names(physeq)[contamrow], conc="gDNA")

# kleen, selected based on lm fit
ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, physeq)

#######################
# Deep clean all data #
#######################
# taxa to be gone forever bye bye
# common contaminants to almost always remove for GOM2-2 project
gross <- c("Endozoicomonadaceae","Chloroplast","Mitochondria","Staphylococcaceae",
           "Enterobacteriaceae","Staphylococcales", "Enterobacterales",
           "Streptococcaceae","Micrococcaceae","d__Eukaryota","Enterococcaceae","Nonlabens","Acholeplasma",
           "AEGEAN-169_marine_group","Blastopirellula","Candidatus_Riegeria",
           "Cerasicoccus","Coraliomargarita","Coxiella","Diplosphaera","Erysipelothrix",
           "Peptostreptococcales-Tissierellales","NS4_marine_group","NS2b_marine_group","Bdellovibrio",
           "OM60(NOR5)_clade","Myxococcaceae","Prochlorococcus_MIT9313","Rubidimonas","Synechococcus_CC9902",
           "Tenacibaculum","Trichodesmium_IMS101","Vibrio","Moraxellaceae","Phormidiaceae")
tax <- tax_table(ps.noncontam)
taxa_remove <- rownames(tax)[tax[, "Family"] %in% gross | tax[, "Order"] %in% gross]

# take care of SAR
contains_sar_pattern <- function(row) {
  any(grepl("SAR[0-9]+", row))
}
sar_remove <- rownames(tax)[apply(tax, 1, function(row) contains_sar_pattern(row["Family"]))]
allremove <- c(taxa_remove,sar_remove)

# kill shot!
ps.noncontam <- prune_taxa(!rownames(tax) %in% allremove, ps.noncontam)

##############################
# Visualize the contaminants #
##############################
# sourcetracker2 output from qiime
sourcetracker <- read.table("./pseudopool/st2output_noLOO/mix_prop_st2.tsv",sep="\t",
                            header=TRUE,row.names = 1)

# organize the df
st <- sourcetracker %>%
  cor_gather(drop.na=TRUE)
names(st)[1] <- "cont"
names(st)[2] <- "sample.id"

st$sample.id <- gsub("\\.", "-", st$sample.id)
st <- merge(st,metadata,by="sample.id")
st$cont <- gsub("Unknown", "True Sample", st$cont)
st <- st[st$sample.id != "Undetermined_S0", ]
relabel <- c("1H4C" = "4")

# visualize
st2 <- ggplot(st, aes(x=fct_reorder(name, desc(depthoncenter)),y=cor,fill=cont)) +
  geom_bar(position="fill", stat="identity",color = "black") +
  theme_bw() +
  scale_fill_manual(name = "", values = colors2,
                    labels = c("DNA Extraction", "Library Prep", "Lab Air",
                               "Ship Air","Core Rinds/Pairings","Drilling Mud",
                               "Seawater","True Sediment")) +
  scale_x_discrete(labels=st$depthoncenter) +
  labs(x = "", y = "Taxa Relative Abundance (%)", fill="",title = "") +
  coord_flip()+
  theme_bw(base_size = 18) +
  theme(legend.position  = "bottom", 
        legend.box       = "horizontal",
        legend.direction = "horizontal",
        panel.background = element_blank())
print(st2)
#ggsave("~/Desktop/contamination_source.png", width = 4, height = 4, dpi = 200)

##############################
# Shallow Subsurface Samples # HIGH
##############################
pshigh <-subset_samples(ps.noncontam,depthoncenter < 300)

# extract dfs for more stuff
otu_tablehigh <- (otu_table(pshigh))
# rotate otu matrix layout
rownames(otu_tablehigh) <- factor(rownames(otu_tablehigh), levels = rownames(otu_tablehigh))
otu_mathigh <- as.data.frame(t(otu_tablehigh[rowSums(otu_tablehigh)!=0, ]))
# remove any rows or columns with only 0s
otu_mathigh <- otu_mathigh[, colSums(otu_mathigh !=0)>0]
otu_mathigh <- otu_mathigh[rowSums(otu_mathigh[])>0,]
# match metadata and taxa if any rows were removed
metadatahigh <- metadata[metadata$sample.id %in% rownames(otu_mathigh), ]
tax_tablehigh <- tax_table[tax_table$taxa %in% colnames(otu_mathigh), ]

#tree stuff
tree_name <- taxa_names(phy_tree(pshigh))
common_taxa <- intersect(tree_name, colnames(otu_mathigh))
treeh <- phy_tree(pshigh)
tree_high <- prune_taxa(common_taxa,treeh)

###########################
# Deep Subsurface Samples # LOW
###########################
pslow <-subset_samples(ps.noncontam,depthoncenter > 300)

# extract dfs for more stuff
otu_tablelow <- (otu_table(pslow))
# rotate otu matrix layout
rownames(otu_tablelow) <- factor(rownames(otu_tablelow), levels = rownames(otu_tablelow))
otu_matlow <- as.data.frame(t(otu_tablelow[rowSums(otu_tablelow)!=0, ]))
# remove any rows or columns with only 0s
otu_matlow <- otu_matlow[, colSums(otu_matlow !=0)>0]
otu_matlow <- otu_matlow[rowSums(otu_matlow[])>0,]
# match metadata and taxa if any rows were removed
metadatalow <- metadata[metadata$sample.id %in% rownames(otu_matlow), ]
tax_tablelow <- tax_table[tax_table$taxa %in% colnames(otu_matlow), ]

#tree stuff
tree_name2 <- taxa_names(phy_tree(pslow))
common_taxa2 <- intersect(tree_name2, colnames(otu_matlow))
treel <- phy_tree(pslow)
tree_low <- prune_taxa(common_taxa2,treel)

###########################################
# Save fully clean asv, tax, and metadata #
###########################################
write.csv(otu_matlow, 'otu_deep.csv')
write.csv(tax_tablelow, 'tax_deep.csv')
write.csv(metadatalow, 'metadata_deep.csv')
write.tree(tree_low, 'tree_deep.nwk')

write.csv(otu_mathigh, 'otu_shallow.csv')
write.csv(tax_tablehigh, 'tax_shallow.csv')
write.csv(metadatahigh, 'metadata_shallow.csv')
write.tree(tree_high, 'tree_shallow.nwk')

save(otu_mathigh,otu_matlow,tax_tablehigh,tax_tablelow,metadatalow,metadatahigh,tree_low,tree_high,file="clean_GOM.RData")
