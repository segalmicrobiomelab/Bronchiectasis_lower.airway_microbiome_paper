# Bronchiectasis Cohort 16S analysis 

#packages needed for the script

library(microViz)
library(microbiome)
library(stringr)
library(qiime2R)
library(metagMisc)
library(YesSiR)
library(patchwork)
library(phyloseq)
library(vegan)
library(pheatmap)
library(ggplot2)
library(ade4)
library(ggforce)
library("matrixStats")
library(RColorBrewer)
library(data.table)
library(microbiome)
library(knitr)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(decontam)
library(DESeq2)
library(microbiomeMarker)
library(stringr)
library(ggprism)
library(forcats)
library(rstatix)
library(patchwork)
library(rmarkdown)
library(ggtext)
library(glue)
library(tidyverse)
#-------------------------------------------------------------------------------------------------------
#set your working directory 
setwd("")

# Load R data object for this analysis
load(file="NTM.analysis.Rdata")

# create phyloseq object
physeq<-qza_to_phyloseq(
  features="no-miss-table-dada2.ntm.qza",
  tree="rooted-tree_quality.ntm.qza",
  taxonomy="taxonomy.ntm.qza",
  metadata = "NTM_Master.file.with_NO_MRN.txt")

#define nromalization for relative table creation 
normalizeSample = function(x) {x/sum(x)}

#save as NTM OTU table
NTM.OTU.Table <- physeq

#subset samples without moc, blank. etc. 
Subset.OTU.Table = subset_samples(NTM.OTU.Table, Description  %in% c("BAL", "BKG", "Sup"))

###KEEP ONLY patients which are for final analysis
Subset.OTU.Table <- subset_samples(Subset.OTU.Table, Cohort_yn == "1")

#keep only baseline BAL samples, this table will be used for BAL only analysis 
Subset.OTU.Table.baseline <- subset_samples(Subset.OTU.Table, Baselines =="1")

#get bal baseline table 
BAL.baseline.table <- subset_samples(Subset.OTU.Table.baseline, Description=="BAL")
#transform for relative tables 
BAL.baseline.table.rel <- transform_sample_counts(BAL.baseline.table, normalizeSample)
NTM.OTU.Table.rel = transformSampleCounts(NTM.OTU.Table, normalizeSample)
Subset.OTU.Table.baseline.rel <- transform_sample_counts(Subset.OTU.Table.baseline, normalizeSample)


##################################  Topographical Analysis of all samples Supp figure 2 ######################

#use Subset.OTU.Table
#relative table
Subset.OTU.Table.rel <- transformSampleCounts(Subset.OTU.Table, normalizeSample)

####ddPCR (supp. figure 2A)

Subset.OTU.Table.metadata <- data.frame(sample_data(Subset.OTU.Table))
Subset.OTU.Table.metadata$Description <- factor(Subset.OTU.Table.metadata$Description, levels = c("BKG", "Sup", "BAL"))
Subset.OTU.Table.metadata$ddPCR_value <- as.numeric(Subset.OTU.Table.metadata$ddPCR_value)

pdf(file = "Figures/All_samples_16S_dPCR.pdf", height = 10, width = 7)
ggplot(data = Subset.OTU.Table.metadata, aes(x=Description, y=ddPCR_value, fill=Description))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "purple", "dodgerblue"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  labs(x= '', y='16S rRNA Gene Copies/uL')+
  scale_x_discrete(labels = c("BKG", "Sup", "BAL"))+
  scale_y_continuous(trans = "log10")+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("BKG", "Sup"), c("Sup", "BAL"),c("BKG", "BAL")), 
                     method = "wilcox.test", label = "p.value", size = 10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))
dev.off()

#############Sequence Detph (Supp figure 2B)
sample_data(Subset.OTU.Table)$Description

# Calculating sequence depth and plot - Read counts

#remove all taxa with 0 abundance
otu.table.sub = subset_taxa(Subset.OTU.Table, rowSums(otu_table(Subset.OTU.Table)) != 0) 

sample_sums(otu.table.sub)
sort(sample_sums(otu.table.sub))

sample_data(otu.table.sub)$Description

sequence_depth.sub <- sample_sums(otu.table.sub)
Sample_Type.sub <- sample_data(otu.table.sub)$Description

sequence.per.sample.sub <- data.frame(Sample_Type.sub, sequence_depth.sub)
sequence.per.sample.sub$Sample_Type.sub <- factor(sequence.per.sample.sub$Sample_Type.sub, levels = c("BKG", "Sup", "BAL"))

#get stats 
all_seq_data <- sequence.per.sample.sub %>% 
  summarise(mean_seq=mean(sequence_depth.sub), 
            median_seq=median(sequence_depth.sub), 
            q1_seq=quantile(sequence_depth.sub, 0.25), 
            q3_seq=quantile(sequence_depth.sub, 0.75))
#save 
write.csv(all_seq_data, file = "Results/sequence_depth_data.csv")

seq_table_sample_type <- sequence.per.sample.sub %>% 
  group_by(Sample_Type.sub) %>% 
  summarise(mean_seq=mean(sequence_depth.sub), 
            median_seq=median(sequence_depth.sub), 
            q1_seq=quantile(sequence_depth.sub, 0.25), 
            q3_seq=quantile(sequence_depth.sub, 0.75))

#save it 
write.csv(seq_table, file = "Results/sequence_depth_data_sample_type.csv")
#stats 
my.comparisons <- compare_means(sequence_depth.sub~Sample_Type.sub, data = sequence.per.sample.sub )

# Colored boxplot according to sample type
pdf(file = "Figures/Sequence.depth.All.samples.pdf"
    , height = 10, width = 7)
ggplot(data=sequence.per.sample.sub, mapping = aes(x=Sample_Type.sub, y=sequence_depth.sub, fill=Sample_Type.sub))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "purple", "dodgerblue"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  labs(x= '', y='Total Reads')+
  scale_x_discrete(labels = c("BKG", "Sup", "BAL"))+
  scale_y_continuous(trans = "log10")+
  guides(fill="none")+
  stat_compare_means(comparisons = list(c("BKG", "Sup"), c("Sup", "BAL"),c("BKG", "BAL")), 
                     method = "wilcox.test", label = "p.value", size = 10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))
dev.off()


################# Alpha Diversity (Supp figure 2C)

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(Subset.OTU.Table, measures = "Shannon"), 
                             "Observed" = estimate_richness(Subset.OTU.Table, measures = "Observed"),
                             "Simpson" = estimate_richness(Subset.OTU.Table, measures = "Simpson"),
                             "Sample_Type" = sample_data(Subset.OTU.Table)$Description)
alpha.measures <- alpha.measures %>% 
  mutate(Sample_Type = factor(Sample_Type, levels = c("BKG", "Sup", "BAL")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~Sample_Type, data = alpha.measures)

# plot it
pdf(file = "Figures/Shannon.All.samples.baseline.and.long.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=Sample_Type, y=Shannon, fill=Sample_Type))+
  geom_boxplot()+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkgrey", "purple", "dodgerblue"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("BKG", "Sup", "BAL"))+
  stat_compare_means(comparisons = list(c("BKG", "Sup"), c("Sup", "BAL"),c("BKG", "BAL")), 
                     method = "wilcox.test", label = "p.value", size = 10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

########Beta Diversity( Supp figure 2D)

#Create Distance Matrix
vegdist = vegdist(t(otu_table(Subset.OTU.Table.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(Subset.OTU.Table.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Description,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Description",suffixes=c("",".centroid"))

newResults$Description <- factor(newResults$Description, levels = c("BKG", "Sup", "BAL"))
centroids$Description <- factor(centroids$Description, levels = c("BKG", "Sup", "BAL"))

#sats
adonis2(vegdist ~ sample_data(Subset.OTU.Table.rel)$Description)

pdf(file = "Figures/Beta.Diversity.Bray.All.Samples.baseline.and.long.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Description)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("darkgrey", "purple", "dodgerblue")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Description), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Description)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("BAL", "BKG", "Sup")), size=10) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.001")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


###############Decontamination Analysis (supp figure 3) ####################
####Functions are defined in separate script 

decontaminant_subplot_KW(input_phyloseq = Subset.OTU.Table.baseline, 
                         sample_type_var_name = "Description", 
                         sample_types = c("BKG", "BAL", "Sup"), 
                         sample_type_color = c("darkgrey","purple","dodgerblue"), 
                         sample_type_color_2nd = c("darkgrey","purple","dodgerblue"), 
                         negative_sample_type = "BKG", 
                         compare_type = c("BAL_and_Sup"),
                         stat_option = "mean", 
                         method_type = "prev",
                         test_threshold = 0.25, 
                         graph_option=c("boxplot"))
                         

#run it so it can produce mean and SE bars and not IQR 
decontaminant_subplot_KW(input_phyloseq = Subset.OTU.Table.baseline, 
                         sample_type_var_name = "Description", 
                         sample_types = c("BKG", "BAL", "Sup"), 
                         sample_type_color = c("darkgrey","purple","dodgerblue"), 
                         sample_type_color_2nd = c("darkgrey","purple","dodgerblue"), 
                         negative_sample_type = "BKG", 
                         compare_type = c("BAL_and_Sup"),
                         stat_option = "mean", 
                         method_type = "prev",
                         test_threshold = 0.25, 
                         graph_option=c("mean_SE"))

#Final Contamlist is the file you get as an output from decontam analysis 

contamlist <- read.csv("contaminant_preval_negC_BKG_compare_BAL_and_Supthres_0.25_OP.csv")
contamlist<- contamlist %>% 
  dplyr::select(Row.names, contaminant) %>% 
  dplyr::mutate(asv= str_sub(Row.names, start = -32L, end = -1L))

#save your contamlist 
write.csv(contamlist, file = "Results/final_contamlist.csv")

##############################Differential Analysis BAL vs upper airway (Supp Figure 4)####

#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 

#define function
phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if (identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts = x, group = group, genes = taxonomy, remove.zeros = TRUE, 
              ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method = method)
  # Check for division by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop("Something wrong with edgeR::calcNormFactors on this data,\n non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}


#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = Subset.OTU.Table.baseline, group = "Description", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
res$logFC <- res$logFC*(-1)

# make taxa name column and name asv column 
res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
NTM.OTU.no.BKG.rel.BAL <- subset_samples(Subset.OTU.Table.baseline.rel, Description=="BAL")
NTM.OTU.no.BKG.rel.Sup <- subset_samples(Subset.OTU.Table.baseline.rel, Description=="Sup")

RA.df.1 <- data.frame(otu_table(NTM.OTU.no.BKG.rel.pruned.BAL))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.BAL <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(NTM.OTU.no.BKG.rel.pruned.Sup))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Sup <- meanRA.save

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 

#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "BAL", "Sup"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.BAL, abundance.Sup, category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results.BAL_vs_Sup_Pruned.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="BAL", res$abundance.BAL, res$abundance.Sup)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

# Create variable for color gradient of the dots in the Volcano Plot
cols <- densCols(res$logFC, res$sig)
cols[res$PValue ==0] <- "purple"
cols[res$logFC > 1 & res$FDR < 0.2 ] <- "dodgerblue"
cols[res$logFC < -1 & res$FDR < 0.2 ] <- "purple"

#Bubble plot of top 25 taxa from each group 
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"dodgerblue", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "purple","grey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>%
  dplyr::mutate(color=ifelse(res.bubble.filtered$contaminant=="FALSE", "black", "red")) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("dodgerblue", "purple"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_BAL_vs_Sup_Pruned_with_contam_labels.pdf", height = 11, width = 12)
bubble.plot+
  guides(size="none")
dev.off()

#get separate legend of relative abundance 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_BAL_vs_Sup_Pruned_with_contam_labels_Legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 



#########################Prune the Data and continue analysis for lower airway samples only ########################
#subset a table of Sup and BAL only 
NTM.no.BKG <- subset_samples(Subset.OTU.Table, Description %in% c("BAL", "Sup"))

#make sure you have only baseline and cohor 
NTM.no.BKG.baseline <- subset_samples(NTM.no.BKG, Cohort_yn == "1")

#make sure you get only baselines 
NTM.no.BKG.baseline <- subset_samples(NTM.no.BKG, Baselines == "1")

#create relative table 
NTM.no.BKG.baseline.rel <- transform_sample_counts(NTM.no.BKG.baseline, normalizeSample)

#prune the data and use them for all analysis that combine both types of samples later 
pruned.taxa <- genefilter_sample(NTM.no.BKG.baseline.rel, filterfun_sample(function(x) x> 0.002), 
                                 A = 0.003 * nsamples(NTM.OTU.no.BKG.rel))

#prune relative table 
NTM.no.BKG.baseline.rel.pruned <- prune_taxa(pruned.taxa, NTM.no.BKG.baseline.rel)

#see how much taxa you have 
NTM.no.BKG.baseline.rel.pruned

#prune count table 
NTM.no.BKG.baseline.pruned <- prune_taxa(pruned.taxa, NTM.no.BKG.baseline)

#get BAL table only 
BAL.baseline.table.pruned <- subset_samples(NTM.no.BKG.baseline.pruned, Description=="BAL")

#get relative table 
BAL.baseline.table.pruned.rel <- transform_sample_counts(BAL.baseline.table.pruned, normalizeSample)



#################################################################################
############################ Comparison BAL  NTM pos vs Neg (Figure 1A-D) ################
#############################################################################

#choose baseline BAL only (pruned already)
BAL.baseline.table.pruned@sam_data$Primary.lobe
BAL.baseline.table.pruned.rel
#choose most involved samples (should get 200 samples)
BAL.baseline.table.pruned.prim <- subset_samples(BAL.baseline.table.pruned, Primary.lobe!="N.A." )

#obtain relative table 
BAL.baseline.table.pruned.prim.rel <- transform_sample_counts(BAL.baseline.table.pruned.prim, normalizeSample)

#remove samples without data on NTM pos or neg (if any)
BAL.baseline.table.pruned.prim <- subset_samples(BAL.baseline.table.pruned.prim, NTM_pos_neg != "N.A.")
BAL.baseline.table.pruned.prim.rel <- subset_samples(BAL.baseline.table.pruned.prim.rel, NTM_pos_neg != "N.A.")

#variable to compare : NTM pos vs neg 
BAL.baseline.table.pruned.prim@sam_data$NTM_pos_neg

#########Alpha diversity

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(BAL.baseline.table.pruned.prim, measures = "Shannon"), 
                             "Observed" = estimate_richness(BAL.baseline.table.pruned.prim, measures = "Observed"), 
                             "Simpson" = estimate_richness(BAL.baseline.table.pruned.prim, measures = "Simpson"),
                             "NTM_status" = sample_data(BAL.baseline.table.pruned.prim)$NTM_pos_neg)
alpha.measures <- alpha.measures %>% 
  mutate(NTM_status = factor(NTM_status, levels = c("0", "1")))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~NTM_status, data = alpha.measures)

# plot it
pdf(file = "Figures/Shannon.BAL.samples.NTM.pos.vs.neg.baseline_Inv_max.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=NTM_status, y=Shannon, fill=NTM_status))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("skyblue", "orangered"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Negative", "Positive"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  guides(fill="none")+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#####Beta Divesity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(BAL.baseline.table.pruned.prim.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(BAL.baseline.table.pruned.prim.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ NTM_pos_neg,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="NTM_pos_neg",suffixes=c("",".centroid"))

newResults$NTM_pos_neg <- factor(newResults$NTM_pos_neg, levels = c("0", "1"))
centroids$NTM_pos_neg <- factor(centroids$NTM_pos_neg, levels = c("0", "1"))

#stats 
adonis2(vegdist ~ sample_data(BAL.baseline.table.pruned.prim.rel)$NTM_pos_neg)
#plot
pdf(file = "Figures/Beta.Diversity.Bray.BAL.samples.NTM.pos.vs.neg.baseline_Inv_max.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= NTM_pos_neg)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("skyblue", "orangered")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= NTM_pos_neg), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= NTM_pos_neg)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Negative", "Positive")), size=10) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.418")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

####ddPCR figure 
pdf(file = "Figures/ddPCR.BAL.samples.NTM.pos.vs.neg.baseline_Inv_max.pdf", width=7, height=10)
ggplot(newResults, aes(x=NTM_pos_neg, y=as.numeric(ddPCR_value), fill=NTM_pos_neg))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("skyblue", "orangered"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("16S_rRNA Copies/uL")+
  scale_x_discrete(labels = c("Negative", "Positive"))+
  scale_y_continuous(trans = "log10")+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  guides(fill="none")+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()


#############################Differential Analysis
ps <- BAL.baseline.table.pruned.prim
ps.rel <- BAL.baseline.table.pruned.prim.rel
    
#prepare data by converting phyloseq to edgeR and performing analysis.
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = BAL.baseline.table.pruned.prim, group = "NTM_pos_neg", method = "TMM")
    
#run test 
ET <- exactTest(dge)
    
#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)
    
# make taxa name column and name asv column 
    
res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                              ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                     ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                            ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                                   ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                          ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                                 ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))
    
#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                                ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                       ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                              ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                     ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                            ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                                   ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))
    
#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
      dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
    res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))
    
#get relative abundance data 
ps.rel.1 <- subset_samples(BAL.baseline.table.pruned.prim.rel, NTM_pos_neg=="1")
ps.rel.2 <- subset_samples(BAL.baseline.table.pruned.prim.rel, NTM_pos_neg=="0")
    
RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.NTM.pos <- meanRA.save
    
#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.NTM.neg <- meanRA.save
    
#clean the results df 
    
#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]
    
#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))
    
# drop any NA values 
res<- res %>% 
      drop_na(., FDR) %>% 
      dplyr::arrange(desc(FDR))
    
#get abundance of each category alone to save it in the table 
    
#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "0"))
    
res <- inner_join(res, contamlist, by="asv")
    
#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.NTM.pos, abundance.NTM.neg, category, contaminant) %>% 
  arrange(., desc(category))
    
# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_NTM_pos_vs_neg_BAL_Pruned.csv")
    
    
#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.NTM.pos, res$abundance.NTM.neg)
    
#Bubble plot#
res.bubble <- res
        
#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)
        
#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)
#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)
        
#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"skyblue", 
                                          ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orangered","grey"))
        
#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                                   paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))
        
res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))
        
#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")
        
bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
          geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
                       aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
          scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
          scale_fill_manual(values = c("skyblue", "orangered"))+
          theme(panel.background = element_blank(),
                panel.border=element_rect(fill=NA),
                panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
                panel.grid.minor = element_blank(),
                strip.background=element_blank(),
                axis.title=element_text(size=20,face="bold"),
                axis.text.x=element_text(colour="black", size=18, face="bold"),
                axis.text.y=element_markdown(face="bold",size=16),
                axis.ticks=element_line(colour="black"),
                legend.background = element_rect(color=NA), 
                plot.title=element_text(size=23, face="bold"))+
          xlab("") +
          ylab("")+
          ggtitle("EdgeR")+
          geom_vline(xintercept=0, color="red",linetype="dashed")+
          guides(fill="none")
        
        pdf(file="Figures/edgeR_Bubble_NTM_pos_vs_neg_BAL_Pruned_with_contam_labels.pdf", height = 11, width = 8)
        bubble.plot+
          guides(size="none")
        dev.off()
        
        #get seperate legend plot 
        my_legend <- get_legend(bubble.plot)
        #save it 
        pdf(file="Figures/edgeR_Bubble_NTM_pos_vs_neg_BAL_Pruned_with_contam_labels_Legend.pdf", 
            height = 4, width = 2)
        as_ggplot(my_legend)
        dev.off() 
        


####### DMM clustering analysis ###### 


library(DirichletMultinomial)
library(lattice)
library(parallel)

# table will be used. without BKG

BAL.baseline.table.pruned.prim

#starting DMM
#convert table to a matrix of counts 
dat <- otu_table(BAL.baseline.table.pruned.prim)
count <- t(dat)
nrow(count)
colnames(count) <- count[1,]

#change class to numeric
class(count) <- "numeric"

#calculate DMM 
fit.all <- lapply(1:7, dmn, count=count, verbose=TRUE)

#check model fit based on laplace
lplc <- base::sapply(fit.all, DirichletMultinomial::laplace)

#plot to find best model 
num.of.components <- c(1:7)
plot.df <- data.frame(lplc, num.of.components)

#plot laplace figure (Supp figure 5)
pdf(file="Figures/DMM.fit.model.pdf", height = 6, width = 6)
ggplot(plot.df, aes(x=num.of.components, y=lplc))+
  geom_point()+
  geom_line()+
  labs(x="Number of Drichlet Components", y="Model Fit (Laplace)")+
  theme_bw()
dev.off()

#find best model 
best <- fit.all[[which.min(unlist(lplc))]]

#extract best model to use later 
dmm_3 <- fit.all[[3]]

#sample component assigments 
ass <- apply(mixture(best), 1, which.max)

#starting a data frame of groups assignments
ass.df <- data.frame(ass)

colnames(ass.df)[1] <- "dmm_cluster"

#### add dmm cluseter to your table 
BAL.baseline.table.pruned.prim@sam_data$dmm_cluster <- ass.df$dmm_cluster


############### DMM clusters analysis in all samples (Supp figure 6)####


ps <- BAL.baseline.table.pruned.prim
ps.rel <- transform_sample_counts(ps, normalizeSample)    

#ddpcr differences 

dp.dat <- data.frame(sample_data(ps))
dp.dat$dmm_cluster <- as.factor(dp.dat$dmm_cluster)
dp.dat$ddPCR_value <- as.numeric(dp.dat$ddPCR_value)
#calculate stats for ddPCR 
my.comparisons <- compare_means(ddPCR_value~dmm_cluster, data = dp.dat)
#save it
write.csv(my.comparisons, file = "Results/ddPCR.stats.all_samples_dmm_cluster.csv")

# plot it
pdf(file = "Figures/ddPCR.all_samples_.DMM.pdf", width=7, height=10)
ggplot(dp.dat, aes(x=dmm_cluster, y=ddPCR_value, color=dmm_cluster))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_color_manual(values = c("orange", "red2", "blue" ))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("16S rRNA Gene Copies")+
  scale_y_log10()+
  scale_x_discrete(labels = c("1", "2", "3"))+
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")
dev.off()


#########Alpha diversity##### 

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "dmm_cluster" = sample_data(ps)$dmm_cluster)
alpha.measures <- alpha.measures %>% 
  mutate(dmm_cluster= factor(dmm_cluster))


#calculate stats for shannon
my.comparisons <- compare_means(Shannon~dmm_cluster, data = alpha.measures)
#save it
write.csv(my.comparisons, file = "Results/Shannon.stats.all_samples.dmm_cluster.csv")

# plot it
pdf(file = "Figures/Shannon.all_samples.dmm_cluster.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=dmm_cluster, y=Shannon, color=dmm_cluster))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_color_manual(values = c("orange", "red2", "blue" ))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("1", "2", "3"))+
  stat_compare_means(comparisons = list(c("1", "2"), c("2", "3"), c("1", "3")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(color="none")
dev.off()



###Beta Diversity####

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ dmm_cluster,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="dmm_cluster",suffixes=c("",".centroid"))

newResults$dmm_cluster <- factor(newResults$dmm_cluster, levels = c("1", "2", "3"))
centroids$dmm_cluster <- factor(centroids$dmm_cluster, levels = c("1", "2", "3"))

pdf(file = "Figures/Beta.Diversity.Bray.BAL.samples.dmm_cluster_all_samples.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= dmm_cluster)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("orange", "red2", "blue" )) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= dmm_cluster), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= dmm_cluster)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("1", "2", "3"), size=10)) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.005")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "darkgrey"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "darkgrey", size = rel(0.75)),axis.text.y=element_text(colour = "darkgrey", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

#sats
write.table(adonis2(vegdist ~ newResults$dmm_cluster), 
            file = "Results/Beta.Diversity.Bray.BAL.samples.dmm_cluster_all_samples.txt", sep = "\t", row.names = T)


#1 vs 2
ps.re.1 <- subset_samples(ps.rel, dmm_cluster %in% c("1", "2"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$dmm_cluster), file = "Bray.Stats.all_samples_dmm_cluster_1_vs_2.csv")

#2 vs 3 
ps.re.1 <- subset_samples(ps.rel, dmm_cluster %in% c("2", "3"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$dmm_cluster), file = "Bray.Stats.all_samples_dmm_cluster_2_vs_3.csv")

#1 vs 3
ps.re.1 <- subset_samples(ps.rel, dmm_cluster %in% c("1", "3"))
vegdist   = vegdist(t(otu_table(ps.re.1)), method = "bray")
write.csv(adonis2(vegdist~sample_data(ps.re.1)$dmm_cluster), file = "Bray.Stats.all_samples_dmm_cluster_1_vs_3.csv")


############################Differential - only edgeR #####################################
#cluster 1 vs 2 

#set tables to use 
ps <- subset_samples(BAL.baseline.table.pruned.prim, dmm_cluster != "3")

ps.rel <- transform_sample_counts(ps, normalizeSample)    

ps
ps.rel
################################################################################
#####################edgeR ##########################################

#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "dmm_cluster", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here#############


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)


# make taxa name column and name asv column 

res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, dmm_cluster=="1")
ps.rel.2 <- subset_samples(ps.rel, dmm_cluster=="2")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance_dmm_1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance_dmm_2 <- meanRA.save

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 

#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "2", "1"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance_dmm_1, abundance_dmm_2,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_dmm_cluster_all_samples.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance_dmm_1, 
                        res$abundance_dmm_2)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"red", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange","darkgrey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "orange2", "red2"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_all_samples_cluster_2_vs_1_category_BAL_Pruned_with_contam_labels.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_all_samples_cluster_2_vs_1_category_BAL_Pruned_with_contam_labels_legeneds.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 


##### cluster  3 vs 2####

#set tables to use 
ps <- subset_samples(BAL.baseline.table.pruned.prim, dmm_cluster != "1")

ps.rel <- transform_sample_counts(ps, normalizeSample)    

ps
ps.rel

#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "dmm_cluster", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here#############


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)


# make taxa name column and name asv column 

res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, dmm_cluster=="2")
ps.rel.2 <- subset_samples(ps.rel, dmm_cluster=="3")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance_dmm_2 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance_dmm_3 <- meanRA.save

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 

#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "3", "2"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance_dmm_2, abundance_dmm_3,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_dmm_cluster_3_vs_2_all_samples.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="2", res$abundance_dmm_2, 
                        res$abundance_dmm_3)

###### volcano plot####
# create color variable 
res$sig <- -log10(res$FDR)


#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"blue", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "red","darkgrey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "blue", "red"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_all_samples_cluster_3_vs_2_category_BAL_Pruned_with_contam_labels.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_all_samples_cluster_3_vs_2_category_BAL_Pruned_with_contam_labels_legened.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 



#########  
##### cluster  1 vs 3####

#set tables to use 
ps <- subset_samples(BAL.baseline.table.pruned.prim, dmm_cluster  != "2")

ps.rel <- transform_sample_counts(ps, normalizeSample)    

ps
ps.rel

#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "dmm_cluster", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

############continue edgeR here#############


res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)


# make taxa name column and name asv column 

res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, dmm_cluster=="1")
ps.rel.2 <- subset_samples(ps.rel, dmm_cluster=="3")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance_dmm_1 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance_dmm_3 <- meanRA.save

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 

#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "3", "1"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance_dmm_1, abundance_dmm_3,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_dmm_cluster_1vs3_all_samples.csv")


#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance_dmm_1, 
                        res$abundance_dmm_3)


#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"blue", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "orange","darkgrey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="grey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "blue", "orange"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_all_Samples_dmm_cluster_1vs3_category_BAL_Pruned_with_contam_labels.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_all_Samples_dmm_cluster_1vs3_category_BAL_Pruned_with_contam_labels_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 







#######Figure 2 Biplot ####
library(microViz)
#plot NTM pos only 

#define table to use
ps <- subset_samples(BAL.baseline.table.pruned.prim, NTM_pos_neg=="1")

#fix taxa
ps_analysis <- tax_fix(ps)
ps_analysis <- phyloseq_validate(ps_analysis, remove_undetected = TRUE)

# first we make a function that replaces any unwanted "_" in our taxa labels with spaces
renamer <- function(x) str_replace(x, pattern = "_", replacement = " ")

#generate the plot and save it 
p <- ps_analysis %>%
  ps_mutate(
    NET=as.numeric(NET_Assay_Primary_lobes_catg=="High"),
    Macrophages = as.numeric(Macrophages_catg == "High_Macro"), 
    Neutrophils=as.numeric(Neutrophils_catg=="High_Neu")
  ) %>%
  tax_fix(unknowns = c("Bacillus", "Clostridium")) %>% 
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("Macrophages", "Neutrophils", "NET"),
    method = "RDA",
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    size = 8, alpha = 0.5, color="orangered",
    auto_caption = NA, # remove the helpful automatic caption
    plot_taxa = 1:7, taxon_renamer = renamer, # renamer is the function we made earlier
    tax_vec_length = 6, # this value is a scalar multiplier for the biplot score vectors
    tax_lab_length = 7, # this multiplier moves the labels, independently of the arrowheads
    tax_lab_style = tax_lab_style(size = 6, alpha = 0.5), # create a list of options to tweak the taxa labels' default style
    constraint_vec_length = 6, # this adjusts the length of the constraint arrows, and the labels track these lengths by default
    constraint_vec_style = vec_constraint(6, alpha = 0.5), # this styles the constraint arrows
    constraint_lab_style = constraint_lab_style(size = 6) # this styles the constraint labels
  ) +
  # the functions below are from ggplot2:
  # You can pick a different colour scale, such as a color_brewer palette
  #scale_color_manual(values = "orangered")+
  # You can set any scale's values manually, such as the shapes used
  #scale_shape_manual(values = c(
  # active = "circle", mild = "circle cross",
  #inactive = "circle open", control = "square open"
  #) +
  # this is how you add a title and subtitle
  #ggtitle(
  # label = "[Insert your exciting interpretations here?]",
#  subtitle = "RDA with clr-transformed genera: constraints in red, taxa in black"
#) +
# and this is how you make your own caption
#labs(caption = "91 samples, 178 genera. Type 2 scaling.") +
# this is one way to set the aspect ratio of the plot
coord_fixed(ratio = 1, clip = "off")+
  #labs(color = "NTM Status")+
  xlab("RDA1")+ylab("RDA2")+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        axis.title=element_text(size=30,face="bold"),
        axis.ticks=element_blank(),
        plot.margin=unit(c(1,1,1,1),"line"))+
  guides(color="none")

#save it 
pdf(file = "Figures/biplot_NET_cell_count_NTM_pos_samples.pdf", height = 15, width = 20)
p
dev.off()


#repeat the above for every cell type or NET for other NTM groups 



#### Figure 3 - NTM relative abundance 

####### Use Phyloseq function psmelt to get abundance data
abundance.data <- psmelt(BAL.baseline.table.prim.rel)

#arrange the df                 
abundance.data <- abundance.data %>% 
  dplyr::rename(., ASV=OTU) %>% 
  mutate(meanRA=mean(Abundance)) %>% 
  dplyr::arrange(ASV)

#get taxa names and name with ASV 
abundance.data$names <- paste(ifelse(!is.na(abundance.data$Species), paste("g_",abundance.data$Genus,"_s_",abundance.data$Species,sep=""),
                                     ifelse(!is.na(abundance.data$Genus), paste("g_",abundance.data$Genus,sep=""),
                                            ifelse(!is.na(abundance.data$Family), paste("f_",abundance.data$Family,"_g__",sep=""),
                                                   ifelse(!is.na(abundance.data$Order), paste("o_",abundance.data$Order, "_f__g__",sep=""),
                                                          ifelse(!is.na(abundance.data$Class), paste("c_",abundance.data$Class, "_o__f__g__",sep=""),
                                                                 ifelse(!is.na(abundance.data$Phylum), paste("p_",abundance.data$Phylum, "_c__o__f__g__",sep=""),
                                                                        ifelse(!is.na(abundance.data$Domain), paste("k_",abundance.data$Domain, "_p__c__o__f__g__",sep=""), paste(abundance.data$ASV)))))))))


#create another column of nameASV with last three digits from ASV 
abundance.data$nameASV.trim <- paste0(abundance.data$names,paste0(".."), 
                                      paste0(str_sub(abundance.data$ASV, - 3, - 1)))

#create a new column differentiating mycobacterium from other taxa 
abundance.data <- abundance.data %>% 
  dplyr::mutate(myco_or_not = ifelse(Genus == "Mycobacterium", nameASV.trim,"Non-Mycobacterium")) %>% 
  dplyr::mutate(myco_or_not=factor(myco_or_not))

#fill na values with non mycobacterium 
abundance.data$myco_or_not[is.na(abundance.data$myco_or_not)] <- "Non-Mycobacterium"

# Calculate the percentage values
abundance.data$percent_abundance <- abundance.data$Abundance * 100

# change levels of myco or not 
abundance.data$myco_or_not <- factor(abundance.data$myco_or_not, levels = rev(levels(abundance.data$myco_or_not)))

#seperate data frames 
abundance.data.myco <- abundance.data %>% 
  dplyr::filter(myco_or_not !="Non-Mycobacterium") %>% 
  arrange(desc(percent_abundance))

abundance.data.non_myco <- abundance.data %>% 
  dplyr::filter(myco_or_not == "Non-Mycobacterium")

#number of patietns with positive mycobacterium sequence - overall 40 patients, 35 out of them in NTM positive group
abundance.data.myco %>% 
  distinct(SubjID, .keep_all = TRUE) %>% 
  dplyr::filter(percent_abundance!=0) %>% 
  filter(NTM_pos_neg=="1")

#now combine the two datafrmaes again 
abundance.data.sorted <- rbind(abundance.data.myco, abundance.data.non_myco)

#crewate seperate figures for NTM pos and NTM neg 
abundance.data.NTM.pos <- abundance.data.sorted %>% 
  dplyr::filter(NTM_pos_neg=="1")

#create order column 
abundance.data.NTM.pos <- abundance.data.NTM.pos %>% tibble::rownames_to_column(var="RowID") 
abundance.data.NTM.pos$Sample <- factor(abundance.data.NTM.pos$Sample, levels = unique(abundance.data.NTM.pos$Sample))

# Create the bar chart with percent scale for NTM pos 
p1 <- ggplot(abundance.data.NTM.pos, aes(x =fct_reorder(Sample, RowID), y= percent_abundance, fill=myco_or_not)) +
  geom_bar(stat = "identity", alpha=0.5)+
  scale_fill_manual(values = c("Non-Mycobacterium"="lightgrey",
                               "g_Mycobacterium..d3a" ="red", 
                               "g_Mycobacterium..c6e" ="blue", 
                               "g_Mycobacterium..a6b" = "skyblue", 
                               "g_Mycobacterium..914" = "cyan3", 
                               "g_Mycobacterium..852" = "orange", 
                               "g_Mycobacterium..69e" = "orangered", 
                               "g_Mycobacterium..647" = "pink1", 
                               "g_Mycobacterium..54a" = "lightgreen",
                               "g_Mycobacterium..356" = "green", 
                               "g_Mycobacterium..1ee" = "brown4", 
                               "g_Mycobacterium..1a8" = "tomato2",
                               "g_Mycobacterium_s_vaccae..687" = "yellow",
                               "g_Mycobacterium_s_llatzerense..8c2" ="yellow4",
                               "g_Mycobacterium_s_gordonae..700" = "purple1",
                               "g_Mycobacterium_s_celatum..f72" = "gold1"
  ))+
  ylab("Relative Abundance")+xlab("")+
  scale_y_continuous(breaks = seq(0,100,by=10))+
  theme_classic()+
  ggtitle("NTM positive")+
  guides(fill="none")+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold", angle = 90),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))

# Create the bar chart with percent scale for NTM neg 
abundance.data.NTM.neg <- abundance.data.sorted %>% 
  dplyr::filter(NTM_pos_neg=="0")

#create order column 
abundance.data.NTM.neg <- abundance.data.NTM.neg %>% tibble::rownames_to_column(var="RowID") 
abundance.data.NTM.neg$Sample <- factor(abundance.data.NTM.neg$Sample, levels = unique(abundance.data.NTM.neg$Sample))

p2 <- ggplot(abundance.data.NTM.neg, aes(x =fct_reorder(Sample, RowID), y= percent_abundance, fill=myco_or_not)) +
  geom_bar(stat = "identity", alpha=0.5) +
  scale_fill_manual(values = c("Non-Mycobacterium"="lightgrey",
                               "g_Mycobacterium..d3a" ="red", 
                               "g_Mycobacterium..c6e" ="blue", 
                               "g_Mycobacterium..a6b" = "skyblue", 
                               "g_Mycobacterium..914" = "cyan3", 
                               "g_Mycobacterium..852" = "orange", 
                               "g_Mycobacterium..69e" = "orangered", 
                               "g_Mycobacterium..647" = "pink1", 
                               "g_Mycobacterium..54a" = "lightgreen",
                               "g_Mycobacterium..356" = "green", 
                               "g_Mycobacterium..1ee" = "brown4", 
                               "g_Mycobacterium..1a8" = "tomato2",
                               "g_Mycobacterium_s_vaccae..687" = "yellow",
                               "g_Mycobacterium_s_llatzerense..8c2" ="yellow4",
                               "g_Mycobacterium_s_gordonae..700" = "purple1",
                               "g_Mycobacterium_s_celatum..f72" = "gold1"
  ))+
  ylab("Relative Abundance")+xlab("")+
  scale_y_continuous(breaks = seq(0,100,by=10))+
  theme_classic()+
  ggtitle("NTM Negative")+
  labs(fill="")+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold", angle = 90),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"), 
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.text = element_text(size=16)) #change legend text font size

#put the two figures together 
pdf(file = "Figures/relative.abundance.of.mycobacterium.pdf", width = 60, height = 20)
ggarrange(p1, p2, ncol = 2, nrow = 1)
dev.off()

#if you want to add zoom in on specefic values 
p1+  facet_zoom(ylim=c(0,10))
p2+ facet_zoom(ylim = c(0,10))

#figure with zoom 
p1_zoomed <- p1+facet_zoom(ylim = c(0,10))
p2_zoomed <- p2+facet_zoom(ylim = c(0,10))+theme(legend.position = "top", 
                                                 legend.key.size = unit(0.5, 'cm'), #change legend key size
                                                 legend.text = element_text(size=12))
p2_zoomed_no_legend <- p2_zoomed+guides(fill="none")

#export 
pdf(file = "Figures/Relateive_abundance_zoomed.pdf", width = 50, height = 40)
ggarrange(p1_zoomed, p2_zoomed_no_legend, ncol = 1, nrow = 2)
dev.off()

#get seperate legend plot 
my_legend <- get_legend(p2_zoomed)
#save it 
pdf(file="Figures/Relateive_abundance_zoomed_legend.pdf", 
    height = 6, width = 10)
as_ggplot(my_legend)
dev.off() 







######### Figure 4A-D###########

cav.NTM.pos.table <- subset_samples(BAL.baseline.table.pruned.prim, NTM_pos_neg=="1")
cav.NTM.pos.table.rel <- transform_sample_counts(cav.NTM.pos.table, normalizeSample)    

ps <- cav.NTM.pos.table
ps.rel <- cav.NTM.pos.table.rel

#ddpcr differences 
dp.dat <- data.frame(sample_data(ps))
dp.dat$Cavitary <- as.factor(dp.dat$Cavitary)
#calculate stats for ddPCR 
my.comparisons <- compare_means(ddPCR_value~Cavitary, data = dp.dat)

# plot it
pdf(file = "Figures/ddPCR.NTM_pos.Cavitary.pdf", width=7, height=10)
ggplot(dp.dat, aes(x=Cavitary, y=ddPCR_value, fill=Cavitary))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("azure4", "aquamarine4"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("16S rRNA Gene Copies")+
  scale_y_log10()+
  scale_x_discrete(labels = c("Non Cavitary", "Cavitary"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#########Alpha diversity
#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "Cavitary" = sample_data(ps)$Cavitary.)
alpha.measures <- alpha.measures %>% 
  mutate(Cavitary= factor(Cavitary))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~Cavitary, data = alpha.measures)

# plot it
pdf(file = "Figures/Shannon.NTM_pos.Cavitary.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=Cavitary, y=Shannon, fill=Cavitary))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("azure4", "aquamarine4"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Non Cavitary", "Cavitary"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()



###Beta Diversity 

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")

#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))

#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"

#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Cavitary.,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Cavitary.",suffixes=c("",".centroid"))

newResults$Cavitary <- factor(newResults$Cavitary., levels = c("0", "1"))
centroids$Cavitary <- factor(centroids$Cavitary., levels = c("0", "1"))

#stats
adonis2(vegdist ~ newResults$Cavitary)
#plot
pdf(file = "Figures/Beta.Diversity.Bray.BAL.samples.Cavitary_NTM_pos.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Cavitary.)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("azure4", "aquamarine4")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Cavitary.), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Cavitary.)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non Cavitary", "Cavitary")), size=10) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.013")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

###Differential

#set tables to use 
ps
ps.rel
#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "Cavitary", method = "TMM")

#create results table 
ET <- exactTest(dge)
#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")

res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

# make taxa name column and name asv column 
res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, Cavitary=="0")
ps.rel.2 <- subset_samples(ps.rel, Cavitary=="1")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Cav_catg_0 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Cav_catg_1 <- meanRA.save

#clean the results df 
#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 
#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "0"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.Cav_catg_0, abundance.Cav_catg_1,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_cavitary_NTM_pos.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="0", res$abundance.Cav_catg_0, 
                        res$abundance.Cav_catg_1)

#Bubble plot
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"aquamarine4", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "azure4","grey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "aquamarine4", "azure4"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_NTM_pos_Cavitary_category_BAL_Pruned_with_contam_labels.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_NTM_pos_Cavitary_BAL_Pruned_with_contam_labels_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 

####repeat in NTM negative
cav.NTM.neg.table <- subset_samples(BAL.baseline.table.pruned.prim, NTM_pos_neg=="0")
cav.NTM.neg.table.rel <- transform_sample_counts(cav.NTM.neg.table, normalizeSample)    

ps <- subset_samples(cav.NTM.neg.table, sample_from_contralateral_noncavitary_lobe =="n")
ps <- subset_samples(ps,Cavitary != "NA")

ps.rel <- subset_samples(cav.NTM.neg.table.rel, sample_from_contralateral_noncavitary_lobe =="n")
ps.rel <- subset_samples(ps.rel,Cavitary != "NA")

#ddpcr differences 
dp.dat <- data.frame(sample_data(ps))
dp.dat$Cavitary <- as.factor(dp.dat$Cavitary)
#calculate stats for ddPCR 
my.comparisons <- compare_means(ddPCR_value~Cavitary, data = dp.dat)

# plot it
pdf(file = "Figures/ddPCR.NTM_neg.Cavitary.pdf", width=7, height=10)
ggplot(dp.dat, aes(x=Cavitary, y=ddPCR_value, fill=Cavitary))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("grey", "magenta3"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("16S rRNA Gene Copies")+
  scale_y_log10()+
  scale_x_discrete(labels = c("Non Cavitary", "Cavitary"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

##Alpha diversity
#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "Cavitary" = sample_data(ps)$Cavitary)
alpha.measures <- alpha.measures %>% 
  mutate(Cavitary= factor(Cavitary))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~Cavitary, data = alpha.measures)

# plot it
pdf(file = "Figures/Shannon.NTM_neg.Cavitary.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=Cavitary, y=Shannon, fill=Cavitary))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("grey", "magenta3"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("Non Cavitary", "Cavitary"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ Cavitary,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="Cavitary",suffixes=c("",".centroid"))
newResults$Cavitary <- factor(newResults$Cavitary, levels = c("0", "1"))
centroids$Cavitary <- factor(centroids$Cavitary, levels = c("0", "1"))
#stats 
adonis2(vegdist ~ newResults$Cavitary)
#plot
pdf(file = "Figures/Beta.Diversity.Bray.BAL.samples.Cavitary_NTM_neg.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= Cavitary)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("grey", "magenta3")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= Cavitary), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= Cavitary)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("Non Cavitary", "Cavitary")), size=10) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.009")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()



#####Differential

#set tables to use 
ps
ps.rel
#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "Cavitary", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res <- TT$table
#Reverse Directionality if you need to
#res$logFC <- res$logFC*(-1)

# make taxa name column and name asv column 
res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, Cavitary=="0")
ps.rel.2 <- subset_samples(ps.rel, Cavitary=="1")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Cav_catg_0 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.Cav_catg_1 <- meanRA.save

#clean the results df 
#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 
#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "0"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.Cav_catg_0, abundance.Cav_catg_1,
                category, contaminant) %>% 
  arrange(., desc(category))
# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_cavitary_NTM_neg.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="0", res$abundance.Cav_catg_0, 
                        res$abundance.Cav_catg_1)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"magenta3", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "grey","grey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 
bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c( "grey", "magenta3"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_NTM_neg_Cavitary_category_BAL_Pruned_with_contam_labels.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_NTM_neg_Cavitary_BAL_Pruned_with_contam_labels_legend.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 














#### Figure 5A-5D
BAL.NTM.pos.Table <- subset_samples(BAL.baseline.table.pruned.prim, NTM_pos_neg=="1")
BAL.NTM.pos.Table.rel <- transform_sample_counts(BAL.NTM.pos.Table, normalizeSample)

####categories used: 0, 1 or more exacerbations 
#add exacerbation category column 
BAL.NTM.pos.Table@sam_data$exac_catg_0_1 <- ifelse(BAL.NTM.pos.Table@sam_data$exacerbations_num >= 1, "1", "0")

#add to the relative table 
BAL.NTM.pos.Table.rel@sam_data$exac_catg_0_1 <- ifelse(BAL.NTM.pos.Table.rel@sam_data$exacerbations_num >= 1, "1", "0")

#exclude patients with exacerbations in the last 4 weeks 
ps <- subset_samples(BAL.NTM.pos.Table, ABx.4wk.prior.BAL != "1")
ps.rel <- subset_samples(BAL.NTM.pos.Table.rel, ABx.4wk.prior.BAL != "1")

###ddpcr

dp.dat <- data.frame(sample_data(ps))
dp.dat$exac_catg_0_1 <- as.factor(dp.dat$exac_catg_0_1)
dp.dat$ddPCR_value <- as.numeric(dp.dat$ddPCR_value)
#calculate stats for ddPCR 
my.comparisons <- compare_means(ddPCR_value~exac_catg_0_1, data = dp.dat)

# plot it
pdf(file = "Figures/ddPCR.NTM.pos.exac_catg_0_1.pdf", width=7, height=10)
ggplot(dp.dat, aes(x=exac_catg_0_1, y=ddPCR_value, fill=exac_catg_0_1))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("cyan4", "brown1"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("16S rRNA Gene Copies")+
  scale_y_log10()+
  scale_x_discrete(labels = c("0", "1 or more"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "exacerbations_num" = sample_data(ps)$exacerbations_num)
alpha.measures <- alpha.measures %>% 
  mutate(exac_catg = ifelse(exacerbations_num >=1, "1", "0")) %>% 
  mutate(exac_catg= factor(exac_catg))

#calculate stats for shannon
my.comparisons <- compare_means(Shannon~exac_catg, data = alpha.measures)

# plot it
pdf(file = "Figures/Shannon.NTM.pos.exacerbation_0_vs_1.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=exac_catg, y=Shannon, fill=exac_catg))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("cyan4", "brown1"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("0", "1 or more"))+
  stat_compare_means(comparisons = list(c("1", "0")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()


###Beta Diversity

#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ exac_catg_0_1,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="exac_catg_0_1",suffixes=c("",".centroid"))
#sats
adonis2(vegdist ~ newResults$exac_catg_0_1)
#plot
pdf(file = "Figures/Beta.Diversity.Bray.BAL.samples.NTM.POS.Exacerbation_category_0_1.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= exac_catg_0_1)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("cyan4", "brown1")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= exac_catg_0_1), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= exac_catg_0_1)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("0", "1 or more")), size=10) +
  ggtitle("Beta Diversity, Bray, All samples, p=0.059")+
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()


#####Differential

#set tables to use 
ps 
ps.rel
#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "exac_catg_0_1", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res <- TT$table
#Reverse Directionality if you need to
res$logFC <- res$logFC*(-1)

# make taxa name column and name asv column 
res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, exac_catg_0_1=="1")
ps.rel.2 <- subset_samples(ps.rel, exac_catg_0_1=="0")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.exac_catg_0 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.exac_catg_1 <- meanRA.save

#clean the results df 
#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]
#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))
# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "0"))
res <- inner_join(res, contamlist, by="asv")
#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.exac_catg_0, abundance.exac_catg_1,
                category, contaminant) %>% 
  arrange(., desc(category))

# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_Exacerbations_category_BAL_Pruned_NTM.pos.exac_categ_0_1.csv")

#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.exac_catg_1, 
                        res$abundance.exac_catg_0)

#Bubble plot
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"brown1", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "cyan4","grey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))

#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 
bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("brown1", "cyan4"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_NTM_pos_Exacerbation_category_BAL_Pruned_with_contam_labels_exac_categ_0_1.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_NTM_pos_Exacerbation_category_BAL_Pruned_with_contam_labels_legend_exac_categ_0_1.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 



###########################
######repeat for NTM - 
BAL.NTM.neg.Table <- subset_samples(BAL.baseline.table.pruned.prim, NTM_pos_neg=="0")
BAL.NTM.neg.Table.rel <- transform_sample_counts(BAL.NTM.neg.Table, normalizeSample)

#add exacerbation category column 
BAL.NTM.neg.Table@sam_data$exac_catg_0_1 <- ifelse(BAL.NTM.neg.Table@sam_data$exacerbations_num >= 1, "1", "0")
#add to the relative table 
BAL.NTM.neg.Table.rel@sam_data$exac_catg_0_1 <- ifelse(BAL.NTM.neg.Table.rel@sam_data$exacerbations_num >= 1, "1", "0")

#exclude patients with exacerbations in the last 4 weeks 
ps <- subset_samples(BAL.NTM.neg.Table, ABx.4wk.prior.BAL != "1")
ps.rel <- subset_samples(BAL.NTM.neg.Table.rel, ABx.4wk.prior.BAL != "1")


###ddpcr
dp.dat <- data.frame(sample_data(ps))
dp.dat$exac_catg_0_1 <- as.factor(dp.dat$exac_catg_0_1)
dp.dat$ddPCR_value <- as.numeric(dp.dat$ddPCR_value)
#calculate stats for ddPCR 
my.comparisons <- compare_means(ddPCR_value~exac_catg_0_1, data = dp.dat)
# plot it
pdf(file = "Figures/ddPCR.NTM.neg.exac_catg_0_1.pdf", width=7, height=10)
ggplot(dp.dat, aes(x=exac_catg_0_1, y=ddPCR_value, fill=exac_catg_0_1))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkslategray3", "chocolate"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("16S rRNA Gene Copies")+
  scale_y_log10()+
  scale_x_discrete(labels = c("0", "1 or more"))+
  stat_compare_means(comparisons = list(c("0", "1")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()

#########Alpha diversity
#calculate alpha diversity 
alpha.measures <- data.frame("Shannon" = estimate_richness(ps, measures = "Shannon"), 
                             "Observed" = estimate_richness(ps, measures = "Observed"), 
                             "Simpson" = estimate_richness(ps, measures = "Simpson"),
                             "exacerbations_num" = sample_data(ps)$exacerbations_num)
alpha.measures <- alpha.measures %>% 
  mutate(exac_catg = ifelse(exacerbations_num >=1, "1", "0")) %>% 
  mutate(exac_catg= factor(exac_catg))
#calculate stats for shannon
my.comparisons <- compare_means(Shannon~exac_catg, data = alpha.measures)
# plot it
pdf(file = "Figures/Shannon.NTM.neg.exacerbation_0_vs_1.pdf", width=7, height=10)
ggplot(alpha.measures, aes(x=exac_catg, y=Shannon, fill=exac_catg))+
  geom_boxplot(outlier.shape = NA)+
  stat_boxplot(geom = "errorbar", width= 0.15)+
  scale_fill_manual(values = c("darkslategray3", "chocolate"))+
  geom_jitter(color="black", size=0.4, alpha=0.7)+
  xlab("")+ylab("Shannon")+
  scale_x_discrete(labels = c("0", "1 or more"))+
  stat_compare_means(comparisons = list(c("1", "0")), 
                     method = "wilcox.test", label = "p.value", size=10)+
  theme_classic()+
  theme(text = element_text(size = 16), 
        axis.title.x = element_text(size = 24, face = "bold"), 
        axis.text.x = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20, face = "bold"), 
        axis.title.y = element_text(size = 24, face = "bold"))+
  guides(fill="none")
dev.off()
###Beta Diversity
#Create Distance Matrix
vegdist = vegdist(t(otu_table(ps.rel)), method = "bray")
#Formulate principal component co-ordinates for PCOA plot, k as the choice of PCs
CmdScale <- cmdscale(vegdist, k =10)
#calculated Sample variance for each PC
vars <- apply(CmdScale, 2, var)
#Create Variable with the Percent Variance
percentVar <- round(100 * (vars/sum(vars)))
#Merge PC Data with MetaData
newResults <- merge(x = CmdScale, y = sample_data(ps.rel), by = "row.names", all.x = TRUE)
#Rename Variables for PC1 and PC2
colnames(newResults)[colnames(newResults)=="V1"] <- "PC1"
colnames(newResults)[colnames(newResults)=="V2"] <- "PC2"
colnames(newResults)[colnames(newResults)=="Row.names"] <- "name"
#Calculate the Centroid Value
centroids <- aggregate(cbind(PC1,PC2)~ exac_catg_0_1,data= newResults, mean)
#Merge the Centroid Data into the PCOA Data
newResults <- merge(newResults,centroids,by="exac_catg_0_1",suffixes=c("",".centroid"))
#perform adonis and add p value 
x=adonis2(vegdist ~ newResults$exac_catg_0_1)
p_value<-x$`Pr(>F)`[[1]]
subtitle_output<- paste0("Adonis, p=",p_value)
#plot
pdf(file = "Figures/Beta.Diversity.Bray.BAL.samples.NTM.neg.Exacerbation_category_0_1.pdf", 
    height = 20, width = 18)
ggplot(newResults, aes(PC1, PC2, color= exac_catg_0_1)) + # Graph PC1 and PC2
  geom_point(size=5) + # Set the size of the points
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + #Label PC1
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + #Label PC2
  #Set colors for each category, should be the same number of categories
  scale_color_manual(values=c("darkslategray3", "chocolate")) + 
  #plot point and lines from centroid
  geom_point(data=centroids, aes(x=PC1, y=PC2, color= exac_catg_0_1), size=0) +
  geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color= exac_catg_0_1)) + 
  #labels centroids should be same number of categories
  geom_label_repel(data = centroids, aes(x=PC1, y=PC2, label=c("0", "1 or more")), size=10) +
  labs(title="Beta Diversity", 
       subtitle=subtitle_output) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_line(linetype = "dashed", size = 0.5, colour = "grey80"),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.title=element_text(size=30,face="bold"),axis.text.x=element_text(colour = "grey80", size = rel(0.75)),axis.text.y=element_text(colour = "grey80", size = rel(0.75)),axis.ticks=element_blank(),plot.margin=unit(c(1,1,1,1),"line"), legend.position="none")
dev.off()

######Differential
#set tables to use 
ps 
ps.rel 
#prepare data by converting phyloseq to edgeR and performing analysis. This function from Physloeq package should give the solution 
#run the function and get edgeR object
dge <- phyloseq_to_edgeR(physeq = ps, group = "exac_catg_0_1", method = "TMM")

#create results table 
#run test 
ET <- exactTest(dge)

#get results 
TT <- topTags(ET,  n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
res <- TT$table
#Reverse Directionality if you need to
res$logFC <- res$logFC*(-1)
# make taxa name column and name asv column 
res$names <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",res$Species,sep=""),
                          ifelse(!is.na(res$Genus), paste("g_",res$Genus,sep=""),
                                 ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",sep=""),
                                        ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",sep=""),
                                               ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",sep=""),
                                                      ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",sep=""),
                                                             ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(res))))))))))

#repeat for nameASV 
res$nameASV <- paste(ifelse(!is.na(res$Species), paste("g_",res$Genus,"_s_",rownames(res),sep=""),
                            ifelse(!is.na(res$Genus), paste("g_",res$Genus,"_",rownames(res),sep=""),
                                   ifelse(!is.na(res$Family), paste("f_",res$Family,"_g__",rownames(res),sep=""),
                                          ifelse(!is.na(res$Order), paste("o_",res$Order, "_f__g__",rownames(res),sep=""),
                                                 ifelse(!is.na(res$Class), paste("c_",res$Class, "_o__f__g__",rownames(res),sep=""),
                                                        ifelse(!is.na(res$Phylum), paste("p_",res$Phylum, "_c__o__f__g__",rownames(res),sep=""),
                                                               ifelse(!is.na(res$Domain), paste("k_",res$Domain, "_p__c__o__f__g__",rownames(res),sep=""), paste(rownames(res))))))))))

#define factors 
res <- res %>% 
  dplyr::mutate(names = factor(names, levels= unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV)))
res$nameASV.trim <- paste0(res$names,paste0(".."), paste0(str_sub(rownames(res), - 3, - 1)))

#get relative abundance data 
ps.rel.1 <- subset_samples(ps.rel, exac_catg_0_1=="1")
ps.rel.2 <- subset_samples(ps.rel, exac_catg_0_1=="0")

RA.df.1 <- data.frame(otu_table(ps.rel.1))
meanRA <- rowMeans(RA.df.1)
meanRA.save <- meanRA[rownames(res)]
res$abundance.exac_catg_0 <- meanRA.save

#repeat for the other group
RA.df.2 <- data.frame(otu_table(ps.rel.2))
meanRA <- rowMeans(RA.df.2)
meanRA.save <- meanRA[rownames(res)]
res$abundance.exac_catg_1 <- meanRA.save

#clean the results df 

#Remove Any Data without LOGFC data
res <- res[!is.na(res$logFC),]

#Convert Important columns to Numeric
res$FDR <-   as.numeric(as.character(res$FDR))
res$logFC <-       as.numeric(as.character(res$logFC))

# drop any NA values 
res<- res %>% 
  drop_na(., FDR) %>% 
  dplyr::arrange(desc(FDR))

#get abundance of each category alone to save it in the table 

#add contamlist 
res<- res %>% 
  mutate(asv=rownames(.)) %>% 
  mutate(category= ifelse(logFC>0, "1", "0"))

res <- inner_join(res, contamlist, by="asv")

#get results to save in the table 
edgeR.to.save <- res %>% 
  dplyr::select(asv, Kingdom, Phylum, Class, Order, Family, Genus, Species, 
                nameASV, logFC, PValue, FDR, abundance.exac_catg_0, abundance.exac_catg_1,
                category, contaminant) %>% 
  arrange(., desc(category))
# save results
write.csv(edgeR.to.save, file = "Results/edgeR.results_Exacerbations_category_BAL_Pruned_NTM.neg.exac_categ_0_1.csv")
#combine abundance into one column to make it easy to plot later
res$abundance <- ifelse(res$category=="1", res$abundance.exac_catg_1, 
                        res$abundance.exac_catg_0)

#Bubble plot#
res.bubble <- res

#get top taxa upregaulated
Sig.Up.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(desc(logFC)) %>% 
  dplyr::slice(., 1:25)

#get top taxa downregulated 
Sig.Down.Taxa <- res.bubble %>% 
  filter(., FDR < 0.2) %>% 
  filter(., logFC != 0) %>% 
  arrange(logFC) %>% 
  dplyr::slice(., 1:25)

#combine two dataframes to plot 
res.bubble.filtered <- dplyr::bind_rows(Sig.Up.Taxa, Sig.Down.Taxa)

#update color vector 
res.bubble.filtered$col <- ifelse(res.bubble.filtered$logFC>0&res.bubble.filtered$FDR<0.2,"chocolate", 
                                  ifelse(res.bubble.filtered$logFC<0&res.bubble.filtered$FDR<0.2, "darkslategray3","grey"))

#Now add ASV trimmed 
res.bubble.filtered$nameASV.trim <- paste0(res.bubble.filtered$names,paste0(".."), 
                                           paste0(str_sub(res.bubble.filtered$asv, - 3, - 1)))

res.bubble.filtered <- res.bubble.filtered %>% 
  dplyr::mutate(names = factor(names, levels = unique(names))) %>% 
  dplyr::mutate(nameASV = factor(nameASV, levels = unique(nameASV))) %>% 
  dplyr::mutate(nameASV.trim = factor(nameASV.trim, levels = unique(nameASV.trim))) %>% 
  dplyr::mutate(nameASV.trim.colored=nameASV.trim) %>% 
  dplyr::arrange(., logFC) %>% 
  mutate(ord = order(logFC))


#add color of potential contaminants 
res.bubble.filtered$color <- ifelse(res.bubble.filtered$contaminant == "TRUE", "red", "black")
res.bubble.filtered$nameASV.trim.colored <- paste0("<span style=\"color: ", res.bubble.filtered$color, "\">", res.bubble.filtered$nameASV.trim.colored, "</span>")

#if statistically significant --> add segment line 
#size based on relative abundance. Color based on category 

bubble.plot <- ggplot(res.bubble.filtered, mapping = aes(x=logFC, y=fct_reorder(nameASV.trim.colored, ord), fill=col, size=abundance))+
  geom_point(color="black",alpha=0.8,shape=21)+
  geom_segment(data = res.bubble.filtered[res.bubble.filtered$FDR < 0.2,],
               aes(yend=fct_reorder(nameASV.trim.colored, ord)),xend=(-30), color="darkgrey", linetype="solid", size=1)+  
  scale_size_continuous(name = "Relative Abundance", range = c(5,20))+
  scale_fill_manual(values = c("chocolate", "darkslategray3"))+
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major.y = element_line(colour = "#EBEBEB",linetype="dashed"),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x=element_text(colour="black", size=18, face="bold"),
        axis.text.y=element_markdown(face="bold",size=16),
        axis.ticks=element_line(colour="black"),
        legend.background = element_rect(color=NA), 
        plot.title=element_text(size=23, face="bold"))+
  xlab("") +
  ylab("")+
  ggtitle("EdgeR")+
  geom_vline(xintercept=0, color="red",linetype="dashed")+
  guides(fill="none")

pdf(file="Figures/edgeR_Bubble_NTM_neg_Exacerbation_category_BAL_Pruned_with_contam_labels_exac_categ_0_1.pdf", height = 11, width = 8)
bubble.plot+
  guides(size="none")
dev.off()

#get seperate legend plot 
my_legend <- get_legend(bubble.plot)
#save it 
pdf(file="Figures/edgeR_Bubble_NTM_neg_Exacerbation_category_BAL_Pruned_with_contam_labels_legend_exac_categ_0_1.pdf", 
    height = 4, width = 2)
as_ggplot(my_legend)
dev.off() 










