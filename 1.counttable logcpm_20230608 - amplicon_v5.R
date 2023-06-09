#Clears the workspace
rm(list = ls())

setwd('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_amplicon_newaligned/')

#loads the required libraries, but suppresses the output to the screen 
suppressMessages(lapply(c('edgeR','MASS' ,'stringr', 'ggplot2', 'pheatmap', 'writexl',
                          'cowplot', 'reshape2'), require, character.only = TRUE))

#writes the session information to an output file
sink("sessionInfo.txt")
sessionInfo()
sink()

#read in the counts matrix 
countData <- read.delim('all_sample_featureCounts.tsv/all_sample_featureCounts.tsv',row.names = 'feature') #read rawcounts


#fix the colnames
colnames(countData) <- gsub('.counts_fC', '',colnames(countData) )
colnames(countData) <- gsub('sample_', '',colnames(countData) )
colnames(countData)[c(2,5:12)] <- paste0('0', colnames(countData)[c(2,5:12)])
countData <- countData[,order(colnames(countData))]

###########
countData <- countData[,1:10] #remove the control samples
colnames(countData) <- paste0(colnames(countData), 'R')
colnames(countData) <- paste0('0', colnames(countData))
###############


#read in metadata and check it's matched
meta <- read.delim('../input/metadata_ampliseq_v2.txt') #read in metadata


#check if the meta and count are matched and then add the IDs
colnames(countData) == meta[,1]
colnames(countData) <- meta[,6]
meta$TotalCounts <- colSums(countData) * 1e-3
meta$TotalGenes <- colSums(countData > 1) #set to minimum 10 reads/sample


#check if the genes are known
dir.create('output')

#make barplots of expression and number of genes of all reads
write.table(as.matrix(colSums(countData)), 'output/20230608_counts_all.txt', sep = '\t')
write.table(countData, 'output/20230608_countData_all.txt', sep = '\t', row.names = T)


pdf('output/20230608_barplots counts_allreads.pdf', width = 8, height = 5)
barplot(colSums(countData) * 1e-3, names = colnames(countData), ylab = "Library Size (thousands)")
dev.off()



pdf('output/20230512_barplots genes_allreads.pdf', width = 8, height = 5)
barplot(colSums(countData > 1), names = colnames(countData), ylab = "Genes detected")
dev.off()


#filter on genes to keep
#make an empty frame for number of genes
genes_retained <- matrix(ncol=3, nrow=10)
genes_retained[,1] <- seq(0, 0.9, 0.1)


#make a loop to go through the options
j = 1 #initialize the row
minCount = 10 #minimum number of counts needed to be kept
for(i in seq(0, 0.9, 0.1)){
    
high <- countData > minCount #Is there a gene present
ret <- sum((rowSums(high)/ncol(countData)) > i) # genes that are expressed in XX% of all samples
idx <- (rowSums(high)/ncol(countData)) > i #save the genes that are expressed in >20% of all samples
countData2 <- countData[which(idx),]

genes_retained[j,2] <- ret #store in the matrix
genes_retained[j,3] <- sum(countData2) #store in the matrix
j = j+1

}

#make into a plot
genes_retained <- as.data.frame(genes_retained)
colnames(genes_retained) <- c('Cutoff', 'Genes', 'Reads')
ggplot(genes_retained, aes(x=Cutoff*100, y=Genes)) + 
    geom_bar(stat = 'identity') + 
    ylab('Number of detected genes (thousands)') + 
    xlab('Present in % of samples')
ggsave('output/20230512_genes_filtered_allreads.pdf', width = 8, height = 5)


#make a new count table with retained genes
co1 <- 0.1
co2 <- 0.1
minCount = 100
high1 <- countData[,1:6] > minCount #Is there a gene present
high2 <- countData[,7:10] > minCount #Is there a gene present

id1 <- which((rowSums(high1)/6) > co1) #save the genes that are expressed in >20% of all samples
id2 <- which((rowSums(high2)/4) > co2) #save the genes that are expressed in >20% of all samples

idx <- unique(c(id1, id2))
countData2 <- countData[idx,]

nrow(countData2) #number of genes after filtering


setwd('output')
pdf('barplots counts_co_filteredreads.pdf', width = 8, height = 5)
barplot(colSums(countData2) * 1e-3, names = colnames(countData2), ylab = "Library Size (thousands)")
dev.off()


pdf('barplots genes detect_filteredreads.pdf', width = 8, height = 5)
barplot(colSums(countData2 > 0), names = colnames(countData), ylab = "Genes detected")
dev.off()


#write the count data version 2
write.table(countData2, '20230608_countData_filtered.txt', sep = '\t', row.names = T)


#make into logcpm and remove the ones not in ampliseq
all2 <- logcpm <- cpm(countData2, prior.count=1, log=TRUE)


#do overall PCA on the data
data.pca <- prcomp(t(logcpm), center=T, scale=T)
summary(data.pca)                 


meta.pca <- cbind(meta, data.pca$x[,1:2])
ggplot(data=meta.pca, aes(x=PC1, y=PC2)) +
    geom_point(aes(colour=Group)) + 
  geom_text(aes(label = ID), nudge_y = 3)  + 
  theme_bw()
ggsave('PCA_newfilteredreads_based on counts.pdf', width = 8, height = 7)
ggsave('PCA_newfilteredreads_based on counts.pdf', width = 3, height = 2)



write.table(logcpm, 'logCPM_amplicon_sel.txt', sep = '\t', row.names = T)


dev.off()
pdf('heatmap_logcpm.pdf', width=5, height = 15)
pheatmap(as.matrix(logcpm))
dev.off()


#get data per gene out
comparisons <- matrix(nrow = nrow(all2), ncol=2)
rownames(comparisons) <- rownames(all2)

for(i in 1:nrow(all2)){
    
    comparisons[i,1] <- mean(as.numeric(all2[i,meta$Group == 'case'])) - 
        mean(as.numeric(all2[i,meta$Group == 'control']))  
    comparisons[i,2] <- t.test(as.numeric(all2[i,meta$Group == 'case']),  
        as.numeric(all2[i,meta$Group == 'control']))[[3]]  
}

comparisons <- data.frame(comparisons)
colnames(comparisons) <- c('Log2FC_case_ctrl', 'pval')
comparisons$padj <- p.adjust(comparisons$pval, method = 'fdr')
comparisons$Sig <- 'No'
comparisons$Sig[comparisons$padj < 0.05] <- ' Yes'
comparisons$name2 <- comparisons$name <- rownames(comparisons)
comparisons$name[comparisons$padj > 0.05] <- NA
comparisons$name[comparisons$Log2FC_case_ctrl > 0 & 
                   comparisons$padj > 0.001  ] <- NA


ggplot(comparisons, aes(x=Log2FC_case_ctrl, y = -log10(padj))) + 
    geom_point(aes(colour = Sig)) + 
    geom_label(aes(label=name), nudge_y = 0.1) + 
  scale_colour_manual(values = c('red', 'grey')) + 
  theme_bw()  + geom_vline(xintercept = 0, linetype  = 'dashed')+ 
  geom_hline(yintercept = -log10(0.05), linetype  = 'dashed')
ggsave('volcano_plot_2.pdf', width = 7, height = 7)


table(comparisons$Log2FC_case_ctrl>0 & comparisons$padj < 0.05)
table(comparisons$Log2FC_case_ctrl<0 & comparisons$padj < 0.05)



ggplot(comparisons, aes(x=Log2FC_case_ctrl, y = -log10(padj))) + 
  geom_point(aes(colour = Sig)) + 
  scale_colour_manual(values = c('red', 'grey')) + 
  theme_bw()  + geom_vline(xintercept = 0, linetype  = 'dashed')+ 
  geom_hline(yintercept = -log10(0.05), linetype  = 'dashed')
ggsave('volcano_plot_2.pdf', width = 5, height = 5)

write_xlsx(comparisons, 'comparisons_volcano.xlsx')
    
#select top genes
names <- comparisons$name2[comparisons$padj < 0.002 & abs(comparisons$Log2FC_case_ctrl)>2]

names1 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\BTM_source_ES.xlsx', 1)
names11 <- names1[names1$pathway == 'enriched in monocytes (II) (M11.0)',]$leadingEdge2
names12 <-do.call(c,strsplit(names11, split = '/'))
names2 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\BTM_source_ES.xlsx', 2)
names21 <- names2[names2$pathway == 'enriched in monocytes (II) (M11.0)',]$leadingEdge2
names22 <-do.call(c,strsplit(names21, split = '/'))
names3 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\BTM_source_ES.xlsx', 3)
names31 <- names3[names3$pathway == 'enriched in monocytes (II) (M11.0)',]$leadingEdge2
names32 <-do.call(c,strsplit(names31, split = '/'))
names4 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\BTM_source_ES.xlsx', 4)
names41 <- names4[names4$pathway == 'enriched in monocytes (II) (M11.0)',]$leadingEdge2
names42 <-do.call(c,strsplit(names41, split = '/'))
names_all <- c(names12, names22, names32, names42)
names_mono <- names(which(table(names_all) > 3))

select <- all2[rownames(all2) %in% names_mono,]
select_melt <- melt(as.matrix(select), id = 0)
colnames(select_melt) <- c('Gene', 'Sample', 'LogCPM')
select_melt2 <- merge(select_melt, meta, by.x='Sample', by.y='ID')

ggplot(select_melt2, aes(x=Group, y=LogCPM))  + 
    geom_boxplot(outlier.color = NA, aes(fill = Group)) + 
    geom_jitter(width = 0.1, height = 0) + 
    facet_wrap(~Gene, ncol=8) +
  scale_colour_brewer(palette="Set3") + theme_bw()
ggsave('boxplots_selectedgenes_mono.pdf', width = 8, height = 2)




names1 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\Reactome_source.xlsx', 1)
names11 <- names1[names1$pathway == 'Complement cascade',]$leadingEdge2
names12 <-do.call(c,strsplit(names11, split = '/'))
names2 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\Reactome_source.xlsx', 2)
names21 <- names2[names2$pathway == 'Complement cascade',]$leadingEdge2
names22 <-do.call(c,strsplit(names21, split = '/'))
names3 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\Reactome_source.xlsx', 3)
names31 <- names3[names3$pathway == 'Complement cascade',]$leadingEdge2
names32 <-do.call(c,strsplit(names31, split = '/'))
names4 <- readxl::read_xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_ssGSEA_ampliseq\\Reactome_source.xlsx', 4)
names41 <- names4[names4$pathway == 'Complement cascade',]$leadingEdge2
names42 <-do.call(c,strsplit(names41, split = '/'))
names_all <- c(names12, names22, names32, names42)
names_comp <- names(which(table(names_all) > 3))

select <- all2[rownames(all2) %in% names_comp,]
select_melt <- melt(as.matrix(select), id = 0)
colnames(select_melt) <- c('Gene', 'Sample', 'LogCPM')
select_melt2 <- merge(select_melt, meta, by.x='Sample', by.y='ID')

ggplot(select_melt2, aes(x=Group, y=LogCPM))  + 
  geom_boxplot(outlier.color = NA, aes(fill = Group)) + 
  geom_jitter(width = 0.1, height = 0) + 
  facet_wrap(~Gene, ncol=11) +
  scale_colour_brewer(palette="Set3") + theme_bw()
ggsave('boxplots_selectedgenes_complement.pdf', width = 11, height = 2)



