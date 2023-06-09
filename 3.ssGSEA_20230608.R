rm(list = ls())

library('reshape2')
library('ggplot2')
library('Hmisc')
library('fgsea')

setwd('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230523_ssGSEA_ampliseq')

counts <- read.delim('../20230608_amplicon_newaligned/output/logCPM_amplicon_sel.txt')
meta <- read.delim('../input/metadata_ampliseq_v2.txt')
rownames(meta) <- meta$ID

identical(meta$ID, colnames(counts))

#read in the BTMs and gene annotation keys
genes <- read.delim('../input/genes.txt')
genes <- genes[,c('Geneid_clean', 'GeneSymbol')]
BTMs <- gmtPathways("../input/BTM_Li_2013.gmt")
Reactome <- gmtPathways("../input/ReactomePathways.gmt")
BP <- gmtPathways("../input/c5.bp.v6.2.symbols.gmt")
KEGG <- gmtPathways("../input/c2.cp.kegg.v2023.1.Hs.symbols.gmt")

setwd('../20230608_ssGSEA_ampliseq/')

#obtain mean expression of genes for all baseline samples
ctrl_0 <- subset(meta, Group == 'control') 
counts_c0 <- counts[,colnames(counts) %in% ctrl_0$ID]
control <- rowMeans(counts_c0)

#make the FC per sample to all baseline samples
post <- subset(meta, Group =='case')
counts_post <- counts[,colnames(counts) %in% post$ID]
FC <- counts_post - control 


#make a GSEA analysis per sample for either BTM or reactome
################################
#BTM
##########################################
GSEA_list = list()
for(i in 1:ncol(FC)){
  sample <- data.frame('Fold' = FC[,i])
  sample$Gene <- rownames(FC)
  sample <- merge(sample, genes, by.x='Gene', by.y='GeneSymbol')
  sample <- sample[!duplicated(sample$Gene),]
  sample <- sample[order(sample$Fold, decreasing = T),]
  FC_List <- sample[,2]
  names(FC_List) <- sample[,1]
  fgseaRes <- fgseaMultilevel(BTMs, stats=FC_List, minSize = 10)
  GSEA_list[[i]] <- fgseaRes
  assign(paste0('GSEA_', i), GSEA_list[[i]])
}

#make the list into a dataframe
gseacolumns <- do.call(cbind, GSEA_list)
keep <- which(grepl('NES', names(gseacolumns)))

#extract the NES scores per sample and add the sample ID
NES <- gseacolumns[,c(1,5,13,21,29)]
colnames(NES) <- c('Pathway', as.character(post$ID))
NES.m <- melt(NES)
NES.m$Abs <- abs(NES.m$value)
ggplot(NES.m, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))

#extract the p-values per sample and add sample ID as column name
Padj <- gseacolumns[,c(1,3,11,19,27)] 
colnames(Padj) <- c('Pathway', as.character(post$Sample))

#filter on significant pathways that are present in at least 2 per timepoint
#arrange in correct order
NES[Padj>0.05] <- 0 #set all the non-significant ones to 0
NES$Pathway <- Padj$Pathway #re-add the pathways

NES_s <- NES[rowSums(Padj < .05, na.rm = T)>1,] #filter on at least two significant 1


#melt and create absolute value for size of point
NES.m2 <- melt(NES_s)
NES.m2$Abs <- abs(NES.m2$value)
NES.BTM <- NES.m2

#make and save plot
ggplot(NES.m2, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))+ 
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",high = "red")+
  xlab('Donor') + 
  theme_bw()
ggsave('ssGSEA_BTM_2023_06_08_ES.pdf', width=10, height=10)

for(i in 1:length(GSEA_list)){
  tip <- GSEA_list[[i]]
  tip$leadingEdge2 <- 'NA'
  for(j in 1:nrow(tip)){
      tip$leadingEdge2[j] <- paste(unlist(tip$leadingEdge[j]), collapse='/')
  }
  GSEA_list[[i]] <- tip[,-8]
  }

#save source data
writexl::write_xlsx(GSEA_list, 'BTM_source_ES.xlsx')



for(i in 1:ncol(FC)){
  sample <- data.frame('Fold' = FC[,i])
  sample$Gene <- rownames(FC)
  sample <- merge(sample, genes, by.x='Gene', by.y='GeneSymbol')
  sample <- sample[!duplicated(sample$Gene),]
  sample <- sample[order(sample$Fold, decreasing = T),]
  FC_List <- sample[,2]
  names(FC_List) <- sample[,1]
  pl1 <- plotEnrichment(BTMs[["enriched in monocytes (II) (M11.0)"]],
                        FC_List) + labs(title="enriched in monocytes (II) (M11.0)")
  hits <- data.frame(x = which(sample$Gene %in% BTMs[["enriched in monocytes (II) (M11.0)"]]))
  hits$y <- -0.1
  assign(paste0('hits', i), hits)
  if(i==1){
  dat <- pl1$data
  dat$sample <- colnames(FC)[i]
  }else{
    dat2 <- pl1$data
    dat2$sample <- colnames(FC)[i]
    dat <- rbind(dat, dat2)
  }
}
ggplot(dat, aes(x=x, y=y, colour = sample)) + geom_line() + 
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y), data=hits1, 
            height = 0.025, colour = 'red') +
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y-0.04), data=hits2, 
            height = 0.025, colour = 'darkgreen')+ 
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y-0.08), data=hits3, 
            height = 0.025, colour = 'cyan') +
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y-0.12), data=hits4, 
            height = 0.025, colour = 'purple') + 
  theme_bw() + ggtitle('enriched in monocytes (II) (M11.0)') + 
  scale_colour_manual(values = c('red', 'darkgreen', 'cyan', 'purple')) + 
  xlab('Gene rank') + ylab('enrichment score')
ggsave('BTM_monocyte enrichment.pdf', width=5, height=5)





#make a GSEA analysis per sample for either BTM or reactome
################################
#Reactome
##########################################
GSEA_list = list()
for(i in 1:ncol(FC)){
  sample <- data.frame('Fold' = FC[,i])
  sample$Gene <- rownames(FC)
  sample <- merge(sample, genes, by.x='Gene', by.y='GeneSymbol')
  sample <- sample[!duplicated(sample$Gene),]
  sample <- sample[order(sample$Fold, decreasing = T),]
  FC_List <- sample[,2]
  names(FC_List) <- sample[,1]
  fgseaRes <- fgseaMultilevel(Reactome, stats=FC_List, minSize = 10)
  GSEA_list[[i]] <- fgseaRes
  assign(paste0('GSEA_', i), GSEA_list[[i]])
}

#make the list into a dataframe
gseacolumns <- do.call(cbind, GSEA_list)
keep <- which(grepl('NES', names(gseacolumns)))

#extract the NES scores per sample and add the sample ID
NES <- gseacolumns[,c(1,5,13,21,29)]
colnames(NES) <- c('Pathway', as.character(post$ID))
NES.m <- melt(NES)
NES.m$Abs <- abs(NES.m$value)
ggplot(NES.m, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))

#extract the p-values per sample and add sample ID as column name
Padj <- gseacolumns[,c(1,3,11,19,27)] 
colnames(Padj) <- c('Pathway', as.character(post$Sample))

#filter on significant pathways that are present in at least 2 per timepoint
#arrange in correct order
NES[Padj>0.05] <- 0 #set all the non-significant ones to 0
NES$Pathway <- Padj$Pathway #re-add the pathways

NES_s <- NES[rowSums(Padj < .05, na.rm = T)>1,] #filter on at least two significant 1


#melt and create absolute value for size of point
NES.m2 <- melt(NES_s)
NES.m2$Abs <- abs(NES.m2$value)
NES.reac <- NES.m2


#make and save plot
ggplot(NES.m2, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))+ 
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",high = "red")+
  xlab('Donor') + 
  theme_bw()
ggsave('ssGSEA_Reactome_2023_06_08.pdf', width=5, height=5)


NES.all <- rbind(NES.BTM, NES.reac)
NES.all$Pathway <- gsub('Regulation of Complement cascade', ' Regulation of Complement cascade', NES.all$Pathway)
NES.all$Pathway <- gsub('Complement cascade', ' Complement cascade', NES.all$Pathway)


ggplot(NES.all, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))+ 
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",high = "red")+
  xlab('Donor') + 
  theme_bw()
ggsave('ssGSEA_Reactome_KEGG_2023_06_08.pdf', width=10, height=5)





for(i in 1:length(GSEA_list)){
  tip <- GSEA_list[[i]]
  tip$leadingEdge2 <- 'NA'
  for(j in 1:nrow(tip)){
    tip$leadingEdge2[j] <- paste(unlist(tip$leadingEdge[j]), collapse='/')
  }
  GSEA_list[[i]] <- tip[,-8]
}

#save source data
writexl::write_xlsx(GSEA_list, 'Reactome_source.xlsx')



for(i in 1:ncol(FC)){
  sample <- data.frame('Fold' = FC[,i])
  sample$Gene <- rownames(FC)
  sample <- merge(sample, genes, by.x='Gene', by.y='GeneSymbol')
  sample <- sample[!duplicated(sample$Gene),]
  sample <- sample[order(sample$Fold, decreasing = T),]
  FC_List <- sample[,2]
  names(FC_List) <- sample[,1]
  pl1 <- plotEnrichment(Reactome[["Complement cascade"]],
                        FC_List) + labs(title="enriched in monocytes (II) (M11.0)")
  hits <- data.frame(x = which(sample$Gene %in% Reactome[["Complement cascade"]]))
  hits$y <- -0.6
  assign(paste0('hits', i), hits)
  if(i==1){
    dat <- pl1$data
    dat$sample <- colnames(FC)[i]
  }else{
    dat2 <- pl1$data
    dat2$sample <- colnames(FC)[i]
    dat <- rbind(dat, dat2)
  }
}
ggplot(dat, aes(x=x, y=y, colour = sample)) + geom_line() + 
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y), data=hits1, 
            height = 0.025, colour = 'red') +
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y-0.04), data=hits2, 
            height = 0.025, colour = 'darkgreen')+ 
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y-0.08), data=hits3, 
            height = 0.025, colour = 'cyan') +
  geom_tile(inherit.aes=F, mapping = aes(x=x, y=y-0.12), data=hits4, 
            height = 0.025, colour = 'purple') + 
  theme_bw() + ggtitle('Complement cascade (Reactome)') + 
  scale_colour_manual(values = c('red', 'darkgreen', 'cyan', 'purple')) + 
  xlab('Gene rank') + ylab('enrichment score')
ggsave('Complement cascade.pdf', width=5, height=5)







#make a GSEA analysis per sample for either BTM or reactome
################################
#KEGG
##########################################
GSEA_list = list()
for(i in 1:ncol(FC)){
  sample <- data.frame('Fold' = FC[,i])
  sample$Gene <- rownames(FC)
  sample <- merge(sample, genes, by.x='Gene', by.y='GeneSymbol')
  sample <- sample[!duplicated(sample$Gene),]
  sample <- sample[order(sample$Fold, decreasing = T),]
  FC_List <- sample[,2]
  names(FC_List) <- sample[,1]
  fgseaRes <- fgseaMultilevel(KEGG, stats=FC_List, minSize = 10)
  GSEA_list[[i]] <- fgseaRes
  assign(paste0('GSEA_', i), GSEA_list[[i]])
}

#make the list into a dataframe
gseacolumns <- do.call(cbind, GSEA_list)
keep <- which(grepl('NES', names(gseacolumns)))

#extract the NES scores per sample and add the sample ID
NES <- gseacolumns[,c(1,6,14,22,30)]
colnames(NES) <- c('Pathway', as.character(post$ID))
NES.m <- melt(NES)
NES.m$Abs <- abs(NES.m$value)
ggplot(NES.m, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))

#extract the p-values per sample and add sample ID as column name
Padj <- gseacolumns[,c(1,3,11,19,27)] 
colnames(Padj) <- c('Pathway', as.character(post$Sample))

#filter on significant pathways that are present in at least 2 per timepoint
#arrange in correct order
NES[Padj>0.05] <- 0 #set all the non-significant ones to 0
NES$Pathway <- Padj$Pathway #re-add the pathways

NES_s <- NES[rowSums(Padj < .05, na.rm = T)>0,] #filter on at least one significant 1


#melt and create absolute value for size of point
NES.m2 <- melt(NES_s)
NES.m2$Abs <- abs(NES.m2$value)

#make and save plot
ggplot(NES.m2, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))+ 
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",high = "red")+
  xlab('Donor') + 
  theme_bw()
ggsave('ssGSEA_KEGG_2023_05_23.pdf', width=5, height=5)

for(i in 1:length(GSEA_list)){
  tip <- GSEA_list[[i]]
  tip$leadingEdge2 <- 'NA'
  for(j in 1:nrow(tip)){
    tip$leadingEdge2[j] <- paste(unlist(tip$leadingEdge[j]), collapse='/')
  }
  GSEA_list[[i]] <- tip[,-8]
}

#save source data
writexl::write_xlsx(GSEA_list, 'KEGG_source.xlsx')








#make a GSEA analysis per sample for either BTM or reactome
################################
#BP
##########################################
GSEA_list = list()
for(i in 1:ncol(FC)){
  sample <- data.frame('Fold' = FC[,i])
  sample$Gene <- rownames(FC)
  sample <- merge(sample, genes, by.x='Gene', by.y='GeneSymbol')
  sample <- sample[!duplicated(sample$Gene),]
  sample <- sample[order(sample$Fold, decreasing = T),]
  FC_List <- sample[,2]
  names(FC_List) <- sample[,1]
  fgseaRes <- fgseaMultilevel(BP, stats=FC_List, minSize = 10)
  GSEA_list[[i]] <- fgseaRes
  assign(paste0('GSEA_', i), GSEA_list[[i]])
}

#make the list into a dataframe
gseacolumns <- do.call(cbind, GSEA_list)
keep <- which(grepl('NES', names(gseacolumns)))

#extract the NES scores per sample and add the sample ID
NES <- gseacolumns[,c(1,6,14,22,30)]
colnames(NES) <- c('Pathway', as.character(post$ID))
NES.m <- melt(NES)
NES.m$Abs <- abs(NES.m$value)
ggplot(NES.m, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))

#extract the p-values per sample and add sample ID as column name
Padj <- gseacolumns[,c(1,3,11,19,27)] 
colnames(Padj) <- c('Pathway', as.character(post$ID))

#filter on significant pathways that are present in at least 2 per timepoint
#arrange in correct order
NES[Padj>0.05] <- 0 #set all the non-significant ones to 0
NES$Pathway <- Padj$Pathway #re-add the pathways

NES_s <- NES[rowSums(Padj < .05, na.rm = T)>0,] #filter on at least one significant 1


#melt and create absolute value for size of point
NES.m2 <- melt(NES_s)
NES.m2$Abs <- abs(NES.m2$value)

#make and save plot
ggplot(NES.m2, aes(x=variable, y=Pathway))+
  geom_point(aes(size = Abs, colour=value))+ 
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",high = "red")+
  xlab('Donor') + 
  theme_bw()
ggsave('ssGSEA_BP_2023_05_23.pdf', width=5, height=5)

for(i in 1:length(GSEA_list)){
  tip <- GSEA_list[[i]]
  tip$leadingEdge2 <- 'NA'
  for(j in 1:nrow(tip)){
    tip$leadingEdge2[j] <- paste(unlist(tip$leadingEdge[j]), collapse='/')
  }
  GSEA_list[[i]] <- tip[,-8]
}

#save source data
writexl::write_xlsx(GSEA_list, 'BP_source.xlsx')



