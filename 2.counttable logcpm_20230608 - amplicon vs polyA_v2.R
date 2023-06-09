#Clears the workspace
rm(list = ls())

setwd('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\RNA-Seq\\Liver samples Jutte\\20230608_compare_amplicon_polyA')

#loads the required libraries, but suppresses the output to the screen 
suppressMessages(lapply(c('edgeR','MASS' ,'stringr', 'ggplot2', 'pheatmap', 'writexl',
                          'cowplot', 'reshape2'), require, character.only = TRUE))

#writes the session information to an output file
sink("sessionInfo.txt")
sessionInfo()
sink()


#laod in the data
poly <- read.delim('../20221201_newfilter/logcpm_selected.txt')
ampli <- read.delim('../20230608_amplicon_newaligned/output/logCPM_amplicon_sel.txt')

#load in metadata
meta <- read.delim('../input/metadata_ampliseq_v2.txt')


#fix the column names to match the samples, removing the unmatched samples
ampli <- ampli[,-which(colnames(ampli) == 'Co07')]
poly <- poly[,-which(colnames(poly) == 'X003R')]

#make the column names the same
colnames(poly) <- colnames(ampli) 


#select the genes in both
poly2 <- poly[rownames(poly) %in% rownames(ampli),]
ampli2 <- ampli[rownames(ampli) %in% rownames(poly),]

#put in same order and check
poly3 <- poly2[order(rownames(poly2)),]
ampli3 <- ampli2[order(rownames(ampli2)),]
identical(rownames(poly3), rownames(ampli3))


#stick both in ggplot
poly3$Gene <- rownames(poly3)
ampli3$Gene <- rownames(ampli3)
colnames(ampli3) <- colnames(poly3)


poly4 <- melt(poly3, id = c('Gene'))
ampli4 <- melt(ampli3, id = c('Gene'))


all <- cbind(poly4, ampli4[,3])
colnames(all) <- c('Gene', 'Patient', 'PolyA', 'AmpliSeq')

ggplot(all, aes(x=PolyA, y=AmpliSeq))  + geom_point()

cor.test(all$PolyA, all$AmpliSeq)


all2 <- all[all$PolyA > min(all$PolyA) & 
              all$AmpliSeq > min(all$AmpliSeq),  ]

ggplot(all2, aes(x=PolyA, y=AmpliSeq))  + geom_point() + 
  stat_smooth(method='lm') + theme_bw() +  
  ggtitle('correlation: rho = 0.54, p = <2.2x10^-16') + 
  theme(aspect.ratio = 1)
ggsave('correlations_normalized_detected.pdf', width = 4, height = 4)
cor.test(all2$PolyA, all2$AmpliSeq)




sigs <- readxl::read_xlsx('../20230608_amplicon_newaligned/output/comparisons_volcano.xlsx', 1)
signame <- as.character(sigs$name[sigs$Log2FC_case_ctrl != 0]) #choose either upregulated or all
signame <- as.character(sigs$name[sigs$Log2FC_case_ctrl > 0])
signame2 <- signame[!is.na(signame)]


table(rownames(poly3) %in% signame2)



all3 <- all
all3$Sig <- 'No'
all3$Sig[all3$Gene %in% signame2] <- 'Yes'
all3$Type <- 'Control'
all3$Type[all3$Patient %in% c('Ca05', 'Ca01', 'Ca02', 'Ca03' )] <-'Patient'


ggplot(all3, aes(x=PolyA, y=AmpliSeq))  + geom_point() + 
  stat_smooth(method='lm') + theme_bw() +  
  ggtitle('correlation: rho = 0.45, p = <2.2x10^-15') + 
  theme(aspect.ratio = 1) + facet_grid(.~Sig)


all4 <- all3[all3$PolyA > min(all3$PolyA) & 
              all3$AmpliSeq > min(all3$AmpliSeq),  ]

ggplot(all4, aes(x=PolyA, y=AmpliSeq))  + geom_point(aes(colour=Type)) + 
  stat_smooth(method='lm') + theme_bw() +  
  theme(aspect.ratio = 1) + facet_grid(.~Sig)
ggsave('correlations_normalized_detected_sigbased.pdf', width = 8, height = 4)


all5 <- all3[all3$Sig == 'Yes',]
all5m <- melt(all5, id = colnames(all5)[c(1,2,5,6)])
ggplot(all5m, aes(x=Type, y=value)) + 
  geom_boxplot(outlier.colour = NA) + facet_grid(Gene ~ variable) + 
  geom_jitter(aes(colour = Patient), width = 0.1, height = 0)
ggsave('comparison_sigs_2.pdf', width = 4, height = 40)
