# Figure 4

library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(Seurat)
library(viridis)
library(ggpubr)


seu <- readRDS("TV_MEL006_GT_2510.rds")

str(seu@meta.data[ - grep("RNA_",colnames(seu@meta.data))])

cluster <- read.csv("clusterinfo.csv",header=T,row.names=1) # from CNV pipeline

rownames(cluster) <- sub("T00_", "", rownames(cluster))
rownames(cluster) <- sub("T04_", "", rownames(cluster))
rownames(cluster) <- sub("T28_", "", rownames(cluster))
rownames(cluster) <- sub("DNA_", "RNA_", rownames(cluster))

rownames(cluster) <- ifelse(grepl("RNA_S", rownames(cluster)), 
       sub("(.*)_(.*)", "\\1", rownames(cluster)), 
       rownames(cluster))


rownames(ABC@meta.data) <- ifelse(grepl("RNA_S", rownames(ABC@meta.data)), 
       sub("(.*)_(.*)", "\\1", rownames(ABC@meta.data)), 
       rownames(ABC@meta.data))

rownames(ABC@meta.data) <- sub("_T04", "", rownames(ABC@meta.data))
meta.data <- merge(ABC@meta.data, cluster, by ="row.names", all.x = T, all.y=F)
meta.data <- meta.data[ - grep("RNA_",colnames(meta.data))]
rownames(meta.data) <- colnames(ABC)

MEL006 <- read.csv("MEL006_sh.txt",sep="\t") # states
rownames(MEL006) <- MEL006$RNA_nr
MEL006 <- merge(meta.data, MEL006, by ="row.names", all= T)




Idents(ABC) <- 'cluster'


T0seurat <- subset(ABC, subset = timepoint == 'T0')
T4seurat <- subset(ABC, subset = timepoint == 'T04')
T28seurat <- subset(ABC, subset = timepoint == 'T28')


markersAC_T0 <- FindMarkers(T0seurat, ident.1 = "A", ident.2 = "C", min.pct = 0.3)
markersBC_T0 <- FindMarkers(T0seurat, ident.1 = "B", ident.2 = "C", min.pct = 0.3)

markersAC_T28 <- FindMarkers(T28seurat, ident.1 = "A", ident.2 = "C", min.pct = 0.3)
markersBC_T28 <- FindMarkers(T28seurat, ident.1 = "B", ident.2 = "C", min.pct = 0.3)

markersAC_T4 <- FindMarkers(T4seurat, ident.1 = "A", ident.2 = "C", min.pct = 0.3)
markersBC_T4 <- FindMarkers(T4seurat, ident.1 = "B", ident.2 = "C", min.pct = 0.3)

markersAB_T0 <- FindMarkers(T0seurat, ident.1 = "A", ident.2 = "B")
markersAB_T28 <- FindMarkers(T28seurat, ident.1 = "A", ident.2 = "B")
markersAB_T4 <- FindMarkers(T4seurat, ident.1 = "A", ident.2 = "B")


markersAC_T0$timepoint <- "T0"
markersAB_T0$timepoint <- "T0"
markersBC_T0$timepoint <- "T0"

markersAC_T0$comparison <- "A versus C"
markersAB_T0$comparison <- "A versus B"
markersBC_T0$comparison <- "B versus C"

markersAC_T0$genes <- rownames(markersAC_T0)
markersAB_T0$genes <- rownames(markersAB_T0)
markersBC_T0$genes <- rownames(markersBC_T0)

markersAC_T28$timepoint <- "T28"
markersBC_T28$timepoint <- "T28"
markersAB_T28$timepoint <- "T28"

markersAC_T28$comparison <- "A versus C"
markersBC_T28$comparison <- "B versus C"
markersAB_T28$comparison <- "A versus B"

markersAC_T28$genes <- rownames(markersAC_T28)
markersBC_T28$genes <- rownames(markersBC_T28)
markersAB_T28$genes <- rownames(markersAB_T28)

markersAB_T4$timepoint <- "T4"
markersAC_T4$timepoint <- "T4"
markersBC_T4$timepoint <- "T4"

markersAB_T4$comparison <- "A versus B"
markersAC_T4$comparison <- "A versus C"
markersBC_T4$comparison <- "B versus C"

markersAB_T4$genes <- rownames(markersAB_T4)
markersAC_T4$genes <- rownames(markersAC_T4)
markersBC_T4$genes <- rownames(markersBC_T4)



markers <- rbind(markersAB_T0,markersAC_T0,markersBC_T0,
                markersAB_T28,markersAC_T28,markersBC_T28,
                markersAB_T4,markersAC_T4,markersBC_T4)



markersDF <- markers[markers$p_val_adj < 0.05,]

genesAnno <- read.csv("genesAnnotation.csv",header = T)

df <- merge(markersDF,genesAnno,by="genes",all.x=F,all.y=F)

GOI_T0 <- df[df$timepoint=="T0",]
GOI_T28 <- df[df$timepoint=="T28",]

T0A_C <- GOI_T0[GOI_T0$comparison=="A versus C",]

T0A_C$Size <- pmax(T0A_C$pct.1, T0A_C$pct.2)
rownames(T0A_C) <- T0A_C$genes


# T0
T0A_B <-  GOI_T0[GOI_T0$comparison=="A versus B",]
T0A_B$Size <- pmax(T0A_B$pct.1, T0A_B$pct.2)
rownames(T0A_B) <- T0A_B$genes
#head(T0A_B) 

T0B_C <- GOI_T0[GOI_T0$comparison=="B versus C",]

T0B_C$Size <- pmax(T0B_C$pct.1, T0B_C$pct.2)
rownames(T0B_C) <- T0B_C$genes
#head(T0B_C)



# T28
T28A_C <- GOI_T28[GOI_T28$comparison=="A versus C",]

T28A_C$Size <- pmax(T28A_C$pct.1, T28A_C$pct.2)
rownames(T28A_C) <- T28A_C$genes
#head(T28A_C)


T28A_B <- GOI_T28[GOI_T28$comparison=="A versus B",]

T28A_B$Size <- pmax(T28A_B$pct.1, T28A_B$pct.2)
rownames(T28A_B) <- T28A_B$genes
#head(T28A_B) # empty

T28B_C <- GOI_T28[GOI_T28$comparison=="B versus C",]

T28B_C$Size <- pmax(T28B_C$pct.1, T28B_C$pct.2)
rownames(T28B_C) <- T28B_C$genes
#head(T28B_C)


fullT0 <- rbind(T0A_C,T0B_C,T0A_B)
fullT28 <- rbind(T28A_C,T28B_C)
full <- rbind(fullT0,fullT28)


rownames(fullT0) <- NULL
rownames(fullT28) <- NULL
rownames(full) <- NULL

full <- full[order(full$timepoint,decreasing = T),]


full$StateMarker <- as.factor(full$StateMarker)
full$StateMarker <- factor(full$StateMarker,levels=c("NCSC","SMC","Invasive","hyperdiff","Other"))




options(repr.plot.width=25, repr.plot.height=12)

dotchart <- ggdotchart(full, x = "genes", y = "avg_log2FC",
           color = "geneLocation",  # Color by groups
           shape="StateMarker",
           palette = c("deepskyblue","orange","darkblue","red", "grey"), # Custom color palette
           sorting = "none",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically                          
           dot.size = "Size",                                 # Large dot size
           #label = round(fullT28$avg_logFC,1),                        # Add mpg values as dot labels
           #font.label = list(color = "black", size = 10, 
           #                    vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr(),
           #repel = T,
           legend = "right"
           ) + xlab("") + ylab("Average fold change") +
        facet_grid( ~ timepoint + comparison , scales = "free_y",space='free') +
        #facet_grid( ~ comparison , scales = "free_y",space='free') +
        theme(panel.spacing = unit(2.5, "lines")) +
        theme(strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15),
                legend.title=element_text(size=18),
                legend.text=element_text(size=15),
                axis.text.y = element_text(size=10)) +
        geom_hline(yintercept=0, linetype="dashed", 
                color = "grey", size=1) +
        scale_shape_manual(values=c(15, 17, 6, 18, 19)) + 
        labs(color = "Gene location", size = "Percentage cells expressed",shape="Cell-state marker") +
          guides(color = guide_legend(order = 1, override.aes = list(size = 6)),
                 size = guide_legend(order = 3),
                 shape=guide_legend(order = 2))




# B - D

cnv <- read.csv("reanalyseCNV_10kb/CNseg.SPCF.35.txt",header=T,row.names=1)
cnv <- filter(cnv, chr %in% c(5,9,11,13,22))
cnv <- filter(cnv, (chr == 5 & start > 33000000 & end < 34000000) |
            (chr == 9 & start > 5800000 & end < 6000000) |
            (chr == 11 & start > 88800000 & end < 89100000) |
            (chr == 22 & start > 21000000 & end < 22000000) |
            (chr == 13 & start > 36000000 & end < 39000000) |
            (chr == 13 & start > 77000000 & end < 95200000))


segGenes <- t(cnv)
colnames(segGenes) <- paste0(colnames(segGenes),"_X")
segGenes <- tail(segGenes, -3)
segGenes <- as.data.frame(segGenes)

logcounts <- read.csv("lognormalisedcounts.csv",row.names=1)


genes <- logcounts[rownames(logcounts) %in% c("UFM1","CLN5","DCLK1","GPC5",
                                              "DCT","SLC45A2","TYR","MLANA",
                                              "CRKL","THAP7","LZTR1","PI4KA"
                                              ),]


rownames(genes) <- paste0(rownames(genes),"_Y")
rownames(segGenes) <- sub("T00_", "", rownames(segGenes))
rownames(segGenes) <- sub("T04_", "", rownames(segGenes))
rownames(segGenes) <- sub("T28_", "", rownames(segGenes))
rownames(segGenes) <- sub("DNA_", "RNA_", rownames(segGenes))
rownames(segGenes) <- ifelse(grepl("RNA_S", rownames(segGenes)), 
       sub("(.*)_(.*)", "\\1", rownames(segGenes)), 
       rownames(segGenes))
colnames(genes) <- sub("_T04", "", colnames(genes))
colnames(genes) <- ifelse(grepl("RNA_S", colnames(genes)), 
       sub("(.*)_(.*)", "\\1", colnames(genes)), 
       colnames(genes))


genes <- t(genes)


combined <- merge(segGenes,genes, by='row.names', all.y=T,all.x=F)

combined <- na.omit(combined)

combined <- combined %>% mutate_at(c("UFM1_X","CLN5_X","DCLK1_X",
                                     "GPC5_X","DCT_X","CRKL_X","THAP7_X",
                                     "LZTR1_X","PI4KA_X","SLC45A2_X","TYR_X"), 
                                   as.numeric)


cluster <- read.csv("clusterinfo.csv")
rownames(cluster) <- cluster$X
rownames(cluster) <- ifelse(grepl("DNA_S", rownames(cluster)), 
       sub("(.*)_(.*)", "\\1", rownames(cluster)), 
       rownames(cluster))
rownames(cluster) <- sub("T00_", "", rownames(cluster))
rownames(cluster) <- sub("T04_", "", rownames(cluster))
rownames(cluster) <- sub("T28_", "", rownames(cluster))
rownames(cluster) <- sub("DNA_", "RNA_", rownames(cluster))

combinedC <- merge(combined,cluster, by='row.names', all.x=T,all.y=F)
combinedC$X <- NULL
rownames(combinedC) <- combinedC$Row.names
combinedC$Row.names <- NULL

CLN <- combinedC[colnames(combinedC) %in% c("CLN5_X","CLN5_Y","cluster")]
CLN$gene <- "CLN5"
CLN$chr <- "chr13"
colnames(CLN) <- c("X","Y","Subclone","gene","chr")


DCLK1 <- combinedC[colnames(combinedC) %in% c("DCLK1_X","DCLK1_Y","cluster")]
DCLK1$gene <- "DCLK1"
DCLK1$chr <- "chr13"
colnames(DCLK1) <- c("X","Y","Subclone","gene","chr")

GPC5 <- combinedC[colnames(combinedC) %in% c("GPC5_X","GPC5_Y","cluster")]
GPC5$gene <- "GPC5"
GPC5$chr <- "chr13"
colnames(GPC5) <- c("X","Y","Subclone","gene","chr")

DCT <- combinedC[colnames(combinedC) %in% c("DCT_X","DCT_Y","cluster")]
DCT$gene <- "DCT"
DCT$chr <- "chr13_1"
colnames(DCT) <- c("X","Y","Subclone","gene","chr")

CRKL <- combinedC[colnames(combinedC) %in% c("CRKL_X","CRKL_Y","cluster")]
CRKL$gene <- "CRKL"
CRKL$chr <- "chr22"
colnames(CRKL) <- c("X","Y","Subclone","gene","chr")

THAP7 <- combinedC[colnames(combinedC) %in% c("THAP7_X","THAP7_Y","cluster")]
THAP7$gene <- "THAP7"
THAP7$chr <- "chr22"
colnames(THAP7) <- c("X","Y","Subclone","gene","chr")

LZTR1 <- combinedC[colnames(combinedC) %in% c("LZTR1_X","LZTR1_Y","cluster")]
LZTR1$gene <- "LZTR1"
LZTR1$chr <- "chr22"
colnames(LZTR1) <- c("X","Y","Subclone","gene","chr")

PI4KA <- combinedC[colnames(combinedC) %in% c("PI4KA_X","PI4KA_Y","cluster")]
PI4KA$gene <- "PI4KA"
PI4KA$chr <- "chr22"
colnames(PI4KA) <- c("X","Y","Subclone","gene","chr")

SLC45A2 <- combinedC[colnames(combinedC) %in% c("SLC45A2_X","SLC45A2_Y","cluster")]
SLC45A2$gene <- "SLC45A2"
SLC45A2$chr <- "chr5"
colnames(SLC45A2) <- c("X","Y","Subclone","gene","chr")

TYR <- combinedC[colnames(combinedC) %in% c("TYR_X","TYR_Y","cluster")]
TYR$gene <- "TYR"
TYR$chr <- "chr11"
colnames(TYR) <- c("X","Y","Subclone","gene","chr")

UFM1 <- combinedC[colnames(combinedC) %in% c("UFM1_X","UFM1_Y","cluster")]
UFM1$gene <- "UFM1"
UFM1$chr <- "chr13"
colnames(UFM1) <- c("X","Y","Subclone","gene","chr")

MLANA <- combinedC[colnames(combinedC) %in% c("MLANA_X","MLANA_Y","cluster")]
MLANA$gene <- "MLANA"
MLANA$chr <- "chr9"
colnames(MLANA) <- c("X","Y","Subclone","gene","chr")


df <- rbind(UFM1,CLN,DCLK1,GPC5,DCT,CRKL,THAP7,LZTR1,PI4KA,SLC45A2,TYR,MLANA)


plot <-     plot_grid(ggplot(df[df$chr=="chr13",],aes(x=log(X),y=Y,group=Subclone,color = Subclone))+ 
        geom_point(size=1) + 
        geom_smooth(method = "lm", formula = y ~ x,aes(fill = Subclone)) +
        stat_cor() +
        ylim(0,6) + 
        xlim(0,5) + 
        scale_color_manual(values=c("#00A087FF","#3C5488FF","#F39B7FFF")) +
        scale_fill_manual(values=c("#00A087FF","#3C5488FF","#F39B7FFF")) +
        theme_bw() +
        ylab("log(normalized counts)") +
        xlab("log(copy number)") +
        #labs(col="subclones") + 
        facet_wrap(vars(gene), ncol = 4,scales='free') +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()),
         
        ggplot(df[df$chr=="chr22",],aes(x=log(X),y=Y,group=Subclone,color = Subclone))+ 
            geom_point(size=1) + 
            geom_smooth(method = "lm", formula = y ~ x,aes(fill = Subclone)) +
            stat_cor() +
            ylim(0,6) + 
            xlim(0,5) + 
            scale_color_manual(values=c("#00A087FF","#3C5488FF","#F39B7FFF")) +
            scale_fill_manual(values=c("#00A087FF","#3C5488FF","#F39B7FFF")) +
            theme_bw() +
            ylab("log(normalized counts)") +
            xlab("log(copy number)") +
            #labs(col="subclones") + 
            facet_wrap(vars(gene), ncol = 4,scales='free') +
            theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())  ,

        ggplot(df[df$chr!="chr13" & df$chr!="chr22",],aes(x=log(X),y=Y,group=Subclone,color = Subclone))+ 
            geom_point(size=1) + 
            geom_smooth(method = "lm", formula = y ~ x,aes(fill = Subclone)) +
            stat_cor() +
            ylim(0,6) + 
            xlim(0,5) + 
            scale_color_manual(values=c("#00A087FF","#3C5488FF","#F39B7FFF")) +
            scale_fill_manual(values=c("#00A087FF","#3C5488FF","#F39B7FFF")) +
            theme_bw() +
            ylab("log(normalized counts)") +
            xlab("log(copy number)") +
            #labs(col="subclones") + 
            facet_wrap(vars(gene), ncol = 4,scales='free') +
            theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()),

        nrow=3,
        labels=c("B","C","D"),label_size = 18) 

