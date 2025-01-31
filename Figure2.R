
# Figure 2


library(ComplexHeatmap)
library(factoextra)
library(RColorBrewer)
library(Hmisc)
library(matrixStats)
library(circlize)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(DescTools)
library(tidyr)




## Heatmap (A)

selected_path <- "cnv_results/Matrices/selected.txt" # CNV pipeline
segCN_path <- "cnv_results/Matrices/CNseg.SPCF.10_500K.txt" # CNV pipeline
locationcluster <- "cnv_results/clusterinfo.csv"
colscheme_path <- "colscheme"
locationSampleDescr <- "samplesheetMEL006.csv"
species <- "Human"



selected <- read.csv(selected_path, header = T, row.names = 1)
ploidy <- round(median(selected$X))

segCN <- read.csv(segCN_path, header=T, row.names=1, check.names = T)
segCN$chr<-factor(segCN$chr,levels = c(as.character(1:22),"X",'Y'))


sampleinfo <- read.csv(locationSampleDescr, header=T, sep="\t")
rownames(sampleinfo) <- sampleinfo$sample
sampleinfo <- sampleinfo[rownames(sampleinfo) %in% colnames(segCN),]

colscheme <- read.csv(colscheme_path)
colnames(colscheme) <- c("timepoint", "method", "Empty")

cluster <- read.csv(locationcluster, header=T, row.names=1)
cluster$name <- rownames(cluster)

sampleinfo <- merge(sampleinfo, cluster, by=0)
colnames(sampleinfo)[ncol(sampleinfo)]<-"cluster"
rownames(sampleinfo)<-sampleinfo$Row.names
sampleinfo <- sampleinfo[rownames(sampleinfo)[order(rownames(sampleinfo))],]

int1 <- colnames(colscheme)[1]
int2 <- colnames(colscheme)[2]
int3 <- colnames(colscheme)[3]
collist1 <- colscheme[,1]
collist2 <- colscheme[,2]
collist3 <- colscheme[,3]

    
if(length(rownames(selected))>length(rownames(sampleinfo))){
  selected <- selected[rownames(selected) %in% rownames(sampleinfo),]
  segCN <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(sampleinfo))]
  
}else {
  sampleinfo <- sampleinfo[rownames(sampleinfo) %in% rownames(selected),]
  segCN <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(selected))]
  
}

sampleinfo <- sampleinfo[order(rownames(sampleinfo)),]
selected <- selected[order(rownames(selected)),]





dm <- segCN
interest1 <- int1
interest2 <- int2
interest3 <- int3

#chrs will be used as splits for the heatmap
sections <- c()
for(i in unique(dm$chr)){
    sections <- c(sections, rep(toString(i), sum(dm$chr==i)))
}
split <- data.frame(sections)

CNheatmapmatrix<-dm[,-c(1,2,3)]
CNheatmapmatrix[CNheatmapmatrix<0] <- 0
CNheatmapmatrix <- CNheatmapmatrix[,colnames(CNheatmapmatrix)[order(colnames(CNheatmapmatrix))]]
Sds <- colSds(as.matrix(CNheatmapmatrix))
ids <- which(Sds==0)

first <- ComplexHeatmap::draw(Heatmap(t(as.matrix(CNheatmapmatrix)),column_split = factor(split$sections, levels=unique(dm$chr)), 
                                          cluster_columns = F, show_column_names = F, 
                                          #cluster_rows = F,
                                          show_row_names = F, show_heatmap_legend = T, cluster_column_slices = FALSE,
                                          column_title=unique(dm$chr), column_title_gp = gpar(fontsize=10),
                                          clustering_distance_rows = "canberra", clustering_method_rows = "ward.D2"))
order<-row_order(first)

  #' Next, we use the estimated ploidy to give a color code to the heatmap.
  #' As the range of the CNs is usually very broad, we cap the CN range, so that the
  #' legend of the heatmap does not become too crowded. However, since this changes
  #' the original data, we will need the original clustering from before for
  #' downstream plotting
  if(ploidy==4){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "royalblue2","#0571b0", "#92c5de", "grey96", "lightsalmon", "tomato", "tomato3", "red4", "maroon", "magenta4", "plum3","purple4","black"), names = c(0:(3*ploidy),20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  }else if(ploidy==3){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "#0571b0", "#92c5de", "grey96", "lightsalmon", "tomato", "tomato3", "red4", "maroon", "magenta4", "black"), names = c(0:9,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  } else if(ploidy==2){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "#0571b0", "grey96", "lightsalmon","tomato", "red3", "red4", "black"), names = c(0:6,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  } else if(ploidy==1){
    CNheatmapmatrix[CNheatmapmatrix>5] <- 20
    colorshtmp = structure(c("navyblue", "grey96", "lightsalmon","tomato", "red3", "red4", "black"), names = c(0:5,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:5,20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  }
      
#' the column annotation for the heatmap, each chromosome will get a different color
#' alternating grey and black
col <- list(Chr=rep(c("grey","black"),round(length(unique(dm$chr))/2)))
col$Chr <- col$Chr[1:length(unique(dm$chr))]
names(col$Chr) <- unique(dm$chr)
column_ha = HeatmapAnnotation(Chr=factor(split$sections, levels=unique(dm$chr)), show_legend = F, show_annotation_name = F, col = col)

  #' now for the row annotation, depending on if you want labels for 
  #' one or for more interests, which are specified in the sample description.
  #' If there was no sample description, all flags have been put to empty, and nothing
  #' will happen. If there was, the colorscheme will be applied for each interest
  col_fun2 = colorRamp2(c(min(selected[,1]),ploidy,max(selected[,1])), c("#0571b0","grey96", "tomato"))
  col_fun2(seq(min(selected[,1]), max(selected[,1])))


colclust <- brewer.pal(n=length(unique(sampleinfo$cluster)), name="Set3")
names(colclust) <- as.character(unique(sampleinfo["cluster"][,1]))
colclust <- list(clus=colclust)

row_ha_right=HeatmapAnnotation(#bar=anno_text(selected[,1], location = 0.5, just = "center",
                                #              gp = gpar(fill = col_fun2(selected[,1]), col = "black", border = "black"),
                                #              width = max_text_width(selected[,1])*1.2), 
                                clus=sampleinfo["cluster"][,1],
                                which="row",
                                #col=list(bar=col_fun2, colclust), 
                                col=colclust,
                                show_annotation_name = F, 
                                annotation_legend_param = list(clus=list(title="Cluster")),
                                annotation_name_rot = 90)

  if (interest1 == "Empty") {
    row_ha=NULL
  } else if (interest2 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfo[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfo[interest1][,1])[order(unique(sampleinfo[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    row_ha = rowAnnotation(int1 = sampleinfo[interest1][, 1],show_legend = T,show_annotation_name = F,
                           col = uniqint1, annotation_legend_param = list(int1 = list(title = interest1)),
                           annotation_name_rot = 90)
  } else if (interest3 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfo[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfo[interest1][,1])[order(unique(sampleinfo[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfo[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfo[interest2][,1])[order(unique(sampleinfo[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    row_ha = rowAnnotation(int1 = sampleinfo[interest1][, 1],int2 = sampleinfo[interest2][, 1],
                           show_legend = T,show_annotation_name = F,col = c(uniqint1, uniqint2), 
                           annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2))),
                           annotation_name_rot = 90)
  } else {
    uniqint1 <- collist1[1:length(unique(sampleinfo[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfo[interest1][,1])[order(unique(sampleinfo[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfo[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfo[interest2][,1])[order(unique(sampleinfo[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    uniqint3 <- collist3[1:length(unique(sampleinfo[interest3][,1]))]
    names(uniqint3) <- as.character(unique(sampleinfo[interest3][,1])[order(unique(sampleinfo[interest3][,1]))])
    uniqint3 <- list(int3=uniqint3)
    
    row_ha = rowAnnotation(
      int1 = sampleinfo[interest1][, 1],int2 = sampleinfo[interest2][, 1],int3 = sampleinfo[interest3][, 1],
      show_legend = T,show_annotation_name = F, col = c(uniqint1, uniqint2, uniqint3), 
      annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2)),int3 = list(title = paste0("\n", interest3))),
      annotation_name_rot = 90)
 }


ht <- Heatmap(t(as.matrix(CNheatmapmatrix)), heatmap_legend_param = lgd, 
                row_order = order,
                column_split = factor(split$sections, levels=unique(dm$chr)), cluster_columns = F, show_column_names = F, 
                show_row_names = F, show_heatmap_legend = T, cluster_column_slices = FALSE,
                column_title=unique(dm$chr), column_title_gp = gpar(fontsize=10),top_annotation = column_ha, 
                left_annotation = row_ha, 
                right_annotation = row_ha_right,
                split=cluster$cluster,
                col=colorshtmp, clustering_distance_rows = "canberra", clustering_method_rows = "ward.D2")




## Heatmap - pseudobulk (B)

# to make random list of samples for each group:
# Gtag
A_0 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="A" & sampleinfo$timepoint=="T0",]
A_4 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="A" & sampleinfo$timepoint=="T04",]
A_28 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="A" & sampleinfo$timepoint=="T28",]

A_0 <- sample_n(A_0,14)
A_4 <- sample_n(A_4,14)
A_28 <- sample_n(A_28,14)

B_0 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="B" & sampleinfo$timepoint=="T0",]
B_4 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="B" & sampleinfo$timepoint=="T04",]
B_28 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="B" & sampleinfo$timepoint=="T28",]

B_0 <- sample_n(B_0,14)
B_4 <- sample_n(B_4,14)
B_28 <- sample_n(B_28,14)

C_0 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="C" & sampleinfo$timepoint=="T0",]
C_4 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="C" & sampleinfo$timepoint=="T04",]
C_28 <- sampleinfo[sampleinfo$method=="Gtag" & sampleinfo$cluster=="C" & sampleinfo$timepoint=="T28",]

C_0 <- sample_n(C_0,14)
C_4 <- sample_n(C_4,14)
C_28 <- sample_n(C_28,14)

# Picoplex
A_0 <- sampleinfo[sampleinfo$method=="picoPlex" & sampleinfo$cluster=="A" & sampleinfo$timepoint=="T0",]
A_28 <- sampleinfo[sampleinfo$method=="picoPlex" & sampleinfo$cluster=="A" & sampleinfo$timepoint=="T28",]

A_0 <- sample_n(A_0,14)
A_28 <- sample_n(A_28,14)

B_0 <- sampleinfo[sampleinfo$method=="picoPlex" & sampleinfo$cluster=="B" & sampleinfo$timepoint=="T0",]
B_28 <- sampleinfo[sampleinfo$method=="picoPlex" & sampleinfo$cluster=="B" & sampleinfo$timepoint=="T28",]

B_0 <- sample_n(B_0,14)
B_28 <- sample_n(B_28,14)

C_0 <- sampleinfo[sampleinfo$method=="picoPlex" & sampleinfo$cluster=="C" & sampleinfo$timepoint=="T0",]
C_28 <- sampleinfo[sampleinfo$method=="picoPlex" & sampleinfo$cluster=="C" & sampleinfo$timepoint=="T28",]

C_0 <- sample_n(C_0,14)
C_28 <- sample_n(C_28,14)




# heatmap pseudobulks
selected_path <- "cnv_resultsPseudobulks/Matrices/selected.txt"
segCN_path <- "cnv_resultsPseudobulks/Matrices/CNseg.SPCF.10.txt"
colscheme <- "colscheme"
locationSampleDescr <- "pseudobulks/samplesheet.csv"
species <- "Human"

selected <- read.csv(selected_path, header = T, row.names = 1)
ploidy <- round(median(selected$X))

segCN <- read.csv(segCN_path, header=T, row.names=1, check.names = T)
segCN$chr<-factor(segCN$chr,levels = c(as.character(1:22),"X",'Y'))

sampleinfo <- read.csv(locationSampleDescr, header=T, sep="\t")

colscheme <- read.csv(colscheme)
collist1 <- colscheme[,1]
collist2 <- colscheme[,2]
collist3 <- colscheme[,2]



sampleinfo <- read.csv(locationSampleDescr, header=T, sep="\t")
sampleinfo$sample <- rownames(sampleinfo)

sampleinfoG <- sampleinfo[sampleinfo$method=="Gtag",]
sampleinfoP <- sampleinfo[sampleinfo$method=="picoplex",]


if(length(rownames(selected))>length(rownames(sampleinfoG))){
  selectedG <- selected[rownames(selected) %in% rownames(sampleinfoG),]
  segCNG <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(sampleinfoG))]
  
}else {
  sampleinfoG <- sampleinfo[rownames(sampleinfoG) %in% rownames(selected),]
  segCNG <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(selected))]
  
}

sampleinfoG <- sampleinfoG[order(rownames(sampleinfoG)),]
selectedG <- selectedG[order(rownames(selectedG)),]

dm <- segCNG
interest1 <- "method"
interest2 <- "timepoint"
interest3 <- "subclone"


  #chrs will be used as splits for the heatmap
  sections <- c()
  for(i in unique(dm$chr)){
    sections <- c(sections, rep(toString(i), sum(dm$chr==i)))
  }
  split <- data.frame(sections)
  
  CNheatmapmatrix<-dm[,-c(1,2,3)]
  CNheatmapmatrix[CNheatmapmatrix<0] <- 0
  CNheatmapmatrix <- CNheatmapmatrix[,colnames(CNheatmapmatrix)[order(colnames(CNheatmapmatrix))]]
  Sds <- colSds(as.matrix(CNheatmapmatrix))
  ids <- which(Sds==0)


colnames(CNheatmapmatrix) <- c("A T0","B T0","C T0",
                           "A T28","B T28","C T28",
                           "A T4","B T4","C T4")



    first <- ComplexHeatmap::draw(Heatmap(t(as.matrix(CNheatmapmatrix)),column_split = factor(split$sections, levels=unique(dm$chr)), 
                                          cluster_columns = F, show_column_names = F, 
                                          show_row_names = F, show_heatmap_legend = T, cluster_column_slices = FALSE,
                                          column_title=unique(dm$chr), column_title_gp = gpar(fontsize=10),
                                          clustering_distance_rows = "canberra", clustering_method_rows = "ward.D2"
                                          #cluster_rows = FALSE
                                         ))
    order<-row_order(first)

  if(ploidy==4){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "royalblue2","#0571b0", "#92c5de", "grey96", "lightsalmon", "tomato", "tomato3", "red4", "maroon", "magenta4", "plum3","purple4","black"), names = c(0:(3*ploidy),20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  }else if(ploidy==3){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "#0571b0", "#92c5de", "grey96", "lightsalmon", "tomato", "tomato3", "red4", "maroon", "magenta4", "black"), names = c(0:9,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  } else if(ploidy==2){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "#0571b0", "grey96", "lightsalmon","tomato", "red3", "red4", "black"), names = c(0:6,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  } else if(ploidy==1){
    CNheatmapmatrix[CNheatmapmatrix>5] <- 20
    colorshtmp = structure(c("navyblue", "grey96", "lightsalmon","tomato", "red3", "red4", "black"), names = c(0:5,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:5,20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  }
  #' the column annotation for the heatmap, each chromosome will get a different color
  #' alternating grey and black
  col <- list(Chr=rep(c("grey","black"),round(length(unique(dm$chr))/2)))
  col$Chr <- col$Chr[1:length(unique(dm$chr))]
  names(col$Chr) <- unique(dm$chr)
  column_ha = HeatmapAnnotation(Chr=factor(split$sections, levels=unique(dm$chr)), show_legend = F, show_annotation_name = F, col = col)
  
  #' now for the row annotation, depending on if you want labels for 
  #' one or for more interests, which are specified in the sample description.
  #' If there was no sample description, all flags have been put to empty, and nothing
  #' will happen. If there was, the colorscheme will be applied for each interest
  col_fun2 = colorRamp2(c(min(selectedG[,1]),ploidy,max(selectedG[,1])), c("#0571b0","grey96", "tomato"))
  col_fun2(seq(min(selectedG[,1]), max(selectedG[,1])))
  
  row_ha_right=HeatmapAnnotation(bar=anno_text(selectedG[,1], location = 0.5, just = "center",
                                               gp = gpar(fill = col_fun2(selectedG[,1]), col = "black", border = "black"),
                                               width = max_text_width(selectedG[,1])*1.2), 
                                  which="row", show_annotation_name = T, show_legend = T,
                                 annotation_legend_param = list(title="\nSelected"),
                                 annotation_name_rot = 90)
  
  if (interest1 == "Empty") {
    row_ha=NULL
  } else if (interest2 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfoG[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfoG[interest1][,1])[order(unique(sampleinfoG[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    row_ha = rowAnnotation(int1 = sampleinfoG[interest1][, 1],show_legend = F,show_annotation_name = F,
                           col = uniqint1, annotation_legend_param = list(int1 = list(title = interest1)),
                           annotation_name_rot = 90)
  } else if (interest3 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfoG[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfoG[interest1][,1])[order(unique(sampleinfoG[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfoG[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfoG[interest2][,1])[order(unique(sampleinfoG[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    row_ha = rowAnnotation(int1 = sampleinfoG[interest1][, 1],int2 = sampleinfoG[interest2][, 1],
                           show_legend = F,show_annotation_name = F,col = c(uniqint1, uniqint2), 
                           annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2))),
                           annotation_name_rot = 90)
  } else {
    uniqint1 <- collist1[1:length(unique(sampleinfoG[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfoG[interest1][,1])[order(unique(sampleinfoG[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfoG[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfoG[interest2][,1])[order(unique(sampleinfoG[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    uniqint3 <- collist3[1:length(unique(sampleinfoG[interest3][,1]))]
    names(uniqint3) <- as.character(unique(sampleinfoG[interest3][,1])[order(unique(sampleinfoG[interest3][,1]))])
    uniqint3 <- list(int3=uniqint3)
    
    row_ha = rowAnnotation(
      int1 = sampleinfoG[interest1][, 1],int2 = sampleinfoG[interest2][, 1],int3 = sampleinfoG[interest3][, 1],
      show_legend = F,show_annotation_name = F, col = c(uniqint1, uniqint2, uniqint3), 
      annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2)),int3 = list(title = paste0("\n", interest3))),
      annotation_name_rot = 90)
  }
  

htG <- Heatmap(t(as.matrix(CNheatmapmatrix)), heatmap_legend_param = lgd, 
                row_order = order,
                column_split = factor(split$sections, levels=unique(dm$chr)), cluster_columns = F, show_column_names = F, 
                show_row_names = T, show_heatmap_legend = F, cluster_column_slices = FALSE,
                column_title=unique(dm$chr), column_title_gp = gpar(fontsize=10),top_annotation = column_ha, 
                #left_annotation = row_ha, 
                #right_annotation = row_ha_right,
                split=sampleinfoG$method,
                col=colorshtmp, 
               clustering_distance_rows = "canberra", clustering_method_rows = "ward.D2"
)



# picoplex

if(length(rownames(selected))>length(rownames(sampleinfoP))){
  selectedP <- selected[rownames(selected) %in% rownames(sampleinfoP),]
  segCNP <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(sampleinfoP))]
  
}else {
  sampleinfoP <- sampleinfo[rownames(sampleinfo) %in% rownames(selected),]
  segCNP <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(selected))]
  
}

sampleinfoP <- sampleinfoP[order(rownames(sampleinfoP)),]
selectedP <- selectedP[order(rownames(selectedP)),]

dm <- segCNP

  #chrs will be used as splits for the heatmap
  sections <- c()
  for(i in unique(dm$chr)){
    sections <- c(sections, rep(toString(i), sum(dm$chr==i)))
  }
  split <- data.frame(sections)
  
  CNheatmapmatrix<-dm[,-c(1,2,3)]
  CNheatmapmatrix[CNheatmapmatrix<0] <- 0
  CNheatmapmatrix <- CNheatmapmatrix[,colnames(CNheatmapmatrix)[order(colnames(CNheatmapmatrix))]]
  Sds <- colSds(as.matrix(CNheatmapmatrix))
  ids <- which(Sds==0)

colnames(CNheatmapmatrix) <- c("A T0","B T0","C T0",
                           "A T28","B T28","C T28")





    first <- ComplexHeatmap::draw(Heatmap(t(as.matrix(CNheatmapmatrix)),column_split = factor(split$sections, levels=unique(dm$chr)), 
                                          cluster_columns = F, show_column_names = F, 
                                          show_row_names = F, show_heatmap_legend = T, cluster_column_slices = FALSE,
                                          column_title=unique(dm$chr), column_title_gp = gpar(fontsize=10),
                                          clustering_distance_rows = "canberra", clustering_method_rows = "ward.D2"
                                          #cluster_rows = FALSE
                                         ))
    order<-row_order(first)


  if(ploidy==4){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "royalblue2","#0571b0", "#92c5de", "grey96", "lightsalmon", "tomato", "tomato3", "red4", "maroon", "magenta4", "plum3","purple4","black"), names = c(0:(3*ploidy),20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  }else if(ploidy==3){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "#0571b0", "#92c5de", "grey96", "lightsalmon", "tomato", "tomato3", "red4", "maroon", "magenta4", "black"), names = c(0:9,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  } else if(ploidy==2){
    CNheatmapmatrix[CNheatmapmatrix>3*ploidy] <- 20
    colorshtmp = structure(c("navyblue", "#0571b0", "grey96", "lightsalmon","tomato", "red3", "red4", "black"), names = c(0:6,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:(3*ploidy),20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  } else if(ploidy==1){
    CNheatmapmatrix[CNheatmapmatrix>5] <- 20
    colorshtmp = structure(c("navyblue", "grey96", "lightsalmon","tomato", "red3", "red4", "black"), names = c(0:5,20))
    lngth <- length(colorshtmp)
    ats <- data.frame(value=c(0:5,20), index=1:length(colorshtmp))
    ats <- ats[ats$value %in% unique(unlist(CNheatmapmatrix)),]
    colorshtmp <- colorshtmp[ats$index]
    if(tail(ats$index,1)==lngth){
      labels <- c(as.character(ats$value)[1:(length(ats$value)-1)], paste0(">=", (3*ploidy)+1))
    } else {
      labels <- c(as.character(ats$value))
    }
    lgd <- list(col_fun = colorshtmp, title = "\nCopy number", at = ats$value, 
                labels = labels)
  }
  #' the column annotation for the heatmap, each chromosome will get a different color
  #' alternating grey and black
  col <- list(Chr=rep(c("grey","black"),round(length(unique(dm$chr))/2)))
  col$Chr <- col$Chr[1:length(unique(dm$chr))]
  names(col$Chr) <- unique(dm$chr)
  column_ha = HeatmapAnnotation(Chr=factor(split$sections, levels=unique(dm$chr)), show_legend = F, show_annotation_name = F, col = col)
  
  #' now for the row annotation, depending on if you want labels for 
  #' one or for more interests, which are specified in the sample description.
  #' If there was no sample description, all flags have been put to empty, and nothing
  #' will happen. If there was, the colorscheme will be applied for each interest
  col_fun2 = colorRamp2(c(min(selectedP[,1]),ploidy,max(selectedP[,1])), c("#0571b0","grey96", "tomato"))
  col_fun2(seq(min(selectedP[,1]), max(selectedP[,1])))
  
  row_ha_right=HeatmapAnnotation(bar=anno_text(selectedP[,1], location = 0.5, just = "center",
                                               gp = gpar(fill = col_fun2(selectedP[,1]), col = "black", border = "black"),
                                               width = max_text_width(selectedP[,1])*1.2), 
                                  which="row", show_annotation_name = T, show_legend = T,
                                 annotation_legend_param = list(title="\nSelected"),
                                 annotation_name_rot = 90)
  
  if (interest1 == "Empty") {
    row_ha=NULL
  } else if (interest2 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfoP[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfoP[interest1][,1])[order(unique(sampleinfoP[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    row_ha = rowAnnotation(int1 = sampleinfoP[interest1][, 1],
                           show_legend = F,
                           show_annotation_name = F,
                           col = uniqint1, annotation_legend_param = list(int1 = list(title = interest1)),
                           annotation_name_rot = 90)
  } else if (interest3 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfoP[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfoP[interest1][,1])[order(unique(sampleinfoP[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[2]
    names(uniqint2) <- as.character(unique(sampleinfoP[interest2][,1])[order(unique(sampleinfoP[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    row_ha = rowAnnotation(int1 = sampleinfoP[interest1][, 1],int2 = sampleinfoP[interest2][, 1],
                           show_legend = F,
                           show_annotation_name = F,col = c(uniqint1, uniqint2), 
                           annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2))),
                           annotation_name_rot = 90)
  } else {
    uniqint1 <- collist1[1:length(unique(sampleinfoP[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfoP[interest1][,1])[order(unique(sampleinfoP[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfoP[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfoP[interest2][,1])[order(unique(sampleinfoP[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    uniqint3 <- collist3[1:length(unique(sampleinfoP[interest3][,1]))]
    names(uniqint3) <- as.character(unique(sampleinfoP[interest3][,1])[order(unique(sampleinfoP[interest3][,1]))])
    uniqint3 <- list(int3=uniqint3)
    
    row_ha = rowAnnotation(
      int1 = sampleinfoP[interest1][, 1],int2 = sampleinfoP[interest2][, 1],int3 = sampleinfoP[interest3][, 1],
      show_legend = F,
      show_annotation_name = F, 
      col = c(uniqint1, uniqint2, uniqint3), 
      annotation_legend_param = list(int1 = list(title = interest1), 
                                     int2 = list(title = paste0("\n", interest2)),
                                     int3 = list(title = paste0("\n", interest3))),
      annotation_name_rot = 90)
  }
  


htP <- Heatmap(t(as.matrix(CNheatmapmatrix)), 
                heatmap_legend_param = lgd, 
                row_order = order,
                column_split = factor(split$sections, 
                levels=unique(dm$chr)), cluster_columns = F, show_column_names = F, 
                show_row_names = T, 
                show_heatmap_legend = F, cluster_column_slices = F,
                column_title=unique(dm$chr), 
                column_title_gp = gpar(fontsize=10),
                top_annotation = column_ha, 
                #left_annotation = row_ha, 
                #right_annotation = row_ha_right,
                split=sampleinfoP$method,
                col=colorshtmp, 
               clustering_distance_rows = "canberra", clustering_method_rows = "ward.D2"
              )





# bin plot - pseudobulks (C)


name <- "T28_C_Gtag"

seglogR <- seglogR[,c("chr","start","end",name)]
logR <- logR[,c("chr","start","end",name)]


pls <- seq(1.2,6,0.01)

#' set up two matrices, 'selected' to store the optimal result for each sample
#' (i.e. the global minimum of the grid)
selected <- matrix(ncol = 1, nrow=1)


scores <- data.frame(ploidygrid=pls, score=NA)

#' find the breakpoints and thus the segments for each sample

brkpts<-findbrkpts(seglogR[,c(1:4)])[[1]]
pens<-c()
for(pl in pls){
    #' for each possible average CN in our 1.2-6 grid,
    #' we transform the segmented logR to CN values, once rounded
    #' and once unrounded.
    CN<-pmax((2^seglogR[,4])*pl,0)
    segCN<-pmax(round((2^seglogR[,4])*pl),0)
    sum<-0
    for(brkpt in 1:nrow(brkpts)){
        #' Then, for each segment, we calculate the sum of squared differences
        #' between the unrounded and the rounded CN, as we expect the CN values
        #' to be as close to integers as possible. If this is not the case,
        #' the penalty for this average CN from the grid will be higher.
        #' Next, we multiply times the length of the segment
        chrseg<-brkpts[brkpt,1]
        beginseg<-brkpts[brkpt,2]
        endseg<-brkpts[brkpt,3]
        beginseg<-which(seglogR$chr==chrseg & seglogR$start==beginseg)
        endseg<-which(seglogR$chr==chrseg & seglogR$end==endseg)
        length<-(endseg-beginseg+1)
        sum<-sum+(sum((segCN[beginseg:endseg]-CN[beginseg:endseg])^2)*length)
    }
    #' Finally, we divide the total penalty value by the length of the entire genome
    pen<-sum/length(CN)
    #' and store this penalty for this specific sample and average CN
    pens<-c(pens,pen)
}


scores$score<-pens
colnames(scores)<-c("x","y")
scores$y <- smooth.spline(scores$y ~ scores$x)$y
a<-modes(scores)

# final <- scores$ploidygrid[which(scores$score==min(scores$score))]
final <- a[which(a$y==min(a$y)),1]
#' we take the global minimum to store in selected,
selected[1,] <- final
a$check <- ifelse(a$x==final, 4,3)
a$col <- ifelse(a$x==final, "tomato","steelblue")




#' we make a CN plot using the global minimum for transformation
CN<-pmax((2^logR[,4])*final,0)
segCN<-pmax(round((2^seglogR[,4])*final),0)
df<-data.frame(chr=seglogR$chr, start=seglogR$start,end=seglogR$end,CN=CN,segCN=segCN)
if(species=="Human"){
df$chr <- factor(df$chr, levels=c(seq(1,22),"X","Y"))
} else if(species=="Drosophila"){
df$chr<-factor(df$chr, levels=c("2L","2R","3L","3R","4","X"))
} else if(species=="Mouse"){
df$chr <- factor(df$chr, levels=c(seq(1,19),"X","Y"))
}


df1 <- df[df$chr %in% c("13","22"),]

q <- ggplot(df1, aes(x=start, y=CN)) + geom_point(size=1)+ylim(0,CNylim)
q <- q + geom_segment(data=df, aes(x = start, xend = end, y = segCN, yend= segCN),colour="red",size=2,alpha=0.5)
q <- q + theme_bw() + ylab("Copy Number") + xlab("Position") + ggtitle(paste0("ploidy = ", final , " sample = ",name))
q <- q + facet_wrap(~ chr,scales="free_x",ncol=2)

pdf("fig2_C.pdf",width=18,height=12)
plot(q)
dev.off()


# bin plot - pseudobulks (D)
name <- "T28_C_Pico"

seglogR <- seglogR[,c("chr","start","end",name)]
logR <- logR[,c("chr","start","end",name)]


pls <- seq(1.2,6,0.01)

#' set up two matrices, 'selected' to store the optimal result for each sample
#' (i.e. the global minimum of the grid)
selected <- matrix(ncol = 1, nrow=1)


scores <- data.frame(ploidygrid=pls, score=NA)

#' find the breakpoints and thus the segments for each sample

brkpts<-findbrkpts(seglogR[,c(1:4)])[[1]]
pens<-c()
for(pl in pls){
    #' for each possible average CN in our 1.2-6 grid,
    #' we transform the segmented logR to CN values, once rounded
    #' and once unrounded.
    CN<-pmax((2^seglogR[,4])*pl,0)
    segCN<-pmax(round((2^seglogR[,4])*pl),0)
    sum<-0
    for(brkpt in 1:nrow(brkpts)){
        #' Then, for each segment, we calculate the sum of squared differences
        #' between the unrounded and the rounded CN, as we expect the CN values
        #' to be as close to integers as possible. If this is not the case,
        #' the penalty for this average CN from the grid will be higher.
        #' Next, we multiply times the length of the segment
        chrseg<-brkpts[brkpt,1]
        beginseg<-brkpts[brkpt,2]
        endseg<-brkpts[brkpt,3]
        beginseg<-which(seglogR$chr==chrseg & seglogR$start==beginseg)
        endseg<-which(seglogR$chr==chrseg & seglogR$end==endseg)
        length<-(endseg-beginseg+1)
        sum<-sum+(sum((segCN[beginseg:endseg]-CN[beginseg:endseg])^2)*length)
    }
    #' Finally, we divide the total penalty value by the length of the entire genome
    pen<-sum/length(CN)
    #' and store this penalty for this specific sample and average CN
    pens<-c(pens,pen)
}


scores$score<-pens
colnames(scores)<-c("x","y")
scores$y <- smooth.spline(scores$y ~ scores$x)$y
a<-modes(scores)

# final <- scores$ploidygrid[which(scores$score==min(scores$score))]
final <- a[which(a$y==min(a$y)),1]
#' we take the global minimum to store in selected,
selected[1,] <- final
a$check <- ifelse(a$x==final, 4,3)
a$col <- ifelse(a$x==final, "tomato","steelblue")




#' we make a CN plot using the global minimum for transformation
CN<-pmax((2^logR[,4])*final,0)
segCN<-pmax(round((2^seglogR[,4])*final),0)
df<-data.frame(chr=seglogR$chr, start=seglogR$start,end=seglogR$end,CN=CN,segCN=segCN)
if(species=="Human"){
df$chr <- factor(df$chr, levels=c(seq(1,22),"X","Y"))
} else if(species=="Drosophila"){
df$chr<-factor(df$chr, levels=c("2L","2R","3L","3R","4","X"))
} else if(species=="Mouse"){
df$chr <- factor(df$chr, levels=c(seq(1,19),"X","Y"))
}


df1 <- df[df$chr %in% c("13","22"),]

q <- ggplot(df1, aes(x=start, y=CN)) + geom_point(size=1)+ylim(0,CNylim)
q <- q + geom_segment(data=df, aes(x = start, xend = end, y = segCN, yend= segCN),colour="red",size=2,alpha=0.5)
q <- q + theme_bw() + ylab("Copy Number") + xlab("Position") + ggtitle(paste0("ploidy = ", final , " sample = ",name))
q <- q + facet_wrap(~ chr,scales="free_x",ncol=2)

pdf("fig2_D.pdf",width=18,height=12)
plot(q)
dev.off()


# focal amplification (E)

selected_path <- "cnv_results10kb/Matrices/selected.txt"
segCN_path <- "cnv_results10kb/Matrices/CNseg.SPCF.35_10kb.txt"
locationSampleDescr <- "samplesheet.csv"
species <- "Human"

FAmergebins <- 10

#' read in the CN data and set the names of the chr accordingly
segCN <- read.csv(segCN_path, header=T, row.names=1, check.names = T)
segCN$chr<-factor(segCN$chr,levels = c(as.character(1:22),"X",'Y'))

selected <- read.csv(selected_path, header = T, row.names = 1)
ploidy <- round(median(selected$X))


sampleinfo <- read.csv(locationSampleDescr, header=T,sep="\t")
sampleinfo <- sampleinfo[sampleinfo$method =="Gtag",]
rownames(sampleinfo) <- sampleinfo$X
sampleinfo <- sampleinfo[rownames(sampleinfo) %in% colnames(segCN),]


sampleinfo <- sampleinfo %>% drop_na(subclones)

sampleinfo$cluster <- sampleinfo$subclones
Subclone <- c("#01a087","#7eeedc","#91d1c2","#3c5488","#7fa2be","#8391b4","#7f6147","#db9467","#f49b7f")
Empty <- c("#01a087","#91d1c2","#3c5488","#8391b4","#7d6148","#f49b7f","#A65628","#F5DEB3","#F5DEB3")

colscheme <- data.frame(Subclone,Empty,Empty)

interest1 <- colnames(colscheme)[1]
interest2 <- colnames(colscheme)[2]
interest3 <- colnames(colscheme)[3]
collist1 <- colscheme[,1]
collist2 <- colscheme[,2]
collist3 <- colscheme[,3]

sampleinfo$Subclone <- paste0(sampleinfo$subclones, " ", sampleinfo$timepoint)
sampleinfo$Subclone <- as.factor(sampleinfo$Subclone)

selected <- selected[rownames(selected) %in% rownames(sampleinfo),]
segCN <- segCN[,colnames(segCN) %in% c("chr","start","end",rownames(sampleinfo))]

#' @author Michiel
#' @description calculates a G-score for each bin. A G-score takes the value of 
#' 1 or 0 (depending on when it reaches the threshold) times the CN of that bin.
#' Thresholds are set as the median for each chromosome arm. 
#' @param dataframe the dataframe containing the CN profiles
#' @return the Gscores per bin
calculateGscore <- function(dataframe){
  GscoreMatrix <- data.frame(dataframe[,c(1,2,3)])
  for(i in 4:ncol(dataframe)){
    sample<-dataframe[,c(1,2,3,i)]
    scores<-c()
    ploidy <- selected[colnames(dataframe)[i],"X"]
    for(j in unique(sample$chr)){
      samplechr<-sample[sample$chr==j,]
      for(row in 2:nrow(samplechr)){
        if(samplechr[row,2]-samplechr[row-1,3]!=1)
          centromere <- row
      }
      if(exists("centromere")){
        pbins <- 1:(centromere-1)
        qbins <- centromere:nrow(samplechr)
      } else {
        pbins <- 1:nrow(samplechr)
        qbins <- NA
      }
      pmedian<-median(samplechr[pbins,4])
      if(is.na(pmedian)){pmedian<-0}
      qmedian<-median(samplechr[qbins,4])
      if(is.na(qmedian)){qmedian<-0}
      samplethreshold <- max(pmedian, qmedian)+ploidy
      identities <- samplechr[,4]>samplethreshold
      scores<-c(scores,samplechr[,4]*identities)
    }
    GscoreMatrix[,i]<-scores
  }
  colnames(GscoreMatrix)<-colnames(dataframe)
  GscoreMatrix$Gscores<-0
  GscoreMatrix$Gscores<-rowSums(GscoreMatrix[,-c(1,2,3)])
  return(GscoreMatrix)
}


GscoreMatrix<-calculateGscore(segCN)
write.csv(GscoreMatrix,"GscoreMatrix.csv")

#' set up a matrix to calculate the significance of the G-scores
pvaluematrix <- data.frame(GscoreMatrix[,c(1,2,3)])
#' for each sample, calculate the relative p-values and store them
for(i in 4:(ncol(GscoreMatrix)-1)){
  finaltest <- data.frame(GscoreMatrix[,i])
  finaltest$pvalue <- 0
  for(j in 1:nrow(finaltest)){
    finaltest$pvalue[j] <- sum(finaltest$GscoreMatrix...i.>=finaltest$GscoreMatrix...i.[j])/nrow(finaltest)
  }
  pvaluematrix[,i]<-finaltest$pvalue
}


pvaluematrix <- pvaluematrix[order(pvaluematrix$chr, as.numeric(as.character(pvaluematrix$start))),]
write.csv(pvaluematrix,"pvaluematrix.csv")


#' Use fishers method to combine the p-values per bin over the samples
fisherpvalues <- data.frame(combinedpvalue=NA)
for(i in 1:nrow(pvaluematrix)){
  fisherpvalues[i,] <-pchisq(-2 * sum(log(as.numeric(pvaluematrix[i,-c(1,2,3)]))),df=2*length(pvaluematrix[i,-c(1,2,3)]),lower=FALSE)
}

write.csv(fisherpvalues,"fisherpvalues.csv")


rownames(fisherpvalues)<-1:nrow(fisherpvalues)
#' peaks (FA) are those bins with a fisher pval below 0.1
peaks <- data.frame(peakpositions=as.numeric(rownames(fisherpvalues)[fisherpvalues$combinedpvalue<0.1]))

#' Before we plot the FA heatmap, we first merge bins that are too close to 
#' eachother, as they are in fact one FA but with a neutral or lower CNV value.
if(dim(peaks)[1]!=0){
  peaks$diff<-0
  for(i in 2:nrow(peaks)){
    peaks$diff[i]<-peaks$peakpositions[i]-peaks$peakpositions[i-1]
  }
  tomerge <- which(peaks$diff %in% c(2:FAmergebins))
  peaks$diff <- NULL
  if(length(tomerge)!=0){
    for(i in 1:length(tomerge)){
      window <- data.frame(peakpositions=peaks$peakpositions[tomerge[i]-1]:peaks$peakpositions[tomerge[i]+1])
      peaks <- rbind(peaks, window)
    }
  }
}

peaks <- unique(peaks)
peaks <- peaks[order(peaks$peakpositions),]
peaks <- data.frame(peakpositions=peaks)


#' The code then continues to extract the determined peaks from the CN profiles,
#' and creates a heatmap of the areas.

  peaks$diff<-0
  for(i in 2:nrow(peaks)){
    peaks$diff[i]<-peaks$peakpositions[i]-peaks$peakpositions[i-1]
  }
  newsect <- c(which(peaks$diff != 1), nrow(peaks)+1)
  ranges<-list()
  sizes <- c()
  sizesString <- c()
  starts <- c()
  ends <- c()
  # truepeaks<-c()
  for(i in 1:(length(newsect)-1)){
    # truepeaks[i] <- floor(median(peaks[newsect[i]:(newsect[i+1]-1),1]))
    ranges[[i]]<-(peaks$peakpositions[newsect[i]]-10):(peaks$peakpositions[newsect[i+1]-1]+10)
    sizes[i] <- pvaluematrix[peaks$peakpositions[newsect[i+1]-1],3]-pvaluematrix[peaks$peakpositions[newsect[i]],2]
    sizesString[i] <- paste0(toString(RoundTo(sizes[i], multiple = 1000)/1000), "kb")
    starts[i] <- segCN[peaks$peakpositions[newsect[i]],2]
    ends[i] <- pvaluematrix[peaks$peakpositions[newsect[i+1]-1],3]
  }



ranges <- ranges[1:9]
ranges[[10]] <- c(7655:7685)
ranges[[11]] <- c(9362:9440)


sizes <- c(1052488,365411,597107,710962,772204,1203865,583332,157316,415539,116075,677697)
sizesString <- c('1052kb\n13q13.3A','365kb\n13q13.3B','597kb\n13q14.2-3','711kb\n13q14.3','772kb\n13q21.2','1204kb\n13q21.33','583kb\n13q22.3','157kb\n13q31.1','416kb\n13q31.3','116kb\n13q32.3','678kb\n22q11.21')
starts <- c(35866758,38623156,50777395,53415279,61400237,68908352,77466797,81477365,91959299,100071299,20804701)
ends <- c(36919246,38988567,51374502,54126241,62172441,70112217,78050129,81634681,92374838,100187374,21482398)



df <- data.frame(starts,ends)
df$midden <- df$starts + (df$ends - df$starts)/2


  filter<-unlist(ranges)
  CNheatmapmatrix<-segCN[filter,]



  sections <- c()
  for(i in 1:length(ranges)){
    sections <- c(sections, rep(toString(i), length(ranges[[i]])))
  }
  split <- data.frame(sections=sections, chr=CNheatmapmatrix$chr)
  #dir.create(paste0(ExportPath,"/Focals"))
  for(i in 1:(ncol(CNheatmapmatrix)-3)){
    x<-data.frame(chr=paste0("chr",CNheatmapmatrix$chr),
                  start=CNheatmapmatrix$start,
                  end=CNheatmapmatrix$end,
                  name=colnames(CNheatmapmatrix)[i+3],
                  copy_number = CNheatmapmatrix[,i+3])
    #write.table(x,paste0(ExportPath,"/Focals/",colnames(CNheatmapmatrix)[i+3],".focals.csv"), 
    #            sep = '\t', col.names = F, row.names=F, quote = F)
  }
  CNheatmapmatrix <- CNheatmapmatrix[,-c(1:3)]
  CNheatmapmatrix <- CNheatmapmatrix[,colnames(CNheatmapmatrix)[order(colnames(CNheatmapmatrix))]]
  



  filter<-unlist(ranges)
  CNheatmapmatrix<-segCN[filter,]
  sections <- c()
  for(i in 1:length(ranges)){
    sections <- c(sections, rep(toString(i), length(ranges[[i]])))
  }
  split <- data.frame(sections=sections, chr=CNheatmapmatrix$chr)
  for(i in 1:(ncol(CNheatmapmatrix)-3)){
    x<-data.frame(chr=paste0("chr",CNheatmapmatrix$chr),
                  start=CNheatmapmatrix$start,
                  end=CNheatmapmatrix$end,
                  name=colnames(CNheatmapmatrix)[i+3],
                  copy_number = CNheatmapmatrix[,i+3])

  }
  CNheatmapmatrix <- CNheatmapmatrix[,-c(1:3)]
  CNheatmapmatrix <- CNheatmapmatrix[,colnames(CNheatmapmatrix)[order(colnames(CNheatmapmatrix))]]


  pft <- c()
  for(i in unique(split$sections)){
    pft <- c(pft,as.character(tail(split[split$sections==i,2],1)))
  }
  lol <- as.data.frame(unique(pft))
  lol$count <- 0
  for(i in 1:nrow(lol)){
    lol$count[i]<-length(which(pft==lol[i,1]))
  }
  if(nrow(lol)<=2){
    xx <- brewer.pal(n=3,"Set2")
    xx <- xx[1:(nrow(lol))]
    lol$col <- xx
  } else {
    lol$col <- brewer.pal(n=nrow(lol),"Set2")
  }

 unique(split$sections)
 pft
 lol$col <- c("#57b4e9","#e79f00")
 lol
  
  blabla <- c()
  for(i in 1:nrow(lol)){
    blabla <- c(blabla, rep(lol$col[i],lol$count[i]))
  }
  names <- c()
  for(i in 1:nrow(lol)){
    names <- c(names, as.character(rep(lol$`unique(pft)`[i],lol$count[i])))
  }
  names <- paste("chr",names)
  column_ha <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = blabla), labels  = names, 
                                                  labels_gp = gpar(col = "black", fontsize = 10)), show_legend = T, show_annotation_name = T)
  
  col_fun2 = colorRamp2(c(min(selected[,1]),ploidy,max(selected[,1])), c("#0571b0","grey96", "tomato"))
  col_fun2(seq(min(selected[,1]), max(selected[,1])))
  
  row_ha_right=HeatmapAnnotation(bar=anno_text(selected[,1], location = 0.5, just = "center",
                                               gp = gpar(fill = col_fun2(selected[,1]), col = "black", border = "black"),
                                               width = max_text_width(selected[,1])*1.2), 
                                 which="row", show_annotation_name = F, show_legend = T,
                                 annotation_legend_param = list(title="\nSelected"),
                                 annotation_name_rot = 90)
  
  if (interest1 == "Empty") {
    row_ha=NULL
  } else if (interest2 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfo[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfo[interest1][,1])[order(unique(sampleinfo[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    row_ha = rowAnnotation(int1 = sampleinfo[interest1][, 1],show_legend = T,show_annotation_name = F,
                           col = uniqint1, annotation_legend_param = list(int1 = list(title = interest1)),
                           annotation_name_rot = 90)
  } else if (interest3 == "Empty") {
    uniqint1 <- collist1[1:length(unique(sampleinfo[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfo[interest1][,1])[order(unique(sampleinfo[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfo[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfo[interest2][,1])[order(unique(sampleinfo[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    row_ha = rowAnnotation(int1 = sampleinfo[interest1][, 1],int2 = sampleinfo[interest2][, 1],
                           show_legend = T,show_annotation_name = F,col = c(uniqint1, uniqint2), 
                           annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2))),
                           annotation_name_rot = 90)
  } else {
    uniqint1 <- collist1[1:length(unique(sampleinfo[interest1][,1]))]
    names(uniqint1) <- as.character(unique(sampleinfo[interest1][,1])[order(unique(sampleinfo[interest1][,1]))])
    uniqint1 <- list(int1=uniqint1)
    
    uniqint2 <- collist2[1:length(unique(sampleinfo[interest2][,1]))]
    names(uniqint2) <- as.character(unique(sampleinfo[interest2][,1])[order(unique(sampleinfo[interest2][,1]))])
    uniqint2 <- list(int2=uniqint2)
    
    uniqint3 <- collist3[1:length(unique(sampleinfo[interest3][,1]))]
    names(uniqint3) <- as.character(unique(sampleinfo[interest3][,1])[order(unique(sampleinfo[interest3][,1]))])
    uniqint3 <- list(int3=uniqint3)
    
    row_ha = rowAnnotation(
      int1 = sampleinfo[interest1][, 1],int2 = sampleinfo[interest2][, 1],int3 = sampleinfo[interest3][, 1],
      show_legend = T,show_annotation_name = F, col = c(uniqint1, uniqint2, uniqint3), 
      annotation_legend_param = list(int1 = list(title = interest1), int2 = list(title = paste0("\n", interest2)),int3 = list(title = paste0("\n", interest3))),
      annotation_name_rot = 90)
  }
  col_fun = colorRamp2(c(0, 5, 5.1,10,20,max(CNheatmapmatrix)), c("grey98","grey98", "grey48","grey48","black","darkred"))
  col_fun(seq(min(CNheatmapmatrix), max(CNheatmapmatrix)))
  lgd <- list(col_fun = col_fun, title = "\nCopy number", at = c(0, 5, 5.1,10,20,max(CNheatmapmatrix)),
              labels = c(0, 5, 6,10,20,max(CNheatmapmatrix)))
  
  print(Heatmap(t(as.matrix(CNheatmapmatrix)), name = "Copy number",heatmap_legend_param = lgd,
                column_split = factor(split$sections, levels=1:length(ranges)), cluster_columns = F, show_column_names = F, 
                show_row_names = F, show_heatmap_legend = T, cluster_column_slices = FALSE,
                column_title = sizesString, column_title_side = "bottom", column_title_gp = gpar(col = "black", fontsize = 10),
                top_annotation = column_ha, left_annotation = row_ha, 
                col=col_fun, column_title_rot = 0,right_annotation = row_ha_right,
                clustering_method_rows = "ward.D2"))

  

no_categories <- Heatmap(t(as.matrix(CNheatmapmatrix)), name = "Copy number",heatmap_legend_param = lgd,
                column_split = factor(split$sections, levels=1:length(ranges)), cluster_columns = F, show_column_names = F, 
                show_row_names = F, show_heatmap_legend = T, cluster_column_slices = FALSE,
                column_title = sizesString, column_title_side = "bottom", column_title_gp = gpar(col = "black", fontsize = 10),
                top_annotation = column_ha, left_annotation = row_ha, 
                col=col_fun, column_title_rot = 0,right_annotation = row_ha_right,
                clustering_method_rows = "ward.D2")



orderNames <- rownames(t(as.matrix(CNheatmapmatrix)))[row_order(no_categories)]

CNheatmapmatrix_2 <- ifelse(CNheatmapmatrix <= 5, "1-5", 
                     ifelse(CNheatmapmatrix > 5 & CNheatmapmatrix <= 15, "6-15",
                            ifelse(CNheatmapmatrix > 15 & CNheatmapmatrix <= 25, "16-25",
                                   ifelse(CNheatmapmatrix > 25 & CNheatmapmatrix <= 35, "26-35", ">35"))))

colorshtmp = structure(c("#fafafa", "#d9d9d9", "#8c8c8c", "#404040", "#000000"), names = c("1-5","6-15","16-25","26-35",">35"))

CNheatmapmatrix_2 <- CNheatmapmatrix_2[,orderNames]

    row_ha = rowAnnotation(int1 = sampleinfo[orderNames,interest1],show_legend = T,show_annotation_name = F,
                           col = uniqint1, annotation_legend_param = list(int1 = list(title = interest1)),
                           annotation_name_rot = 90)



Heatmap(t(as.matrix(CNheatmapmatrix_2)), 
        name = "Copy number",
        #heatmap_legend_param = lgd,
        column_split = factor(split$sections, levels=1:length(ranges)), 
        cluster_columns = F, 
        cluster_rows=T,
        show_column_names = F, 
        show_row_names = F, 
        show_heatmap_legend = T, 
        cluster_column_slices = FALSE,
        column_title = sizesString, 
        column_title_side = "bottom", 
        column_title_gp = gpar(col = "black", fontsize = 10),
        top_annotation = column_ha, 
        left_annotation = row_ha, 
        col=colorshtmp, 
        column_title_rot = 0,
        #split=cluster$cluster
        #right_annotation = row_ha_right,
        clustering_method_rows = "ward.D2"
       )


cat <- Heatmap(t(as.matrix(CNheatmapmatrix_2)), 
        name = "Copy number",
        #heatmap_legend_param = lgd,
        column_split = factor(split$sections, levels=1:length(ranges)), 
        cluster_columns = F, 
        cluster_rows=T,
        show_column_names = F, 
        show_row_names = F, 
        show_heatmap_legend = T, 
        cluster_column_slices = FALSE,
        column_title = sizesString, 
        column_title_side = "bottom", 
        column_title_gp = gpar(col = "black", fontsize = 10),
        top_annotation = column_ha, 
        left_annotation = row_ha, 
        col=colorshtmp, 
        column_title_rot = 0,
        #split=cluster$cluster
        #right_annotation = row_ha_right,
        clustering_method_rows = "ward.D2"
       )



  pdf('focals_Categorical.pdf', height = 6, width = 15)
  cat
  dev.off()

  pdf("focals_NONCategorical.pdf', height = 6, width = 15)
  no_categories
  dev.off()