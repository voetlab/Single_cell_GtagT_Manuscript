# Figure 1: 

library(ggplot2)
library(ggpubr)
library(magrittr)
library(plyr)
library(grid)
library(cowplot)


## Costs (B)
dfc <- data.frame(cost=rep(c("a_seq", "b_DNA", "c_RNA","d_iso"), each=2),
                  method=rep(c("Gtag&T", "G&T"),4),
                  'cost_EUR'=c(3240.00,3240.00,
                               3965.50, 9618.19,
                               5567.48, 5567.48, 
                               1114.84, 1114.84))
dfc


p <- ggbarplot(dfc, x = "method", y = "cost_EUR",
               fill = "cost",
               # color = "cost", 
               palette = c("#CCCCCC","#7A7A7A","white", "black"),
               orientation = "horiz",
               ylab = "Cost for 384 samples (euro)",
               xlab = "\n"
)

q <- p + scale_fill_manual(labels = c("Sequencing",
                                      "DNA",
                                      "RNA",
                                      "Isolation & separation")
                           , values = c("black", "#7A7A7A","#CCCCCC", "white")) 


r <- q + theme(legend.key.size = unit(0.6,"line"))
t <- ggpar(r, font.legend = c(size = 9),
           legend = c(x=0.85, y =0.3)) + rremove("legend.title") +
  scale_y_continuous(breaks = seq(0, 20000, by = 4000)) 
costplot <- t + theme(axis.text=element_text(size=12), axis.title=element_text(size=10), strip.text = element_text( size=10), strip.background = element_rect(size=2))
costplot

dfc <- data.frame(cost=rep(c("a_DNA", "b_RNA", "c_iso"), each=2),
                  method=rep(c("Gtag&T", "G&T"),3),
                  'cost_time'=c(5, 11,
                                14.75, 14.75, 
                                0.5, 0.5 ))
dfc

p <- ggbarplot(dfc, x = "method", y = "cost_time",
               fill = "cost",
               palette = c("#CCCCCC","#7A7A7A", "white"),
               orientation = "horiz",
               ylab = "Time per experiment (hours)",
               xlab = "\n"
)

q <- p + scale_fill_manual(labels = c("Isolation & separation", "RNA","DNA"), values = c("black", "#7A7A7A","#CCCCCC")) 

q <- p + scale_fill_manual(labels = c(
  "DNA",
  "RNA",
  "Isolation & separation")
  , values = c("#7A7A7A","#CCCCCC", "white")) 

r <- q + theme(legend.key.size = unit(0.6,"line"))
t <- ggpar(r, font.legend = c(size = 9),
           legend = c(x=0.88, y =0.3)) + rremove("legend.title") + 
  scale_y_continuous(breaks = seq(0, 25, by = 5)) 
timeplot <- t + theme(axis.text=element_text(size=12), axis.title=element_text(size=10), strip.text = element_text( size=10), strip.background = element_rect(size=2))
timeplot

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(costplot)
grid.draw(mylegend)


combined <- plot_grid(costplot+ theme(legend.position="none"),mylegend,
                      timeplot + theme(legend.position="none"),ncol=2)



pdf(file = "combined.pdf", width = 12, height = 3)
combined
dev.off()

## Expression (C)
countMatrix <- as.data.frame(logcounts(sce))

names_Gtag #sample names
names_Picoplex #sample names

countGenes_G <- countGenes[names_Gtag]
mG <- as.data.frame(rowMeans(countGenes_G))

countGenes_P <- countGenes[names_Picoplex]
mP <- as.data.frame(rowMeans(countGenes_P))

colnames(mP) <- "picoPlex"
colnames(mG) <- "Gtag"

means <- merge(mG, mP, by = 'row.names', all = FALSE) 

mean <- ggplot(means,aes(x=Gtag,y=picoPlex))+ 
  geom_point(size=1)+ 
  geom_smooth(method = "lm", formula = y ~ x,colour="red") +
  stat_cor() +
  theme_bw() +
  ylab("G&T (mean log expression)") +
  xlab("Gtag&T (mean log expression)") +
  theme(axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


## MAPD (D)
mapd <- read.csv("MAPD_samples.txt", header=TRUE, row.names=1) #from CNV pipeline
samples <- read.csv("samples.csv")
mapd <- merge(mapd, samples, on="sample")

mapd %>% 
    dplyr::group_by(method) %>% 
    dplyr::summarise(median = median(mapd, na.rm = TRUE),mean=mean(mapd, na.rm = TRUE))


p <- ggviolin(mapd, x="method",y="mapd",
             fill="method", palette="npg",
             add="boxplot",add.params=list(fill="white"), ylab="\n MAPD \n",trim=T)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) +
    stat_compare_means(label.y=1,label.x=1.45,label="p.signif",size=6) + rremove("xlab")



my_comparisons <- list(c("Gtag", "picoPlex"),c("Gtag", "DNTR"),c("Gtag", "scONE"),c("Gtag","sciL3"))

q <- ggviolin(mapd, x="method",y="mapd",
             fill="method", palette="npg",
             add="boxplot",add.params=list(fill="white"), ylab="\n MAPD \n",trim=T)+
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.3)) +
    stat_compare_means(comparisons=my_comparisons, aes(label=..p.signif..)
                      ) + theme(legend.position = 'none')  + rremove("xlab")



## Lorenz Curve (E)

#   genomeCoverageBed -g Homo_sapiens.GRCh37.dna.primary_assembly.fa -bg -ibam file.name.bam | gzip > tmpLC/name.genomecov.txt.gz
#   zcat tmpLC/name.genomecov.txt.gz | sort -k4,4n -T tmpLC --compress gzip -S 4G | gzip > tmpLC/$name.genomecov.sorted.txt.gz
#   python src/python/compute_cumulative_cov.py tmpLC/name.genomecov.sorted.txt.gz resultsLC/name.cum_cov_stats.txt
#   head -n 1 resultsLC/name.cum_cov_stats.txt > resultsLC/name.cum_cov_sum_stats.txt
#   cat resultsLC/name.cum_cov_stats.txt |  awk '{{printf("%.2f\t %.2f\n ", $1, $2)}}' | sort | uniq | sed 's/ \+//g' | sed '/^$/d' >> resultsLC/name.cum_cov_sum_stats.txt


path_to_files = "lorenzCurve/resultsLC"
files <- list.files(path=path_to_files,pattern=".cum_cov_sum_stats.txt")


samples <- samples[c("sample","method")]

concat_BL = data.frame()
for (file in files){
  #print(file)
  d <- read.delim(paste(path_to_files,"/",file,sep=""))
  d$sample= strsplit(file,"\\.")[[1]][1]
  d$sample=strsplit(d$sample,"_dfp")[[1]][1]
  #d$type=strsplit(d$sample,"_")[[1]][1]
  concat_BL <- rbind(concat_BL,d)
}

concat_BL2 <- merge(concat_BL, samples, on="sample")


colrs <- c("Gtag"="#E64B35FF",
           "picoPlex"="#4DBBD5FF",
           "DNTR"="#00A087FF",
           "scONE"="#3C5488FF",
           "sciL3"="#F39B7FFF")


lcBL <- ggline(concat_BL2, x="cumfrac_bases_genome", y="cumfrac_bases_sequenced", 
               palette = colrs,              
               add = "mean_sd", 
               plot_type = "l",
               color = "method",   #
               error.plot = "linerange",
               add.params = list(size=0.5),
               #facet.by = "method",
               numeric.x.axis =TRUE, xlab = " Cumulative fraction of covered genome", ylab = "\n Cumulative fraction of mapped bases \n"
) + theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=12), 
          legend.key.size = unit(0.5,"line"), 
          strip.text = element_text( size=6, lineheight = 1), 
          strip.background = element_rect(size=1))