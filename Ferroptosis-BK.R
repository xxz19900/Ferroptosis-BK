setwd("/Users/xizhanxu/Learning/Tongren/Projects/bkfk/")

############################################################################===
#####  01. PCA analysis 
############################################################################===

## Load libraries
rm(list = ls())  
options(stringsAsFactors = F)
library(ggplot2)
library(ggstatsplot)
library(patchwork)
library(ggplotify)
library(cowplot)
library("FactoMineR")
library("factoextra")
library(scales)
library(ggforce)


## Data input
load(file = 'data/step1-GSE58291_output.Rdata')
table(pd$title)
pd <- pd[pd$title %in% c("Normal", "BK"), ]
pd$title <- factor(pd$title, levels = c("Normal", "BK"))
GSE58291_anno <- GSE58291_anno[,pd$geo_accession]

## PCA plot
group_list <- factor(pd$title, levels = c("Normal", "BK"))
GSE58291_anno[1:4,1:4]

exp <- GSE58291_anno

exp=t(exp) # row: samples, column: genes
exp=as.data.frame(exp)
dat.pca <- PCA(exp , graph = FALSE)
this_title <- paste0('PCA')

p1 <- fviz_pca_ind(dat.pca,
                   geom="point",
                   geom.ind = "point", # show points only (nbut not "text")
                   col.ind = group_list, # color by groups
                   palette = c(hue_pal()(3)[2], hue_pal()(3)[1]),
                   addEllipses = TRUE, # Concentration ellipses
                   legend.title = "Group",mean.point = FALSE)+
  ggtitle(this_title)+
  theme_cowplot()+
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  coord_equal()

p1
ggsave("figure/BK_vs_Normal_PCA.pdf", p1, width = 5.63, height = 5.69, units = "in")

############################################################################===
#####  02. DEG analysis
#####   Volcano plot (BK vs. Normal)
############################################################################===

## Load library
library(limma)

# Data
group <- factor(pd$title, levels = c("Normal", "BK"))
expr <- GSE58291_anno
contrast <- paste(c("BK","Normal"),collapse = "-")

# limma for DEG analysis
design <- model.matrix( ~ 0 + group) #design
colnames(design) <- levels(group)
contrast.matrix <- makeContrasts(contrast, levels = design)
fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

# results
options(digits = 4)
GSE58291_DEG_BK_Normal <- topTable(fit, coef = 1, n = Inf) %>% 
  rownames_to_column(var = "genesymbol") 

# Selection of DEGs
GSE58291_DEG_p_BK_Normal <- GSE58291_DEG_BK_Normal %>% 
  dplyr::filter(abs(logFC) >= 1, adj.P.Val <= 0.05)

nrDEG=GSE58291_DEG_BK_Normal
head(nrDEG)
attach(nrDEG)
df=nrDEG
df$v= -log10(adj.P.Val) 
# cut_logFC <- with(deg,mean(abs(deg$logFC)) + 2*sd(abs(deg$logFC)) )
# cut_logFC
logFC_t = 1

# Gene state
df$State=ifelse(df$adj.P.Val>0.05,'Stable', 
                ifelse( df$logFC >logFC_t,'Up', 
                        ifelse( df$logFC < -logFC_t,'Down','Stable') ))

# number of DEGs
table(df$State)

# plot title
this_tile <- paste0("BK vs. Normal",
                    '\nThe number of up gene is ',nrow(df[df$State == 'Up',]) ,
                    '\nThe number of down gene is ',nrow(df[df$State == 'Down',])
)

#Volcano plot
library(cowplot)

p2 <- ggplot(data = df, aes(x = logFC, y = v)) +
  geom_point(alpha=0.6, size=1.5, 
             aes(color=State)) +
  ylab("-log10(adjusted P value)")+
  xlab("log2(Fold change)")+
  scale_color_manual(values=c("#34bfb5", "#828586","#ff6633"))+
  geom_vline(xintercept= c(-1,1),lty=4,col="grey",lwd=0.8) +
  geom_hline(yintercept= -log10(0.05),lty=4,col="grey",lwd=0.8) +
  theme_cowplot()+
  ggtitle(this_tile )+
  theme(plot.title = element_text(size=12,hjust = 0.5),
        legend.title = element_blank())+
  coord_equal()

p2
ggsave("figure/GSE58291_BK_normal_volcano.pdf", p2, width = 5.63, height = 5.69, units = "in")

### DEGs
GSE58291_DEG_BK_Normal$State=ifelse(GSE58291_DEG_BK_Normal$adj.P.Val>0.05,'Stable', 
                                    ifelse( GSE58291_DEG_BK_Normal$logFC >logFC_t,'Up', 
                                            ifelse( GSE58291_DEG_BK_Normal$logFC < -logFC_t,'Down','Stable') ))
GSE58291_DEG_p_BK_Normal$State=ifelse(GSE58291_DEG_p_BK_Normal$adj.P.Val>0.05,'Stable', 
                                      ifelse( GSE58291_DEG_p_BK_Normal$logFC >logFC_t,'Up', 
                                              ifelse( GSE58291_DEG_p_BK_Normal$logFC < -logFC_t,'Down','Stable') ))

save(GSE58291_DEG_BK_Normal,GSE58291_DEG_p_BK_Normal, file = 'data/step2-GSE58291_deg_BK_vs_Normal.Rdata')
write.csv(GSE58291_DEG_p_BK_Normal, "data/GSE58291_DEG_p_BK_Normal.csv")

############################################################################===
#####  03. Statistical analysis
############################################################################===
rm(list = ls())
load(file = 'data/step2-GSE58291_deg_BK_vs_Normal.Rdata')
load(file = 'data/step1-GSE58291_output.Rdata')
table(pd$title)
pd <- pd[pd$title %in% c("Normal", "BK"), ]
pd$title <- factor(pd$title, levels = c("Normal", "BK"))
GSE58291_anno <- GSE58291_anno[,pd$geo_accession]

### Ferroptosis gene set
library(openxlsx)
library(tableone)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggsci)
library(scales)

### DEG & Ferroptosis

## Statistical analysis
## Ferroptosis
ferr <- read.xlsx("data/Ferroptosis-gene-list.xlsx")
ferr <- ferr[,c(1,3)]
colnames(ferr) <- c("from", "to")

ferrgene<- ferr$to[1:114]
ferrdb <- GSE58291_anno[ferrgene,]
ferrdb <- t(ferrdb) %>% as.data.frame()
ferrdb$group <- pd$title
ferrdb$ID <- rownames(ferrdb)
ferrdb <- ferrdb[, !(colnames(ferrdb) %in% c("NA", "NA.1","NA.2","NA.3","NA.4","NA.5"))]
ferrdb <- ferrdb[, c(110,109, 1:108)]
colnames(ferrdb)[13] <- "ALOX12"

immunecell_var <- colnames(ferrdb)[2:110]
cat_immunecell_var <- c("group")

## Statistical analysis
library("EasyStat")
library("tidyverse")

norCv = MuiNorCV(data = ferrdb,num = c(3:110), method_cv = "leveneTest")
norCv <- as.data.frame(norCv)
norCv$cor <- as.logical(norCv$cor)
norCv$CV <- as.logical(norCv$CV)

norCv[norCv$cor & norCv$CV,]$DI
cytokine_var_nonmormal <- norCv[!(norCv$cor & norCv$CV),]$DI

tab_immunecell <- CreateTableOne(vars = immunecell_var, strata = "group" , data = ferrdb, factorVars = cat_immunecell_var, addOverall = TRUE)
tab2 <- print(tab_immunecell,  formatOptions = list(big.mark = ","), nonnormal = cytokine_var_nonmormal, exact = "group", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)
summary(tab_immunecell)

write.csv(tab2, file = "data/Ferrdb_stat.csv")

## Heatmap of Ferroptosis genes
library(pheatmap)
ferr <- ferr[1:114,]

immunecell <- ferrdb[c(3,7,10,12,14,16,19,1:2,4:6,8:9,11,13,15,17:18),]
immunecell <- immunecell[,-c(1:2)]
immunecell <- t(immunecell)

### 114 genes
rownames(ferr) <- ferr$to
ferr <- ferr[rownames(immunecell), ]
ferr <- ferr[order(ferr$from),]

ferr <- ferr[c(5:108,1:4),]
immunecell <- immunecell[ferr$to,]
immunecell <- immunecell[, c(8:19,1:7)]

# annotation
annotation_col <- data.frame(Group = c(rep("Control",12), rep("BK",7)), row.names = colnames(immunecell))
annotation_col$Group <- factor(annotation_col$Group, levels = c("Control", "BK"))
annotation_row <- data.frame(Category =  ferr$from, row.names = ferr$to)
annotation_row$Category <- factor(annotation_row$Category, levels = c("Iron metabolism","Lipid metabolism","Oxidant-reductant","ESCRT-III"))

# colors
library(ggsci)
mypal = pal_lancet("lanonc", alpha = 0.6)(9)
show_col(pal_nejm(alpha = 0.6)(8))
library("scales")
show_col(mypal)
mypal

anno_colors <- list(Group = c(Control = "#00BA38", BK = "#F8766D"),
                    Category = c(`Iron metabolism` = "#20854E99", 
                                 `Lipid metabolism` = "#E1872799", `Oxidant-reductant` = "#BC3C2999", 
                                 `ESCRT-III` = "#0072B599"))

##### Circular heatmap

mat1 <- t(scale(t(immunecell)))
split <- annotation_row$Category

library(circlize) # >= 0.4.10
range(mat1)
col_fun1 = colorRamp2(c(-3.7, -0.05, 3.6), c("navy", "white", "firebrick3"))

library(dendextend)
dend_col = structure(1:4, names = levels(split))

circos.par(gap.after = c(2, 2, 2,20))
circos.heatmap(mat1, split = split, col = col_fun1, dend.side = "inside",
               rownames.side = "outside",
               rownames.cex=0.9,
               rownames.font=0.8,
               rownames.col="black",
               dend.track.height = 0.1,
               dend.callback = function(dend, m, si) {
                 # when k = 1, it renders one same color for the whole dendrogram
                 color_branches(dend, k = 1, col = dend_col[si])
               })

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 4) { # the last sector
    circos.rect(CELL_META$cell.xlim[2]+0.1 + convert_x(1, "mm"), 0,
                CELL_META$cell.xlim[2]+0.1 + convert_x(5, "mm"), 7,
                col = "orange", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 3.5,
                "group 1", cex = 0.5, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 7,
                CELL_META$cell.xlim[2] + convert_x(5, "mm"), 19,
                col = "pink", border = NA)
    circos.text(CELL_META$cell.xlim[2] + convert_x(3, "mm"), 13,
                "group 2", cex = 0.5, facing = "clockwise")
  }
}, bg.border = NA)
circos.clear()

library(ComplexHeatmap)
lgd = Legend(title = "mat1", col_fun = col_fun1)
grid.draw(lgd)

### Selected Ferroptosis gene
sigferrgene <- c("ALOX5", "ACSL4", "GPX2", "GPX4", "NCOA4", "SAT1", "TF", "TFRC", "FTL", "IREB2", "SOCS1", "TP53", "SLC7A11", "SLC1A5", "SLC38A1")
ferrcat <- ferr[ferr$to %in% sigferrgene, ]
ferrcat <- ferrcat[order(ferrcat$from),]

sigferr <- GSE58291_anno[ferrcat$to, c(3,7,10,12,14,16,19,1:2,4:6,8:9,11,13,15,17:18)]
sigferr <- t(sigferr)

major_ferr <- sigferr %>% as.data.frame()
major_ferr$ID <- rownames(major_ferr)
major_ferr$Group <- c(rep("BK",7),rep("Control",12))

library(reshape2)
data <- melt(major_ferr,id.vars = c("ID","Group"))
colnames(data)=c("ID","Group","Gene","Expression")
data$Category <- rep(ferrcat$from, each=19)
data$Gene <- factor(data$Gene, levels = ferrcat$to)
data$Category <- factor(data$Category, levels = c("Iron metabolism","Lipid metabolism","Oxidant-reductant"))
data$Group <- factor(data$Group, levels = c("Control", "BK"))
##
library(RColorBrewer)
library(cowplot)

data$Expression <- log2(data$Expression+1)

p3 <- ggplot(data = data, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.6, lwd=0.5) +
  theme_cowplot()+
  scale_fill_manual("Group", values = c("#00BA38","#F8766D"))+
  ylab("Gene expression")+xlab(NULL)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+
  facet_grid(~ Category, space = "free", scales = "free")

ggsave("figure/Ferroptosis-gene-boxplot.pdf", p3, width = 8.7, height = 4.19, units = "in")
