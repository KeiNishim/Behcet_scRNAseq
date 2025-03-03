library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(metap)
library(data.table)
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(icellnet)
library(gridExtra)
library(devtools)
library(circlize)
library(png)

img <- readPNG("Bechedt/Figure/Fig.s4b.png")
plot(1:2, type='n', xlab="", ylab="", xlim=c(0, 1), ylim=c(0, 1))
rasterImage(img, 0, 0, 1, 1)

#For plotting PBMC with ISG/monocyte signature###
all <- FindAllMarkers(subset(s.int,downsample=20000),only.pos = TRUE,
                          logfc.threshold = 0.25, min.pct=0.25, max.cells.per.ident=5000)

write.csv(all, "Genemarkers.csv")

top10 <- all %>% top_n(n = 100, wt = avg_log2FC)

top10 <- all %>% group_by(cluster) %>% top_n(n = -4, wt = p_val_adj)
top10 <- top10 %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)#4 for rhapsody, 3 for chromium
X <- unique(top10['gene'])
Y <- apply(X,2,rev)
pp <- DotPlot(s.int, features=c("TNF","IFNG"), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="clusters")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(1,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=0)+coord_flip()
pp

FeaturePlot(s.int, c("IFNG"),split.by="CRP")
#Read PBMC_ind as s.int
pp <- dittoHeatmap(subset(s.int,downsample=2000),
                   genes=rownames(mar2),
                   annot.by=c("disease","sample"),
                   order.by=c("disease","sample"),
                   scale.to.max=T,
                   slot="data",
                   border_color="black")
pp
ggsave('Fig.S10.png', pp, width=10, height=10)

#Enrichment score analysis
features <- list(c("CYBB","ITGAM","ITGB2","TLR2","TLR7","LILRB2"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="NETscore", slot="data")

features <- list(c("CYBB","IL1B","IL6","ITGAM","ITGB2","MMP9","PTAFR","TLR2","TLR7","TLR8","BST1","CD93","CSF3R","FCGR2B","FPR1","HPSE","LILRB2","PDE4B"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="NETosisscore", slot="data")

features <- list(c(
"B2M","CAMK2A","CAMK2B","CAMK2D","CAMK2G","CD44","CIITA","FCGR1A","FCGR1BP","GBP1","GBP2","GBP3","GBP4","GBP5","GBP6","GBP7","HLA-A","HLA-B","HLA-C","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-E","HLA-F","HLA-G","HLA-H","ICAM1","IFI30","IFNG","IFNGR1","IFNGR2","IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8","IRF9","JAK1","JAK2","MAPK1","MAPK3","MID1","MT2A","NCAM1","OAS1","OAS2","OAS3","OASL","PIAS1","PML","PRKCD","PTAFR","PTPN1","PTPN11","PTPN2","PTPN6","RAF1","SMAD7","SOCS1","SOCS3","SP100","STAT1","SUMO1","TRIM10","TRIM14","TRIM17","TRIM2","TRIM21","TRIM22","TRIM25","TRIM26","TRIM29","TRIM3","TRIM31","TRIM34","TRIM35","TRIM38","TRIM45","TRIM46","TRIM48","TRIM5","TRIM6","TRIM62","TRIM68","TRIM8","VCAM1","YBX1"
))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="T2ISGscore", slot="data", nbin=12)#nbin for monocyte

features <- list(c(
  "A1BG","ABCA13","ACAA1","ACLY","ACPP","ACTR10","ACTR1B","ACTR2","ADAM10","ADAM8","ADGRE3","ADGRE5","ADGRG3","AGA","AGL","AGPAT2","AHSG","ALAD","ALDH3B1","ALDOA","ALDOC","ALOX5","AMPD3","ANO6","ANPEP","ANXA2","AP1M1","AP2A2","APAF1","APEH","APRT","ARG1","ARHGAP9","ARL8A","ARMC8","ARPC5","ARSA","ARSB","ASAH1","ATAD3B","ATG7","ATP11A","ATP11B","ATP6V0A1","ATP6V0C","ATP6V1D","ATP8A1","ATP8B4","AZU1","B2M","B4GALT1","BIN2","BPI","BRI3","BST1","BST2","C16orf62","C1orf35","C3","C3AR1","C5AR1","C6orf120","CAB39","CALML5","CAMP","CAND1","CANT1","CAP1","CAPN1","CAT","CCT2","CCT8","CD14","CD177","CD300A","CD33","CD36","CD44","CD47","CD53","CD55","CD58","CD59","CD63","CD68","CD93","CDA","CDK13","CEACAM3","CEACAM6","CEP290","CFD","CHI3L1","CHIT1","CHRNB4","CKAP4","CLEC12A","CLEC4C","CLEC5A","CMTM6","CNN2","COMMD9","COPB1","COTL1","CPNE1","CPNE3","CPPED1","CR1","CRACR2A","CREG1","CRISP3","CRISPLD2","CST3","CSTB","CTSA","CTSB","CTSC","CTSD","CTSG","CTSH","CTSS","CTSZ","CXCL1","CXCR1","CXCR2","CYB5R3","CYBA","CYFIP1","CYSTM1","DBNL","DDOST","DEFA4","DEGS1","DERA","DGAT1","DIAPH1","DNAJC13","DNAJC3","DNAJC5","DOCK2","DOK3","DPP7","DSC1","DSG1","DSN1","DSP","DYNC1H1","DYNC1LI1","DYNLL1","DYNLT1","EEF1A1","EEF2","ELANE","ENPP4","EPX","ERP44","FABP5","FAF2","FCAR","FCER1G","FCGR2A","FCN1","FGR","FLG2","FOLR3","FPR1","FPR2","FRK","FTH1","FTL","FUCA1","FUCA2","GAA","GALNS","GCA","GDI2","GGH","GHDC","GLB1","GLIPR1","GM2A","GMFG","GNS","GOLGA7","GPI","GPR84","GRN","GSDMD","GSN","GSTP1","GUSB","GYG1","HEBP2","HEXB","HGSNAT","HK3","HMGB1","HMOX2","HP","HPSE","HRNR","HSP90AA1","HSP90AB1","HSPA8","HVCN1","IDH1","IGF2R","ILF2","IMPDH1","IMPDH2","IQGAP1","IQGAP2","IST1","ITGAL","ITGAM","ITGAV","ITGAX","ITGB2","JUP","KCMF1","KCNAB2","KPNB1","KRT1","LAIR1","LAMP1","LAMTOR2","LAMTOR3","LCN2","LGALS3","LPCAT1","LRG1","LRMP","LRRC7","LTA4H","LTF","LYZ","MAN2B1","MANBA","MAPK1","MAPK14","MCEMP1","METTL7A","MGAM","MGST1","MIF","MLEC","MME","MMP25","MMP8","MMP9","MNDA","MPO","MS4A3","MVP","NAPRT","NBEAL2","NCKAP1L","NCSTN","NFAM1","NFASC","NFKB1","NHLRC3","NIT2","NME1-NME2","NPC2","NRAS","OLFM4","ORMDL3","OSCAR","OSTF1","P2RX1","PA2G4","PADI2","PAFAH1B2","PDAP1","PDXK","PECAM1","PFKL","PGAM1","PGLYRP1","PGM1","PGM2","PIGR","PKM","PKP1","PLAC8","PLAUR","PLD1","PLEKHO2","PNP","PPBP","PPIA","PPIE","PRCP","PRDX6","PRG3","PRKCD","PRTN3","PSAP","PSEN1","PSMA2","PSMA5","PSMB1","PSMB7","PSMC2","PSMC3","PSMD1","PSMD11","PSMD12","PSMD13","PSMD14","PSMD2","PSMD3","PSMD6","PSMD7","PTAFR","PTGES2","PTPN6","PTPRB","PTPRC","PTPRJ","PTPRN2","PYCARD","PYGB","PYGL","QPCT","QSOX1","RAB10","RAB14","RAB18","RAB24","RAB27A","RAB31","RAB37","RAB3A","RAB3D","RAB4B","RAB5B","RAB5C","RAB6A","RAB7A","RAC1","RAP1A","RAP1B","RAP2B","RETN","RHOA","RHOF","RHOG","RNASE2","RNASE3","RNASET2","ROCK1","S100A11","S100A12","S100A8","S100A9","S100P","SCAMP1","SDCBP","SELL","SERPINA1","SERPINA3","SERPINB1","SERPINB12","SERPINB6","SIGLEC5","SIGLEC9","SIRPA","SIRPB1","SLC11A1","SLC15A4","SLC27A2","SLC2A3","SLC2A5","SLC44A2","SLCO4C1","SLPI","SNAP23","SNAP25","SNAP29","SPTAN1","SRP14","STBD1","STK10","STK11IP","STOM","SURF4","SVIP","TARM1","TBC1D10C","TCIRG1","TCN1","TIMP2","TLR2","TMC6","TMED7-TICAM2","TMEM173","TMEM179B","TMEM30A","TMEM63A","TNFAIP6","TNFRSF1B","TOLLIP","TOM1","TRAPPC1","TRPM2","TSPAN14","TTR","TUBB4B","TXNDC5","TYROBP","UBR4","UNC13D","VAMP8","VAPA","VAT1","VCL","VCP","VNN1","XRCC5","XRCC6","YPEL5"
))  
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="DEGscore", slot="data")

features <- list(c(
  "ABI1","ABI2","ABL1","ACTB","ACTG1","ACTR2","ACTR3","ARPC1A","ARPC1B","ARPC2","ARPC3","ARPC5","BAIAP2","BRK1","CD247","CD3G","CDC42","CRK","CYFIP1","CYFIP2","DOCK1","ELMO1","ELMO2","FCGR2A","FGR","FYN","GRB2","HCK","HSP90AA1","HSP90AB1","LIMK1","LYN","MAPK1","MAPK3","MYH2","MYH9","MYO10","MYO1C","MYO5A","MYO9B","NCK1","NCKAP1","NCKAP1L","NCKIPSD","NF2","PAK1","PIK3CA","PIK3CB","PIK3R1","PIK3R2","PLA2G6","PLCG1","PLCG2","PLD1","PLD2","PLD3","PLD4","PRKCD","PRKCE","PTK2","RAC1","SRC","SYK","VAV1","VAV2","VAV3","WASF1","WASF2","WASF3","WASL","WIPF1","WIPF2","WIPF3","YES1"
))  
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="Phascore", slot="data")

features <- list(c(
"MPO","NOS2","SLC11A1","ATP6V0A1","ATP6V1H","CYBA","NOS1","NCF4","ATP6V1D","ATP6V0A4","TCIRG1","ATP6V0E1","ATP6V1A","ATP6V1B1","NCF2","ATP6V0B","HVCN1","RAC2","ATP6V1F","ATP6V1E1","ATP6V1G1","ATP6V1C2","ATP6V1B2","ATP6V0D2","ATP6V1G3","ATP6V1C1","NCF1","ATP6V0D1","NOS3","CYBB","LPO","ATP6V0E2","ATP6V0A2","ATP6V0C","ATP6V1G2","ATP6V1G2","ATP6V1G2","ATP6V1G2","ATP6V1G2","ATP6V1G2","ATP6V1G2","ATP6V1E2","ATP6V1G3"
)) 
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="ROSscore", slot="data")


VlnPlot(s.int,c("IFNGscore1","DEGscore1","Phascore1","NETscore1", "ROSscore1", "NETosisscore1"), group.by="map",  ncol=1, pt.size=0)

plot_data <- FetchData(s.int, vars = c("DEGscore1", "map", "clusters"))


ggplot(plot_data, aes(x = clusters, y = DEGscore1, fill = map)) +
  geom_boxplot(width = 0.7, 
               position = position_dodge(width = 0.8),
               outlier.shape = NA,  # ???????????????????????????
               coef = 0.5) +  # ??????????????????????????????????????????????????????1.5???
  scale_fill_manual(values = rainbow(length(unique(plot_data$map)))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "clusters", y = "DEG", fill = "mao")+
  coord_cartesian(ylim = c(-0.08, 0.25))
  
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
cols = rev(brewer.pal(n = 11, name = "RdBu"))

DimPlot(s.int, group.by="clusters")
FeaturePlot(s.int, "ITGAM", split.by="disease")

pp <- DotPlot(s.int, features=rev(c("FCGR1A","FCGR3A","FCER1G","C3AR1")), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="clusters")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(1,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=0)+scale_y_discrete(limits = rev)+coord_flip()
pp

pp <- DotPlot(s.int, features=rev(c("C5AR1","C3AR1","C1QA", "C1QC")), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="combined")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(1,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=0)+scale_y_discrete(limits = rev(levels(factor(s.int$combined)))) +coord_flip()
pp

pp <- DotPlot(s.int, features=rev(c("T2ISGscore1")), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="combined")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(5,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=-0.25)+scale_y_discrete(limits = rev(levels(factor(s.int$combined))))+coord_flip()
pp
S
VlnPlot(s.int, "T2ISGscore1", pt.size=0,cols = c("#EC6464","#EC6464","#EC6464", "#30C030","#30C030","#30C030","#3030B0","#3030B0","#3030B0"))+
  geom_boxplot(width = 0.1)

pp <- DotPlot(s.int, features=c("C5AR1","C3AR1","C1QA","C1QC","CFB","MBL2","C3","C4A","C5"), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="clusters")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(1,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=0)+scale_y_discrete(limits = rev)
pp

pp <- DotPlot(s.int, features=rev(c("IFNG","IL10","IL6","IL1B","TNF","IL27","IL21")), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="seurat_clusters2")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(1,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=0)+scale_y_discrete(limits = rev)+coord_flip()
pp

pp <- DotPlot(s.sub, features=c("TNF"), scale.by="size", 
              scale=T,col.max=5, col.min=-5, group.by="combined")+
  theme(axis.text.x = element_text(angle=90,size =8,vjust=0.5,hjust=0.5))+ 
  theme(axis.text.y = element_text(size = 7.5), axis.title.x = element_blank())+scale_size(range=c(1,5))+
  scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=0)+coord_flip()
pp

pp <- dittoHeatmap(subset(s.int, downsample=1000),
             genes="MPO",
             annot.by=c("disease","sample"),
             order.by=c("disease","sample"),
             scale.to.max=T,
             slot="data",
             border_color="black")
pp

pp <- FeaturePlot(s.int, c("EGR1","IL1B","C3AR1","C5AR1"), max.cutoff="q95", split.by="CRP")
ggsave('Fig.1g.png', pp, width=10.5, height=35)
s.int$clusters_map <- paste(s.int$map, s.int$clusters, sep = "_")

DefaultAssay(s.int)
s.int <- DietSeurat(s.int, graphs="umap")
cds <- as.cell_data_set(s.int, assay="RNA")
cds <- cluster_cells(cds, reduction_method="UMAP", partition_qval=1e-10)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
cds <- learn_graph(cds, use_partition=T,close_loop=F, verbose=T, learn_graph_control=list(ncenter=500, prune_graph=T))
cds <- order_cells(cds)
cds <- choose_graph_segments(cds)

pp <- plot_cells(cds=cds, 
                 color_cells_by="clusters",
                 label_cell_groups = F,
                 label_groups_by_cluster=F,
                 show_trajectory_graph=T,
                 label_leaves=F,
                 label_branch_points=F,
                 label_roots = F,
                 label_principal_points=F,
                 graph_label_size = 2,
                 cell_size=1,
                 trajectory_graph_color="black")
pp
pp <- plot_cells(cds=cds, 
                 color_cells_by="pseudotime",
                 label_groups_by_cluster=F,
                 show_trajectory_graph=T,
                 label_leaves=F,
                 label_branch_points=F,
                 label_roots = F,
                 label_principal_points=F,
                 cell_size=1,
                 trajectory_graph_color="black")
pp

s.int <- AddMetaData(s.int, metadata=cds_sub_graph@principal_graph_aux@listData$UMAP$pseudotime,col.name="pseudotime")
FeaturePlot(s.int, "pseudotime", split.by="disease")
DefaultAssay(s.int) <- "RNA"

FeatureScatter(s.int, feature1="pseudotime2", feature2="PROM1")+  # ?????????
  geom_smooth(method = "loess", se = TRUE, color = "black", level=0.999)

Y <- c("HAVCR1","CDH6","TNFSF10","ANPEP","C4B","AMACR","C4A","VCAM1","SPON2","IL17RB","PLG","HPD","SLC22AB","BHMT","AGMAT","CUBN","PCK1","SLC22A12","SLC16A9","AFM","PFKFB3","PROM1","PAX8","MRC2","CA12","SIM2","EGR3","COL5A1","JAG1","HSD11B2")

# NA????????????????????????
plot_cells(cds_sub_graph,graph_label_size=10)

cds_sub_graph <- choose_graph_segments(cds)
cds_sub_graph <- choose_cells(cds, clear_cds = F, return_list = T)
cds_sub_graph <- cds[,cds_sub_graph]

clean_pseudotime <- pseudotime(cds_sub_graph)[!is.na(pseudotime(cds_sub_graph))]
clean_Y <- Y[!is.na(Y)]

# match???????????????NA?????????
match_result <- match(clean_Y, rownames(rowData(cds_sub_graph)))
valid_indices <- !is.na(match_result)

# NA?????????????????????????????????????????????????????????
pt.matrix <- exprs(cds_sub_graph)[match_result[valid_indices], order(clean_pseudotime)]
#pt.matrix <- exprs(cds)[match(Y,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# pseudotime????????????????????????????????????????????????
column_anno = HeatmapAnnotation(
  Pseudotime = anno_barplot(clean_pseudotime),
  show_annotation_name = TRUE
)

htkm2 <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 5,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  top_annotation               = column_anno  # ?????????????????????????????????
)

htkm2

Y <- c( "S100A9","S100A8","SELL","S100A12","VCAN",
        "EGR1","EGR2","DUSP6","ZFP36","CD74",
        "TXNIP","XIST","CPVL","CXCR4","FKBP5",
        "LINC00486","C1QC","C1QB","LYVE1","FCGBP",
        "C1QA","C3AR1","C5AR1","C3")
Y <- c( "S100A9","S100A8","SELL","S100A12",
        "EGR1","HLA.DPA1","HLA.DRA","HLA.DPB1",
        "TXNIP","XIST","CPVL","CXCR4",
        "LINC00486","C1QC","C1QB","LYVE1",
        "C3AR1","C5AR1","C3","VCAM1","SELL")

s.int <- RenameIdents(s.int, "CD4TNaive"="CD4_Naive","CD4TMemory"="CD4_TCM","CD4TEMRA"="CD4_TEM","Th1"="Th1","Th2"="Th2",   "Treg"="Treg", 
                      "CD8TNaive"="CD8_Naive","GZMK+CD8T"="CD8_GZMK", "GZMB+CD8T"="CD8_GZMB","dnT"="dnT","gdT"="gdT","MAIT"="MAIT",
                      "NK_CD56dim"="NK_CD56dim","NK_CD56high"="NK_CD56bright","NK_Proliferating"="NK_Prol",
                      "Tcell_Proliferating"="Tcell_Prol","Tcell_ISG"="Tcell_ISG",
                      "B_Naive"="B_Naive","B_Memory"="B_Memory","PBPC"="PBPC",
                      "CD14Mono"="CD14Mono","CD14Mono_Act"="CD14Mono_EGR","CD14Mono_HLA"="CD14Mono_HLA","CD14CD16Mono"="CD14CD16Mono","CD16Mono"="CD16Mono",
                      "cDC1"="cDC1","cDC2"="cDC2", "pDC"="pDC", "HSC"="HSC")
