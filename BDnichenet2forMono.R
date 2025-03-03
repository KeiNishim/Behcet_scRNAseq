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

Sys.setenv("http_proxy" = "http://170547:Kn01202584@www-proxy.chugai-pharm.co.jp:80")
Sys.setenv("https_proxy" = "http://170547:Kn01202584@www-proxy.chugai-pharm.co.jp:80")

library(nichenetr)
library(tidyverse)

Idents(s.int) <- "seurat_clusters2"
s.int <- RenameIdents(s.int,
                      "CD4 TRM"="CD4",
                      "Treg"="CD4",
                      "TFH"="CD4",
                      "CD4 CTL"="CD4",
                      "Plasma cell"="PC",
                      "ILC"="NK",
                      "PC_dim"="PC",
                      "TFH"="TFH",
                      "CD8 CTL"="CD8",
                      "CD8 TRM"="CD8",
                      "BAM"="Macrophage",
                      "PVM"="Macrophage",
                      "Monocyte"="Macrophage",
                      "Naive B"="Bcell",
                      "Memory B"="Bcell",
                      "cDC1"="cDC",
                      "cDC2"="cDC")

s.int <- RenameIdents(s.int,
                      "B"="Bcell","Mono"="Macrophage","Neu"="Neutrophil")

DimPlot(s.int)
s.int$seurat_clusters <- Idents(s.int)
DimPlot(s.int, group.by="disease")
DefaultAssay(s.int) <- "SCT"
DefaultAssay(s.sub) <- "SCT"
Ligand =c("CD4","CD8","PC","Bcell", "pDC","Macrophage","NK", "cDC","Neutrophil","Treg","Eosinophil")

Idents(s.int) <- "seurat_clusters3"
s.int <- subset(s.int, idents=Ligand)
s.sub <- subset(s.int, ident="Bcell")

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

expressed_genes_sender = rownames(s.int)
expressed_genes_receiver = rownames(FindMarkers(s.sub, ident.1="BD",ident.2="HD",only.pos = FALSE, min.pct=0.01,logfc.threshold = -Inf, max.cells.per.ident=2000)) #0.001 for Bcell
interested_genes_receiver <- rownames(top_n(
  FindMarkers(s.sub, ident.1="BD",ident.2="HD", logfc.threshold=0.2, min.pct=0.01, only.pos=T),n = -450, wt = p_val)#n=350 for CD4, 450 for Bcell
  )
geneset_oi = interested_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
ligands <-intersect(colnames(ligand_target_matrix),ligands)#
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_ligands
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
best_upstream_ligands = ligand_activities %>% top_n(40, auroc) %>% arrange(-auroc) %>% pull(test_ligand)
best_upstream_ligands

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(auroc) %>% arrange(auroc) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=auroc)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, auroc) %>% pull(pearson))), color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)
nrow(active_ligand_target_links_df)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#Circus output
Idents(s.int) <- "seurat_clusters"
avg_expression_ligands = AverageExpression(s.int, features = best_upstream_ligands, assay="RNA")
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd()*2)
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

CD4_specific_ligands = sender_ligand_assignment$CD4 %>% names() %>% setdiff(general_ligands)
CD8_specific_ligands = sender_ligand_assignment$CD8 %>% names() %>% setdiff(general_ligands)
Bcell_specific_ligands = sender_ligand_assignment$Bcell %>% names() %>% setdiff(general_ligands)
#PC_specific_ligands = sender_ligand_assignment$PC %>% names() %>% setdiff(general_ligands)
Neutrophil_specific_ligands = sender_ligand_assignment$Neutrophil %>% names() %>% setdiff(general_ligands)
Macrophage_specific_ligands = sender_ligand_assignment$Macrophage %>% names() %>% setdiff(general_ligands)
cDC_specific_ligands = sender_ligand_assignment$cDC %>% names() %>% setdiff(general_ligands)
#pDC_specific_ligands = sender_ligand_assignment$pDC %>% names() %>% setdiff(general_ligands)
Eosinophil_specific_ligands = sender_ligand_assignment$Eosinohpil %>% names() %>% setdiff(general_ligands)
Treg_specific_ligands = sender_ligand_assignment$Treg %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$NK %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(
                  rep("CD4-specific", times = CD4_specific_ligands %>% length()),
                  rep("Bcell-specific", times = Bcell_specific_ligands %>% length()),
                  #rep("PC-specific", times = PC_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("Macrophage-specific", times = Macrophage_specific_ligands %>% length()),
                  rep("cDC-specific", times = cDC_specific_ligands %>% length()),
                  rep("CD8-specific", times = CD8_specific_ligands %>% length()),
                  rep("Eosinophil-specific", times = Eosinophil_specific_ligands %>% length()),
                  rep("NK-specific", times = NK_specific_ligands %>% length()),
                  rep("Treg-specific", times = Treg_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  
  ligand = c(CD4_specific_ligands,
             Bcell_specific_ligands,Neutrophil_specific_ligands,
             Macrophage_specific_ligands, Treg_specific_ligands,
             cDC_specific_ligands, CD8_specific_ligands, Eosinophil_specific_ligands, NK_specific_ligands, 
             general_ligands
             ))


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "PC") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
#active_ligand_target_links_df <-active_ligand_target_links_df[-41,]
active_ligand_target_links_df

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.5)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "CD4-specific" = "darkgreen",
                   "Neutrophil-specific" = "steelblue4",
                   "CD8-specific" = "steelblue1",
                   "Bcell-specific"="orange4",
                   "Macrophage-specific"="violetred",
                   
                   "cDC-specific" = "steelblue2",
                   "Treg-specific" = "steelblue3",
                   "Eosinophil-specific"="orange2",
                   "NK-specific"="orange1"
                   
                   
                   
                   )
grid_col_target =c(
  "PC" = "tomato")
grid_col2 <- c(grid_col_ligand, grid_col_target)

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
library(dplyr)
links_circle = circos_links %>% select(ligand,target, weight)
ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
target_order = circos_links$target %>% unique()
ligand_order = c(CD4_specific_ligands,
                 Bcell_specific_ligands,Treg_specific_ligands, Neutrophil_specific_ligands,
                 Macrophage_specific_ligands, 
                 cDC_specific_ligands, CD8_specific_ligands, Eosinophil_specific_ligands, NK_specific_ligands, 
                 general_ligands
) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
width_same_cell_same_ligand_type = 2
width_different_cell = 1
width_ligand_target = 1
width_same_cell_same_target_type = 2
table(circos_links[,"ligand_type"])

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophage-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Bcell-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "cDC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Treg-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "PC") %>% distinct(target) %>% nrow() -1)),width_ligand_target
)
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, 
             link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, 
             diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.6)
  }, bg.border = NA)
legend("bottomright", legend = names(grid_col2), 
       fill = grid_col2, 
       title = "Cell Types", 
       cex = 0.8, 
       bty = "n")
#Save plot image

