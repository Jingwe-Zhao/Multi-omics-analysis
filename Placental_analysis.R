library(Seurat)
library(tidyverse)

input_data_dir <- 'D:/DATA/XMU/bioinfor/data/scrna'
dir <- list.files(input_data_dir)
dir

names(dir) <-dir

scRNAlist <- list()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i],gene.column = 2)
  scRNAlist[[i]] <- CreateSeuratObject(counts,  min.cells=3, min.features = 200)
}


scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:4])
scRNA<-JoinLayers(scRNA)
table(scRNA@meta.data$orig.ident)

scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA","mt_percent"), ncol = 3,raster=FALSE)


maxGene=4000
mt=2
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA  < maxGene & mt_percent < mt)



###Figure1-------
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA","mt_percent"),group.by = 'orig.ident', ncol = 3,pt.size = 0,raster=FALSE)

scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)
DimPlot(scRNA,reduction = 'pca',group.by = "orig.ident")

library(harmony)

scRNA <- RunHarmony(scRNA,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30,reduction.name = "umap")

DimPlot(scRNA, reduction = "umap",group.by = "orig.ident")


scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.6) 
DimPlot(scRNA, reduction = "umap",group.by = "seurat_clusters",label = T)


DotPlot(scRNA,features=c('COL1A1','MYH11','CD14','CSF1R','PSG5','PSG8','TP63','PEG10','VIM','VWF','HLA-G','ASCL2'),cols = c('gray','#1874CD')) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )



scRNA@meta.data$Celltype <- NA
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('0','1','2','3','6','7','9','11')] <- "lSTB"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters  %in% c('4','14')] <- "eSTB"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('8')] <- "EC"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('10')] <- "STR"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('13')] <- "Mac"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('12')] <- "EVT"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('5')] <- "eCTB"
scRNA@meta.data$Celltype[scRNA@meta.data$seurat_clusters %in% c('15')] <- "lCTB"

DimPlot(scRNA, reduction = "umap",group.by = "Celltype",label = T)


My_levels <- c('EVT','EC','lCTB','eCTB','lSTB','eSTB','Mac','STR')
scRNA$Celltype <- factor(scRNA$Celltype, levels = My_levels)
Idents(scRNA)<-scRNA$Celltype

DotPlot(scRNA,features=c('COL1A1','MYH11','CD14','CSF1R','PSG5','PSG8','TP63','PEG10','VIM','VWF','HLA-G','ASCL2'),cols = c('gray','#1874CD')) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )



###Figure2-------
DEG<-FindAllMarkers(scRNA,only.pos = T,min.pct = 0.25)

jjVolcano(diffData = DEG,topGeneN=3,aesCol = c('purple','#CD8C95',size=2,fontface = 'arial'))

head(DEG)

##GO-KEGG
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
gid <-bitr(unique(DEG$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')  
Genelist <- full_join(DEG, gid, by=c('gene' = 'SYMBOL'))
Genelist<-subset(Genelist,Genelist$p_val_adj < 0.001 & Genelist$avg_log2FC > 2)
table(Genelist$cluster)


head(Genelist)

cell_count <- c(EVT=202, EC=281, lCTB=230, eCTB=122, lSTB=63, eSTB=102, Mac=150, STR=188)
df <- data.frame(
  cell_type = names(cell_count),
  gene_num = as.integer(cell_count)
)

df$cell_type <- factor(df$cell_type, levels = df$cell_type[order(df$gene_num)])


ggplot(df, aes(x = gene_num, y = cell_type)) +
  geom_col(fill = "#4682B4", width = 0.7) +
  geom_text(aes(label = gene_num), hjust = -0.3, size = 3.8) +
  labs(
    x = "Number of differential genes",
    y = "Cell Type",
    title = "DEGs ( avg_log2FC > 2, p_val_adj < 0.001 )"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

##GO BP
GOBP = compareCluster(ENTREZID~ cluster,data = Genelist, fun='enrichGO',OrgDb ='org.Hs.eg.db',ont=BP)

dotplot(GOBP, label_format=40,showCategory=3) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle('GO Biological processes')  +
  scale_fill_distiller(palette = "Blues", direction = -1)

##GO CC
GOCC = compareCluster(ENTREZID~ cluster,data = Genelist, fun='enrichGO',OrgDb ='org.Hs.eg.db',ont='CC')

dotplot(GOCC, label_format=40,showCategory=3) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle('GO Cellular components')  +
  scale_fill_distiller(palette = "Blues", direction = -1)

##GO MF
GOMF = compareCluster(ENTREZID~ cluster,data = Genelist, fun='enrichGO',OrgDb ='org.Hs.eg.db',ont='MF')

dotplot(GOMF, label_format=40,showCategory=3) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle('GO Molecular function')  +
  scale_fill_distiller(palette = "Blues", direction = -1)


##KEGG
kegg= compareCluster(ENTREZID~ cluster,data = Genelist, fun='enrichKEGG', organism="hsa")
df<-kegg@compareClusterResult

table(df$Cluster)
head(df)

library(tidyverse)

df_sorted <- df %>%
  group_by(Cluster) %>%
  arrange(desc(zScore), .by_group = TRUE) %>%
  mutate(rank_in_cl = row_number()) %>%
  ungroup()


used_path <- c()
plot_df <- tibble()
cl_order <- c("EVT","EC","lCTB","eCTB","lSTB","eSTB","Mac","STR")

for (cl in cl_order) {
  sub <- df_sorted %>% filter(Cluster == cl)
  avail <- sub %>% filter(!Description %in% used_path)
  pick <- slice_head(avail, n = 3)
  used_path <- c(used_path, pick$Description)
  plot_df <- bind_rows(plot_df, pick)
}


plot_df$Description <- str_wrap(plot_df$Description, width = 50)

ggplot(plot_df, aes(x = Cluster, fill = Cluster)) +
  geom_col(
    aes(y = zScore, group = interaction(Cluster, Description)),
    position = position_dodge(width = 0.9),
    width = 0.75
  ) +
  geom_text(
    aes(
      label = Description,
      group = interaction(Cluster, Description),
      y = 0 
    ),
    position = position_dodge(width = 0.9),
    size = 2.2, 
    color = "black", 
    fontface = "bold",
    angle = 90,
    vjust = 0, 
    hjust = 0
  ) +
  labs(
    x = "Cell Cluster",
    y = "zScore",
    title = "Top3 Unique KEGG Pathways for Each Cell Cluster",
    fill = "Cell Cluster"
  ) +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )



###Figure3-------
input_data_dir <- 'D:/DATA/XMU/bioinfor/data/st'
sample_list <- list.files(input_data_dir)
sample_list


## file path
samples_dir <- sample_list %>% file.path(input_data_dir, .)
samples_dir

##slice id
sample_names <- c("PLA1","PLA2", "PLA3", "PLA4")

library(Seurat)
library(tidyverse)

sample_objects <- purrr::map(1:length(sample_list), function(x) {
  ## read data
  one_dir <- samples_dir[x]
  sample_id <- sample_list[x]
  slice_id <- sample_names[x]
  sample_object <- Load10X_Spatial(
    data.dir = one_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = slice_id,
    filter.matrix = TRUE
  )
  sample_object@project.name <- sample_id
  sample_object@meta.data$orig.ident <- slice_id
  sample_object <- RenameCells(object = sample_object, add.cell.id = slice_id)
  
  return(sample_object)
})



#SCT transform
sample_objects <- lapply(sample_objects, 
                         SCTransform, 
                         assay = "Spatial", 
                         method = "poisson")

ST<-merge(sample_objects[[1]], y = sample_objects[2:4])

DefaultAssay(ST) <- "SCT"

VariableFeatures(ST) <- c(VariableFeatures(sample_objects[[1]]), 
                          VariableFeatures(sample_objects[[2]]),
                          VariableFeatures(sample_objects[[3]]), 
                          VariableFeatures(sample_objects[[4]]))



ST<-RunPCA(ST, assay = "SCT", verbose = FALSE)


library(harmony)
table(ST$orig.ident)

ST <- RunHarmony(ST,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")

#FindClusters
ST <- FindNeighbors(ST, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4)
ST <- RunUMAP(ST, reduction="harmony",dims = 1:30)



library(RColorBrewer)
display.brewer.all()
cols <- colorRampPalette(brewer.pal(9, 'Set1'))(15)

DimPlot(ST, reduction = "umap",cols = cols)+ggtitle('ST Clusters')+
  theme(plot.title = element_text(hjust = 0.5))

SpatialPlot(ST,group.by ='seurat_clusters',ncol = 2,label = F,pt.size.factor = 2.8)& NoLegend() & scale_fill_manual(values = cols)

##RCTD
library(CARD)

sc_count <- GetAssayData(scRNA,layer = 'counts')
sc_meta <- scRNA@meta.data %>% 
  rownames_to_column("cellID") %>%
  dplyr::select(cellID,orig.ident,Celltype) %>% 
  mutate(CB = cellID) %>% 
  column_to_rownames("CB")
head(sc_meta)


CARD_results <- list()
ST_list<-SplitObject(ST,split.by = 'orig.ident')
for(i in 1:4) {
  cat("decov", i, "sample...\n")
  
  spatial_count <- GetAssayData(ST_list[[i]], layer = 'counts')
  cat("Expression matrix dimension:", dim(spatial_count), "\n")
  
  spatial_loca <- GetTissueCoordinates(ST_list[[i]])
  spatial_loca <- spatial_loca[1:2]
  cat("spatial local:", dim(spatial_loca), "\n")
  
  # creat CARD_obj
  CARD_obj <- createCARDObject( 
    sc_count = sc_count, 
    sc_meta = sc_meta, 
    spatial_count = spatial_count, 
    spatial_location = spatial_loca, 
    ct.varname = "Celltype", 
    ct.select = unique(sc_meta$celltype),
    sample.varname = "orig.ident",
    minCountGene = 100,
    minCountSpot = 5)
  
  CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
  
  CARD_results[[i]] <- CARD_obj
}

names(CARD_results) <- c("PLA1","PLA2", "PLA3", "PLA4")


for (i in 1:4) {
  message(paste("draw CARD_results [[", i, "]]...", sep = ""))

  CARD.visualize.pie(
    proportion       = CARD_results[[i]]@Proportion_CARD,
    spatial_location = CARD_results[[i]]@spatial_location
  )
}


P1<-CARD.visualize.Cor(CARD_results[[1]]@Proportion_CARD) 
P2<-CARD.visualize.Cor(CARD_results[[2]]@Proportion_CARD) 
P3<-CARD.visualize.Cor(CARD_results[[3]]@Proportion_CARD) 
P4<-CARD.visualize.Cor(CARD_results[[4]]@Proportion_CARD) 

library(ggpubr)
ggarrange(P1,P2,P3,P4,ncol = 4)



ct.visualize = c("EVT","EC") 

for (i in 1:4) {
  message(paste("draw CARD results [[", i, "]]...", sep = ""))

  current_prop <- CARD_results[[i]]@Proportion_CARD
  current_pos  <- CARD_results[[i]]@spatial_location

  p <- CARD.visualize.prop(
    proportion       = current_prop,        
    spatial_location = current_pos,   
    ct.visualize     = ct.visualize,                  
    colors           = c("lightblue", "lightyellow", "red"),   
    NumCols          = 2,
    pointSize        = 1

}






###Figure4-------
library(CellChat)
data.input <-GetAssayData(scRNA,layer = 'data')
labels <- scRNA$Celltype
meta <-scRNA@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Celltype")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use


# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



cellchat@netP$pathways
netVisual_heatmap(cellchat, signaling = 'PDGF', color.heatmap = "Reds")

netVisual_chord_gene(cellchat, sources.use = c(7), targets.use = 1, legend.pos.x = 15)


netVisual_aggregate(cellchat, signaling = 'SPP1', layout = "chord")

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = 'SPP1', geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair

netVisual_individual(cellchat, signaling = 'SPP1', pairLR.use = LR.show, layout = "chord")


genes_to_plot <- c("SPP1", "ITGA5", "ITGB1")

for (i in 1:4) {
  message(paste("draw CARD_results [[", i, "]] ...", sep = ""))
  
  current_exp <- CARD_results[[i]]@spatial_countMat
  current_pos <- CARD_results[[i]]@spatial_location
  
  CARD.visualize.gene(
    spatial_expression = current_exp,  
    spatial_location   = current_pos,  
    gene.visualize     = genes_to_plot,  
    colors             = NULL, 
    NumCols            = 3
  )
}
