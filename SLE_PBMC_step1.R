rm(list = ls())
setwd("/xywhd/SLE_PBMC")
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(patchwork)
library (readr)
library(data.table)
library(readr)
library(R.utils)
gunzip("GSE189050_final_seurat.RDS.gz", remove = 'TRUE')

scRNA <- readRDS("GSE189050_final_seurat.RDS")

dir.create("cluster")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
ggsave("cluster/VariableFeatures.pdf", plot = plot, width = 8, height = 6) 
ggsave("cluster/VariableFeatures.png", plot = plot, width = 8, height = 6)

##如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
##如果内存不够，可以只对高变基因进行标准化
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

#原始表达矩阵
GetAssayData(scRNA,slot="counts",assay="RNA") 
#标准化之后的表达矩阵              
GetAssayData(scRNA,slot="data",assay="RNA")
#中心化之后的表达矩阵 
GetAssayData(scRNA,slot="scale.data",assay="RNA") 

CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA))

g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(scRNA))
scRNA <- CellCycleScoring(object=scRNA,  g2m.features=g2m_genes,  s.features=s_genes)

scRNAa <- RunPCA(scRNA, features = c(s_genes, g2m_genes))
p <- DimPlot(scRNAa, reduction = "pca", group.by = "Phase")
ggsave("cluster/cellcycle_pca.png", p, width = 8, height = 6)

##如果需要消除细胞周期的影响
#scRNAb <- ScaleData(scRNA, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scRNA))

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
plotc <- plot1+plot2
ggsave("cluster/pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("cluster/pca.png", plot = plotc, width = 8, height = 4)

pc.num=1:20

scRNA <- FindNeighbors(scRNA, dims = pc.num) 
scRNA <- FindClusters(scRNA, resolution = 0.5)
table(scRNA@meta.data$seurat_clusters)
metadata <- scRNA@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster/cell_cluster.csv',row.names = F)

#tSNE
scRNA = RunTSNE(scRNA, dims = pc.num)
embed_tsne <- Embeddings(scRNA, 'tsne')
write.csv(embed_tsne,'cluster/embed_tsne.csv')
plot1 = DimPlot(scRNA, reduction = "tsne") 
ggsave("cluster/tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("cluster/tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNA <- RunUMAP(scRNA, dims = pc.num)
embed_umap <- Embeddings(scRNA, 'umap')
write.csv(embed_umap,'cluster/embed_umap.csv') 
plot2 = DimPlot(scRNA, reduction = "umap") 
ggsave("cluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("cluster/UMAP.png", plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("cluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("cluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)
##保存数据
saveRDS(scRNA, file="cluster/scRNA.rds")

p1 <- DimPlot(scRNA, reduction = "umap", group.by = "fine_cell_type", label = T) + NoLegend()
p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "fine_cell_type", label = T) + NoLegend()
p <- p1+p2
ggsave("tsne+umap_celltype.pdf", p, width = 24, height = 20)

scRNA$group<-recode(scRNA$classification,
                    "Control" = "Control",
                    "SLE ACT" = "SLE",
                    "SLE INACT" = "SLE")
saveRDS(scRNA, file="cluster/scRNA_group.rds")



rm(list = ls())
setwd("/xywhd/SLE_PBMC")
library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(patchwork)
library (readr)
library(data.table)
library(readr)
library(R.utils)
scRNA <- readRDS("cluster/scRNA_group.rds")

dir.create("cell_identify2")
#默认wilcox方法
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top5 = all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(all.markers, "cell_identify2/diff_genes_wilcox.csv", row.names = F)
write.csv(top5, "cell_identify2/top5_diff_genes_wilcox.csv", row.names = F)
##top5基因绘制热图
top5_genes <- read.csv("cell_identify2/top5_diff_genes_wilcox.csv")
top5_genes = CaseMatch(search = as.vector(top5_genes$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top5_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave("cell_identify2/top5_markers.pdf", plot=plot1, width=16, height=12) 
ggsave("cell_identify2/top5_markers.png", plot=plot1, width=16, height=12)

plot1 = DoHeatmap(scRNA, features = top5_genes, group.by = "fine_cell_type", group.bar = T, size = 4)
ggsave("cell_identify2/top5_cell_markers.pdf", plot=plot1, width=16, height=12) 
ggsave("cell_identify2/top5_cell_markers.png", plot=plot1, width=16, height=12)
saveRDS(scRNA,'cell_identify2/scRNA.rds')

dir.create("TNFSF10")
expr <- scRNA@assays$RNA$counts
expr<-as.data.frame(expr)
#Pecam1
gene_expression_TNFSF10 <- expr %>% .['TNFSF10',] %>% as.data.frame() %>% t()
gene_expression_TNFSF10 <- as.data.frame(gene_expression_TNFSF10)
colnames(gene_expression_TNFSF10) <- 'TNFSF10'
gene_expression_TNFSF10$cell <- rownames(gene_expression_TNFSF10)
#每个细胞表达TNFSF10的数值  
gene_expression_TNFSF10$TNFSF10
mean_value <- mean(gene_expression_TNFSF10$TNFSF10)
#mean:0.0266
gene_expression_TRAIL_high <- gene_expression_TNFSF10[which(gene_expression_TNFSF10$TNFSF10 >= 1),]
gene_expression_TRAIL_low <- gene_expression_TNFSF10[which(gene_expression_TNFSF10$TNFSF10 < 1),]
gene_expression_TRAIL_high_cells <- gene_expression_TNFSF10$cell[which(gene_expression_TNFSF10$TNFSF10 >= 1)]
gene_expression_TRAIL_low_cells <- gene_expression_TNFSF10$cell[which(gene_expression_TNFSF10$TNFSF10 < 1)]

write.table(gene_expression_TRAIL_high, file = "TRAIL_high.txt", sep = "\t", row.names = FALSE)
write.table(gene_expression_TRAIL_low, file = "TRAIL_low.txt", sep = "\t", row.names = FALSE)
data1 <- scRNA@meta.data
metadata <- data.frame(barcodes=rownames(data1), data1)

##按rbind()函数可以将数据框按行进行拼接，但它不会合并列
merged_df <- rbind(gene_expression_TRAIL_high, gene_expression_TRAIL_low)

metadata_2 <- merge(metadata,
                    merged_df,
                  by.x='barcodes',
                  by.y='cell',
                  all=T)
metadata_2$TNFSF10[metadata_2$TNFSF10 >= 1] <- 1
metadata_2$TNFSF10[metadata_2$TNFSF10 < 1] <- 0

#第一列作为行名
metadata_2 <- column_to_rownames(metadata_2,var = "barcodes")
sce_info <- AddMetaData(scRNA,metadata_2)
colnames(sce_info@meta.data)
sce_info$TNFSF10_group<-recode(sce_info$TNFSF10,
                    "0" = "low",
                    "1" = "high")
saveRDS(sce_info, file="scRNA_TNFSF10_group.rds")  
sce_info$celltype.group <- paste(sce_info$fine_cell_type, sce_info$TNFSF10_group, sep = "_")
saveRDS(sce_info, file="scRNA_TNFSF10_group.rds")  

rm(list = ls())
setwd("/xywhd/SLE_PBMC")
library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
dir.create("CellChat")
scRNA <- readRDS("scRNA_TNFSF10_group.rds")
#devtools::install_github("sqjin/CellChat")
Sys.setenv(RETICULATE_PYTHON="/home/xywh/anaconda3/bin/python3")
#要根据自己python3的路径来设置，可以在终端使用which python3来查看 

library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)

adole3k.final = scRNA
adole3k.final@commands$FindClusters

#构建Cell Chat对象时，输入的是log后的数据。
cellchat <- createCellChat(object=adole3k.final,group.by = "celltype.group")
#cellchat <- createCellChat(adole3k.final@assays$RNA@data, meta = adole3k.final@meta.data, group.by = "cell_type")
cellchat
summary(cellchat)
str(cellchat)
levels(cellchat@idents)
cellchat <- setIdent(cellchat, ident.use = "celltype.group")
groupSize <- as.numeric(table(cellchat@idents))
#查看每个cluster有多少个细胞，后面画图的时候需要用到这个值
groupSize

#导入配体受体数据库
CellChatDB <- CellChatDB.human
#导入小鼠是CellChatDB <- CellChatDB.mouse
str(CellChatDB) #查看数据库信息
#包含interaction、complex、cofactor和geneInfo这4个dataframe
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
#dev.new() #下一步不出图的时候运行
showDatabaseCategory(CellChatDB)

#在CellChat中，我们还可以先择特定的信息描述细胞间的相互作用，可以理解为从特定的侧面来刻画细胞间相互作用，比用一个大的配体库又精细了许多
unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object

library(future)
library(doFuture)
#预处理
## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 20)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- smoothData(cellchat, adj = PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

#推断细胞通讯网络

#推断配体-受体水平细胞通讯网络
#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
#对于某些函数，每个工作器都需要访问某些全局变量。如果这些变量大于默认限制，您将看到此错误。要解决这个问题，您可以设置
options(future.globals.maxSize = 10000 * 1024^2)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "CellChat/net_lr.csv")

#推断信号通路水平的细胞通讯网络
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "CellChat/net_pathway.csv")

#细胞互作关系展示
#所有细胞群总体观：细胞互作数量与强度统计分析：
#统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
# save as adult/net_number_strength.pdf

#检查每种细胞发出的信号
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as adult/net_number_individual.pdf

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as adult/net_strength_individual.pdf

mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
# save as adult/net_strength_individual.pdf

#单个信号通路或配体-受体介导的细胞互作可视化（层次图、网络图、和弦图、热图）
cellchat@netP$pathways  #查看都有哪些信号通路
pathways.show <- c("TRAIL")  

#层次图
levels(cellchat@idents)    # show all celltype
vertex.receiver = c(1,2,3,4,5,6,7,8,9,
                    10,12,13,14,15,16,17,18,19,
                    20,22,23,24,25,26,27,28) # define a numeric vector （淋系细胞）giving the index of the celltype as targets
#par(mar=c(5.1,4.1,4.1,2.1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# save as adult/TRAIL_hierarchy.pdf

#网络图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# save as adult/TRAIL_circle.pdf

#和弦图
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# save as adult/TRAIL_chord.pdf

#热图（Heatmap）
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# save as adult/TRAIL_heatmap.pdf

#配体-受体层级的可视化
#计算配体受体对选定信号通路的贡献值（在这里就是查看哪条信号通路对TGFb贡献最大）
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.TRAIL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) #提取对TGFb有贡献的所有配体受体 
# save as adult/TRAIL_LR_contribution.pdf

#层次图（ Hierarchy plot）
#提取对这个通路贡献最大的配体受体对来展示（也可以选择其他的配体受体对）
LR.show <- pairLR.TRAIL[1,] 
vertex.receiver = c(1,2,3,4,5,6,7,8,9,
                    10,12,13,14,15,16,17,18,19,
                    20,22,23,24,25,26,27,28) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# save as adult/TRAIL_hierarchy2.pdf

#网络图（Circle plot）
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# save as adult/TRAIL_circle2.pdf

#和弦图（Chord diagram）</p>
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
# save as adult/CXCL_chord2.pdf

#############################自动（批量）保存每个信号通路的互作结果</h5>
# Access all the signaling pathways showing significant communications将所有信号通路找出来
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
#vertex.receiver = c(1,2,4,6) #不画层次图就不需要这一步
dir.create("CellChat/all_pathways_com_circle") #创建文件夹保存批量画图结果
setwd("CellChat/all_pathways_com_circle")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], out.format = c("pdf"),
            vertex.receiver = vertex.receiver, layout = "circle") #绘制网络图
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
         plot=gg, width = 5, height = 2.5, units = 'in', dpi = 300)
}

setwd("/xywhd/SLE_PBMC")

##########7.5多个配体-受体介导的细胞互作关系可视化
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
p = netVisual_bubble(cellchat, sources.use = c(3,4,5,6), 
                     targets.use = c(1,2,7,8,9,
                                     10,12,13,14,15,16,17,18,19,
                                     20,22,23,24,25,26,27,28), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble.pdf", p, width = 8, height = 12) #髓系对淋巴的调节
# save as adult/Mye_Lymph_bubble.pdf  

#比如制定CCL和TRAIL这两个信号通路
netVisual_bubble(cellchat, sources.use = c(3,4,5,6), targets.use = c(1,2,7,8,9,
                                                                     10,12,13,14,15,16,17,18,19,
                                                                     20,22,23,24,25,26,27,28), 
                 signaling = c("CCL","TRAIL"), remove.isolate = FALSE)  

# show all the significant interactions (L-R pairs) based on user's input
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","TRAIL","TGFb"))
netVisual_bubble(cellchat, sources.use = c(3,4,5,6), targets.use = c(1,2,7,8,9,
                                                                     10,12,13,14,15,16,17,18,19,
                                                                     20,22,23,24,25,26,27,28), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE)   

#参与某条信号通路（如TRAIL）的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
## Plot the signaling gene expression distribution
p = plotGeneExpression(cellchat, signaling = "TRAIL")
ggsave("TRAIL_GeneExpression_vln.pdf", p, width = 8, height = 8)
p = plotGeneExpression(cellchat, signaling = "TRAIL", type = "dot")
ggsave("TRAIL_GeneExpression_dot.pdf", p, width = 8, height = 6) 

#通讯网络系统分析
#计算网络中心性权重
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  
#通过计算每个细胞群的网络中心性指标，识别每类细胞在信号通路中的角色/作用C（发送者、接收者、调解者和影响者）
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 15, height = 6, font.size = 10)
# # save as adult/SNA_TRAIL_signalingRole.pdf 

#识别细胞的信号流模式 
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 4)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 4)
ht1 + ht2
# save as adult/SNA_SignalingPattern.pdf  

#非负矩阵分解（NMF）识别细胞的通讯模式    
#信号输出细胞的模式识别
#计算分解成几个因子(pattern)比较合适（这一步运行比较慢 。在使用NMF对细胞进行亚群细分时，如果不测试的话，最好选择比细胞类型多一点的值）
selectK(cellchat, pattern = "outgoing")
# save as adult/NMF_outgoing_selectK.pdf

nPatterns = 3 # 挑选曲线中第一个出现下降的点（从2就开始下降了）
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)
# save as adult/NMF_outgoing_comPattern_heatmap.pdf

#冲击图/河流图（river plot）</p>
netAnalysis_river(cellchat, pattern = "outgoing")
# save as adult/NMF_outgoing_comPattern_river.pdf

#气泡图</p>
netAnalysis_dot(cellchat, pattern = "outgoing")

#######信号输入细胞的模式识别</li>
selectK(cellchat, pattern = "incoming") 
# save as adult/NMF_incoming_selectK.pdf

#热图</p>
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, 
                                          width = 5, height = 9, font.size = 6)
# save as adult/NMF_incoming_comPattern_heatmap.pdf

#冲击图</p>
netAnalysis_river(cellchat, pattern = "incoming")
# save as adult/NMF_incoming_comPattern_river.pdf

#气泡图</p>
netAnalysis_dot(cellchat, pattern = "incoming")
# save as adult/NMF_incoming_comPattern_dotplot.pdf

############信号网络的流行学习与分类
#基于功能相似性的流行学习与分类</li>
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#Error in do_one(nmeth) : NA/NaN/Inf in foreign function call (arg 1)
p = netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
ggsave("Manifold_functional_cluster.pdf", p, width = 8, height = 6)
#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#基于拓扑相似性的流行学习与分类
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
p = netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
ggsave("Manifold_structural_cluster.pdf", p, width = 8, height = 6)

## The end
saveRDS(cellchat, file = "cellchat.rds")

##########################不同分组之间的配对分析
#数据准备，分别创建CellChat对象
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
rm(list = ls())
options(stringsAsFactors = FALSE)

## 创建cellchat对象
### 提取数据子集
setwd("/xywhd/SLE_PBMC")
scRNA <- readRDS("scRNA_TNFSF10_group.rds")
table(scRNA$celltype.group)
Idents(scRNA) <- 'celltype'
#scRNA <- subset(scRNA, idents = c('Monocytes','Endothelial cells','T cells','Obsteoblastic','Fibroblasts','Dendritic cells','B cells'
                                  ,'Osteoclasts','Mast cells'))
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
table(scRNA$group)
Idents(scRNA) <- 'group'
adu <- subset(scRNA, idents = c('SLE'))########实验组
ado <- subset(scRNA, idents = c('Control'))


#创建cellchat对象
cco.adu <- createCellChat(adu@assays$RNA@data, meta = adu@meta.data, group.by = "celltype.group")
cco.ado <- createCellChat(ado@assays$RNA@data, meta = ado@meta.data, group.by = "celltype.group")
setwd("/xywhd/SLE_PBMC/CellChat")
dir.create("vs")
setwd("./vs")
save(cco.adu, cco.ado, file = "cco.rda")
# load("../cco.rda")
# cco.ado <- setIdent(cco.ado, ident.use = "celltype")
# cco.adu <- setIdent(cco.adu, ident.use = "celltype")

#对于某些函数，每个工作器都需要访问某些全局变量。如果这些变量大于默认限制，您将看到此错误。要解决这个问题，您可以设置
options(future.globals.maxSize = 10000 * 1024^2)
cellchat <- cco.ado
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.ado <- cellchat
saveRDS(cco.ado, "cco.ado.rds")

#########分析样本cco.adu的细胞通讯网络
cellchat <- cco.adu
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
cco.adu <- cellchat
saveRDS(cco.adu, "cco.adu.rds")

#合并cellchat对象
cco.aduist <- list(adole=cco.ado, adult=cco.adu)
cellchat <- mergeCellChat(cco.aduist, add.names = names(cco.aduist), cell.prefix = TRUE)

#3. 可视化</h3>
#3.1 所有细胞群总体观：通讯数量与强度对比</h5>
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

#数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# save as Diff_number_strength_net.pdf

#数量与强度差异热图</p>红色是case相对于control上调的，蓝色是下调的
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
# save as Diff_number_strength_heatmap.pdf

#细胞互作数量对比网络图</p>
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.aduist, attribute = c("idents","count"))
for (i in 1:length(cco.aduist)) {
  netVisual_circle(cco.aduist[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.aduist)[i]))
}
# save as Counts_Compare_net.pdf

#########3.2 指定细胞互作数量对比网络图</h5>
par(mfrow = c(1,2))
s.cell <- c('ABCs_high','ABCs_low','T cells','CD14+ monocytes_high','CD14+ monocytes_low','CD16+ monocytes_high','CD16+ monocytes_low',
            'CD4+ T cells_high','CD4+ T cells_low', 'CD8+ cytotoxic T cells_high', 'CD8+ cytotoxic T cells_low', 'CD8+ T cells_high', 'CD8+ T cells_low',
            'cDCs_high', 'cDCs_low', 'Memory B cells_high', 'Memory B cells_low', 'Naive B cells_high', 'Naive B cells_low', 'NK cells_high', 'NK cells_low',
            'pDCs_high', 'pDCs_low', 'Plasmablasts_high', 'Plasmablasts_low', 'Progenitors_high', 'Progenitors_low', 'Transitional B cells_high', 'Transitional B cells_low')
count1 <- cco.aduist[[1]]@net$count[s.cell, s.cell]
count2 <- cco.aduist[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.aduist)[1]))
netVisual_circle(count2, weight.scale = T, label.edge= T, edge.weight.max = weight.max, edge.width.max = 12, 
                 title.name = paste0("Number of interactions-", names(cco.aduist)[2]))
# save as Counts_Compare_select.pdf 10*6.5

####3.3 保守和特异性信号通路的识别与可视化</h5>
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 15, height = 12)   

#####3.4 流行学习识别差异信号通路
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot', type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
p <- rankSimilarity(cellchat, type = "structural") + ggtitle("Structural similarity of pathway")
ggsave("Pathway_Similarity.pdf", p, width = 8, height = 5)

saveRDS(cellchat, "cellchat.rds") 

##3.5 细胞信号模式对比</h5>
library(ComplexHeatmap)
pathway.union <- union(cco.aduist[[1]]@netP$pathways, cco.aduist[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.aduist[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.aduist)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.aduist[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.aduist)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_all.pdf  10*6 

##输出信号模式对比</p>
pathway.union <- union(cco.aduist[[1]]@netP$pathways, cco.aduist[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.aduist[[1]], pattern = "outgoing", signaling = pathway.union, 
                                        title = names(cco.aduist)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.aduist[[2]], pattern = "outgoing", signaling = pathway.union,
                                        title = names(cco.aduist)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_outgoing.pdf  10*6 

#输入信号模式对比</p>
pathway.union <- union(cco.aduist[[1]]@netP$pathways, cco.aduist[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.aduist[[1]], pattern = "incoming", signaling = pathway.union, 
                                        title = names(cco.aduist)[1], width = 8, height = 22)
ht2 = netAnalysis_signalingRole_heatmap(cco.aduist[[2]], pattern = "incoming", signaling = pathway.union,
                                        title = names(cco.aduist)[2], width = 8, height = 22)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# save as Compare_signal_pattern_incoming.pdf  10*6

#特定信号通路的对比</h5>
#normal存在，疾病不存在
#CD40
#CD86
#LIGHT
#NECTIN
#TIGIT
#SEMA4
#网络图</p>
pathways.show <- c("IL16") 
weight.max <- getMaxWeight(cco.aduist, slot.name = c("netP"), attribute = pathways.show) 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.aduist)) {
  netVisual_aggregate(cco.aduist[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.aduist)[i]))
}
# save as Compare_IL16_net.pdf  10*6.5

#热图</p>
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(cco.aduist)) {
  ht[[i]] <- netVisual_heatmap(cco.aduist[[i]], signaling = pathways.show, color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(cco.aduist)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# save as Compare_IL16_heatmap.pdf  12*6.5

#和弦图</p>
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.aduist)) {
  netVisual_aggregate(cco.aduist[[i]], signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05,
                      vertex.label.cex = 0.6, signaling.name = paste(pathways.show, names(cco.aduist)[i]))
}
# save as Compare_IL16_chord.pdf  10*6.5

#######配体-受体对比分析 
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3,4,5,6), targets.use = c(1,2,7,8,9,
                                                                          10,12,13,14,15,16,17,18,19,
                                                                          20,22,23,24,25,26,27,28),  comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble.pdf", p, width = 40, height = 40)

#气泡图展示上调或下调的配体受体对</p>
p1 <- netVisual_bubble(cellchat, sources.use = c(3,4,5,6), targets.use = c(1,2,7,8,9,
                                                                           10,12,13,14,15,16,17,18,19,
                                                                           20,22,23,24,25,26,27,28), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in adult", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(3,4,5,6), targets.use = c(1,2,7,8,9,
                                                                                 10,12,13,14,15,16,17,18,19,
                                                                                 20,22,23,24,25,26,27,28), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in adult", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", pc, width = 40, height = 40)

#和弦图</p>
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(cco.aduist)) {
  netVisual_chord_gene(cco.aduist[[i]], sources.use = c(3,4,5,6), targets.use = c(1,2,7,8,9,
                                                                                  10,12,13,14,15,16,17,18,19,
                                                                                  20,22,23,24,25,26,27,28), signaling = "MIF", 
                       lab.cex = 0.6, legend.pos.x = 10, legend.pos.y = 20,
                       title.name = paste0("Signaling from Treg - ", names(cco.aduist)[i]))
}
# save as Compare_LR_MIF_chord.pdf  10*6.5















