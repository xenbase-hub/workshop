# TOOLBOX
# By Konrad Thorner
# A collection of commonly used analysis steps in R

# RNA-SEQ
#------------------------------------------------------------------------------------------

# Takes a raw counts matrix, a vector of annotations, and a vector of comparisons written as "x vs y"
# Order of annotations must match order of samples in counts table
# Outputs 
run_deseq <- function(counts=NULL,annotations=NULL,comparisons=NULL,fold_change=2.0,p_value=0.05,top_n=NULL,save=F) {
	require(DESeq2)
	data <- round(counts)
	sample_table <- data.frame(annotations)
	rownames(sample_table) <- colnames(counts)
	colnames(sample_table) <- "Group"
	print("IMPORTANT: Check that all annotations are correct!")
	for (x in 1:length(annotations)) {
		print(paste0(rownames(sample_table)[x]," - ",annotations[x]))
	}
	sample_table$Group <- as.factor(sample_table$Group)
	dds <- DESeqDataSetFromMatrix(countData = data,
	                              colData = sample_table,
	                              design = ~ Group)
	dds <- DESeq(dds)
	all_comparisons <- list()
	all_names <- c()
	for (comparison in comparisons) {
		groups <- unlist(strsplit(comparison, " "))
		contrast = c("Group",groups[1],groups[3])
		res <- results(dds,contrast = contrast)
		res <- res[!is.na(res$log2FoldChange),]
		res <- res[abs(res$log2FoldChange) >= fold_change,]
		res <- res[!is.na(res$padj),]
		res <- res[res$padj <= p_value,]
		resOrdered <- res[order(res$padj),]
		resOrdered <- as.data.frame(resOrdered)
		print(paste0(nrow(resOrdered)," total DEGs found for ", comparison))
		if (!is.null(top_n)) {
			resOrdered <- resOrdered[1:min(nrow(resOrdered),top_n),]
		}
		if (save) {
			write.csv(resOrdered, file=paste0(groups[1],"_vs_",groups[3],".csv"),quote = F)
		}
		all_comparisons <- append(all_comparisons,list(resOrdered))
		all_names <- c(all_names,comparison)
	}
	names(all_comparisons) <- all_names
	if (length(all_comparisons) == 1) {
		return(all_comparisons[[1]])
	} else {
		return(all_comparisons)
	}
}

# Takes a counts matrix
# Create a heatmap plot to the current directory
create_heatmap <- function(counts=NULL,genes=NULL,show_names=T,save=F) {
	require(pheatmap)
	data_use <- counts
	sample_names <- colnames(counts)
	sampleTable <- data.frame(name = factor(sample_names))
	rownames(sampleTable) <- sample_names
	dds <- DESeqDataSetFromMatrix(counts, sampleTable, ~1)
	vsd <- vst(dds, blind=TRUE)
	heatmap <- pheatmap(assay(vsd)[genes,], scale="row", cluster_rows=T, cluster_cols=T, main="Heatmap", border_color = "NA",show_rownames = show_names, color=colorRampPalette(c("navy", "white", "red"))(50))
	if (save) {
		png("heatmap.png", width = 850, height = 1000)
		print(heatmap, useSource=TRUE)
		dev.off()
	}
	return(heatmap)
}

# Takes a raw counts matrix
# Saves a PCA plot to the current directory
create_pca <- function(counts=NULL,name="PCA",save=F) {
	require(DESeq2)
	require(ggplot2)
	sample_names <- colnames(counts)
	sampleTable <- data.frame(name = factor(sample_names))
	rownames(sampleTable) <- sample_names
	dds <- DESeqDataSetFromMatrix(counts, sampleTable, ~1)
	vsd <- vst(dds, blind=TRUE)
	pca <- plotPCA(vsd, intgroup=c("name"))
	if (save) {
		pca
		ggsave(paste0(name,".png"),width=8,height=10,limitsize = FALSE)
	}
	return(pca)
}

# Takes a TPM matrix, and an option to include a legend
# Creates a PCA plot and optionally saves to the current directory
create_pca_2 <- function(tpm=NULL,legend=T, name="PCA", save=F) {
	require(plotly)
	require(dplyr)
	require(ggfortify)
	require(ggplot2)
	data_use <- as.matrix(tpm)
	data_use <- data_use[apply(data_use[,-1], 1, function(x) !all(x==0)),]
	data_use <- log2(data_use+1)
	data_use.t <- t(data_use)
	data_use.t <- data_use.t[ , which(apply(data_use.t, 2, var) != 0)]
	pca <- prcomp(data_use.t, center=T, scale. = T)
	pc1 <- round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
	pc2 <- round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
	PC1_use <- paste0("PC1", "(", pc1, "%)")
	PC2_use <- paste0("PC2", "(", pc2, "%)")
	Samples_temp <- rownames(data_use.t)
	Samples <- factor(Samples_temp)
	scores <- data.frame(Samples_temp, pca$x[,1:2])
	MIN_X <- min(scores$PC1)
	Max_X <- max(scores$PC1)
	header <- "PCA"
	if (legend) {
		plot <- qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=16,size=4) + geom_text(aes(label=Samples_temp,size=0), hjust=-0.2, vjust=0.35) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.title=element_text(size=14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header) + theme(legend.position="none") + scale_colour_manual(values = cols)
	} else {
		plot <- qplot(x=PC1, y=PC2, data=scores, colour=Samples, xlim=c(MIN_X-75,Max_X+75)) + xlab(PC1_use) + ylab(PC2_use) + geom_point(shape=16,size=4) + scale_size_area() + theme(axis.text = element_text(size = 14),axis.title=element_text(size=14),axis.line.x = element_line(colour = "black"),axis.line.y = element_line(colour = "black"),legend.key = element_rect(fill = "white"),legend.background = element_rect(fill = "white"),panel.grid.major = element_line(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white")) + ggtitle(header) + scale_colour_manual(values = cols)
	}
	if (save) {
		plot
		ggsave(paste0(name,".png"),width=10,height=10)
	}
	return(plot)
}

# Takes a TPM matrix, a table of markers subset to include gene, annotation, and fold change columns, and an optional name
# Creates a table of estimated cell proprotions
estimate_cells <- function(tpm=NULL,markers=NULL,name="bisque",save=F) {
	require(BisqueRNA)
	require(Biobase)
	require(Seurat)
	require(SeuratObject)
	require(dplyr)
	bulk.eset <- Biobase::ExpressionSet(assayData = as.matrix(tpm))
	res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, data.frame(markers), weighted=T, min_gene = 1)
	bisque <- data.frame(res$bulk.props)
	if (save) {
		save_table(bisque,name)
	}
	return(bisque)
}

# SINGLE CELL
#------------------------------------------------------------------------------------------

# Take a single file location. Must be in .rds format
# Returns a Seurat object
load_seurat <- function(dir=NULL) {
	obj <- readRDS(dir)
	print(paste0(get_seurat_basename(dir), " | ", ncol(obj), " cells, ", nrow(obj), " genes, ~ ", format(object.size(obj),"Gb")))
	return(obj)
}

# Accepts either a directory or a vector of file locations. Must be in .rds format
# Returns a list of named Seurat objects
load_multiple_seurat <- function(dir=NULL,file_list=NULL) {
	if (is.null(file_list)) {
		seurat_files <- list.files(path=dir,pattern=".rds$",recursive=FALSE,full.names=TRUE)
	} else if (is.null(dir)) {
		seurat_files <- file_list
	}
	seurat.objects <- c()
	seurat.names <- c()
	num_files <- length(seurat_files)
	print(paste0("Loading ", num_files, " Seurat objects"))
	for (i in 1:num_files) {
		seurat.object <- load_seurat(seurat_files[i])
		seurat.objects <- c(seurat.objects,seurat.object)
		seurat.names <- c(seurat.names,get_seurat_basename(seurat_files[i]))
	}
	names(seurat.objects) <- seurat.names
	return(seurat.objects)
}

# Takes a single file location. Must be in .rds format
# Returns the file name
get_seurat_basename <- function(dir = NULL) {
	filename <- basename(dir)
	basename <- substr(filename,1,nchar(filename)-4)
	return(basename)
}

# Takes a directory of outs/filtered_feature_bc_matrix from Cellranger
# Returns a unprocessed seurat object
create_seurat <- function(dir=NULL) {
	require(Seurat)
    tryCatch(
        expr = {
			data <- Read10X(data.dir = dir)
        },
        error = function(e){
            data <- Read10X(data.dir = dir, gene.column=1)
        }
    )
	seurat.object <- CreateSeuratObject(counts = data, project = "seurat", min.cells = 3, min.features = 200)
	return(seurat.object)
}

# Takes a single Seurat object
# Returns the object with standard processing applied
process_seurat <- function(seurat=NULL) {
	require(Seurat)
	seurat.object <- NormalizeData(seurat,verbose = FALSE)
	seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
	seurat.object <- ScaleData(seurat.object)
	seurat.object <- RunPCA(seurat.object)
	seurat.object <- FindNeighbors(seurat.object, dims = 1:20)
	seurat.object <- FindClusters(seurat.object, resolution = 0.8)
	seurat.object <- RunUMAP(seurat.object, dims = 1:20)
	return(seurat.object)
}

# Takes a list of Seurat objects
# Returns the CCA integrated Seurat object
integrate_seurat <- function(object_list=NULL) {
	require(Seurat)
	Integration.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:10, reduction="cca")
	Integration.combined <- IntegrateData(anchorset = Integration.anchors, dims = 1:10)
	DefaultAssay(Integration.combined) <- "integrated"
	Integration.combined <- ScaleData(Integration.combined)
	Integration.combined <- RunPCA(Integration.combined, verbose = TRUE)
	Integration.combined <- RunUMAP(Integration.combined, reduction = "pca", dims = 1:10, return.model=TRUE)
	Integration.combined <- FindNeighbors(Integration.combined, reduction = "pca", dims = 1:10)
	Integration.combined <- FindClusters(Integration.combined, resolution = 0.5)
	return(Integration.combined)
}

# Takes a reference and query Seurat object and a column from the reference metadata
# Returns the query Seurat object with labels transferred 
transfer_labels <- function(reference=NULL,query=NULL,meta=NULL) {
	require(Seurat)
	transfer.anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30)
	predictions <- TransferData(anchorset = transfer.anchors, reference = reference, refdata = meta, dims = 1:30)
	seurat.object <- AddMetaData(query, metadata = predictions)
	return(seurat.object)
}

# Takes a reference and query Seurat object and a column from the reference metadata
# Returns the query Seurat object mapped to the reference 
map_query <- function(reference=NULL,query=NULL,meta=NULL) {
	require(Seurat)
	anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:30, reference.reduction = "pca")
	mapped <- MapQuery(anchorset = anchors, reference = reference, query = query, refdata = meta, reference.reduction = "pca", reduction.model = "umap")
	return(mapped)
}

# Takes the output of map_query, a Seurat object with predicitons and reference UMAP applied
# Optionally, assign a name to the file
# Saves a UMAP plot of the query projected onto the reference
plot_projection <- function(seurat=NULL,name="map_query") {
	require(Seurat)
	require(ggplot2)
	DimPlot(mapped, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,label.size = 6, repel = T) + NoLegend()
	ggsave(file=paste0(name,".png"), width = 10, height = 10)
}

# Takes a Seurat object, a column in the metadata, option for labels or legend, and an optional name
# Saves a UMAP plot of the object
save_umap <- function(seurat=NULL,group=NULL,label=FALSE,name="umap") {
	require(Seurat)
	require(ggplot2)
	if (label) {
		DimPlot(seurat, group.by = group, label = T, repel = TRUE) + NoLegend()
	} else {
		DimPlot(seurat, group.by = group, label = F)
	}
	ggsave(file=paste0(name,".png"), width = 10, height = 10)
}

# Takes a Seurat object and a vector of gene names, and an optional name
# Saves a feature plot of the object
save_featureplot <- function(seurat=NULL,features=NULL,name="features") {
	require(Seurat)
	require(ggplot2)
	FeaturePlot(seurat,features=features,pt.size = 1, min.cutoff ="q10", order=TRUE)
	ggsave(file=paste0(name,".png"), width = 10, height = 10)
}

# Takes a Seurat object, a column in the metadata, and an optional name
# Wries to file the number of cells for each category
save_cell_counts <- function(seurat=NULL,ident=NULL,name="cell_counts") {
	require(Seurat)
	cell_counts <- as.data.frame(table(seurat@meta.data[,ident]))
	colnames(cell_counts) <- c(ident,"Frequency")
	print(cell_counts)
	save_table(cell_counts,name,rows=F)
}

# Takes a Seurat object and an optional name
# Writes a file containing the entirety of its metadata
save_metadata <- function(seurat=NULL, name="metadata") {
	require(Seurat)
	metadata <- seurat@meta.data
	save_table(metadata,name)
}

# Takes a Seurat object and an optional name
# Writes a file containing the raw counts data
save_raw_counts <- function(seurat=NULL,name="counts") {
	require(Seurat)
	counts <- seurat@assays$RNA@counts
	save_table(counts,name)
}

# Takes a Seurat object, a column in the metadata, and optional downsampling to X cells per ident
# Returns the markers for the given idents
get_markers <- function(seurat=NULL,idents=NULL,downsample=Inf) {
	require(Seurat)
	Idents(seurat) <- idents
	markers <- FindAllMarkers(seurat,logfc.threshold = 0.5,max.cells.per.ident=downsample,test.use = "MAST")
	return(markers)
}

# Takes the markers output by get_markers, and the amount of markers that should be kept for each cluster
# Returns the subset markers
subset_markers <- function(markers = NULL,amount=100) {
	markers_subset <- markers %>% arrange(desc(avg_log2FC)) %>% group_by(cluster) %>% slice(1:amount)
	return(markers_subset)
}

# Takes a Seurat object, an ident of interest, and whether the data is human or mouse
# Returns a CellChat object
process_cellchat <- function(seurat=NULL,ident=NULL,species="human") {
	require(CellChat)
	require(Seurat)
	require(future)
	cellchat <- createCellChat(object = seurat, group.by = ident)
	if (species=="human") {
		cellchat@DB <- CellChatDB.human
	} else if (species=="mouse") {
		cellchat@DB <- CellChatDB.mouse
	}
	cellchat <- subsetData(cellchat) 
	future::plan("multiprocess", workers = 4) 
	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)
	cellchat <- computeCommunProb(cellchat)
	cellchat <- filterCommunication(cellchat, min.cells = 10)
	cellchat <- computeCommunProbPathway(cellchat)
	cellchat <- aggregateNet(cellchat)
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
	return(cellchat)
}

# Takes a CellChat object from process_cellchat
# Generates a collection of folders with various figures in the current directory
get_cellchat_figures <- function(cellchat=NULL) {
	require(CellChat)
	require(patchwork)
	dir.create(file.path("Pathways/Circle"),recursive=T)
	dir.create(file.path("Pathways/Heatmap"),recursive=T)
	dir.create(file.path("Pathways/Signaling"),recursive=T)
	dir.create(file.path("Pathways/LR"),recursive=T)
	dir.create("Global")

	groupSize <- as.numeric(table(cellchat@idents))
	png(file="Global/all_interactions.png",width=2500, height=1800, res=200)
	par(mfrow = c(1,2), xpd=TRUE)
	netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
	netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
	dev.off()

	mat <- cellchat@net$weight
	for (i in 1:nrow(mat)) {
		png(file=paste0("Global/",rownames(mat)[i], "_interactions.png"),width=800, height=800)
  		mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  		mat2[i, ] <- mat[i, ]
  		netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
		dev.off()
	}

	png(file="Global/signals_outgoing.png",width=1400, height=2000, res=200)
	plot <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",font.size = 10,width=12,height=18)
	print(plot)
	dev.off()

	png(file="Global/signals_incoming.png",width=1400, height=2000, res=200)
	plot <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",font.size = 10,width=12,height=18)
	print(plot)
	dev.off()

	pathways.show.all <- cellchat@netP$pathways
	for (i in 1:length(pathways.show.all)) {
		png(file=paste0("Pathways/Circle/",pathways.show.all[i], "_circle.png"),width=800, height=800)
		plot <- netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "circle",remove.isolate=F,vertex.label.cex=1.5,pt.title = 15)
		print(plot)
		dev.off()
		png(file=paste0("Pathways/Heatmap/",pathways.show.all[i], "_heatmap.png"),width=800, height=900)
		plot <- netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds", font.size=14, font.size.title = 18)
		print(plot)
		dev.off()
		png(file=paste0("Pathways/Signaling/",pathways.show.all[i], "_signaling_heatmap.png"),width=500, height=300)
		plot <- netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 12, height = 4, font.size = 10)
		print(plot)
		dev.off()
		png(file=paste0("Pathways/LR/",pathways.show.all[i], "_L-R_contribution.png"),width=800, height=800)
		plot <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
		print(plot)
		dev.off()
	}
}

#------------------------------------------------------------------------------------------


# GENERAL PURPOSE
#------------------------------------------------------------------------------------------

# Takes a data.table or data.frame, a name, and optional formatting information
# Writes a table to the current directory
save_table <- function(table=NULL, name=NULL, format="tab", rows = T, cols = T) {
	if (format == "tab") {
		write.table(table, file=paste0(name,".txt"), sep="\t", quote=F, row.names = rows, col.names = cols)
	} else if (format == "csv") {
		write.table(table, file=paste0(name,".csv"), sep=",", quote=F, row.names = rows, col.names = cols)
	}
}

# Takes any R object
# Returns the name of the object as a string
get_obj_name <- function(object=NULL) {
	return(deparse(substitute(object)))
}

# Takes a list of genes in the ENSEMBL format, and whether the data is human or mouse
# Returns a new counts matrix with gene symbols as row names, except in the case of duplicates or no matches being found
convert_genes <- function(genes=NULL,species="human") {
	require("org.Hs.eg.db")
	require("org.Mm.eg.db")
	require("AnnotationDbi")
	if (species == "human") {
		symbol_list <- mapIds(org.Hs.eg.db,
							keys=genes,
							column="SYMBOL",
							keytype="ENSEMBL",
							multiVals="first")
	} else if (species == "mouse") {
		symbol_list <- mapIds(org.Mm.eg.db,
							keys=genes,
							column="SYMBOL",
							keytype="ENSEMBL",
							multiVals="first")
	}
	return(symbol_list)
}

# Takes two matrices with the same rows/columns (first has the correct order and second is being matched to the first), and whether you are comparing by row or column
# Returns the unordered matrix with the corrected order
match_order <- function(ordered=NULL,unordered=NULL,by="row") {
	if (by == "row") {
		return(unordered[match(rownames(ordered),rownames(unordered)),])
	} else if (by == "column") {
		return(unordered[,match(colnames(ordered),colnames(unordered))])
	}
}

#------------------------------------------------------------------------------------------
