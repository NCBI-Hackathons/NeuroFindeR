library(Seurat)
set.seed(42)
## Read in the raw 10x data into Seurat

rawData <- Read10X(data.dir = "path/to/file.PREmRNA")()

## Initialize the Seurat object with the raw (non-normalized data).
## We chose to include all genes that are expressed in atleast one cell
## and all cells that have at-least 10 genes expressed
## The idea here was to create low thresholds that will be very inclusive
expRaw <- CreateSeuratObject(raw.data = rawDame, min.cells = 1,
                               min.genes=10, 
                               project = "ProjectName")

############################################################################### 
## LogNormalize, Scale and Center the data using existing methods in Seurat
expRaw <- NormalizeData(object = expRaw, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
## Find the most variable genes
expRaw <- FindVariableGenes(object = expRaw, 
                            mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.3, 
                            x.high.cutoff = 4, 
                            y.cutoff = 0.5)
## Scaling the data
expRawNormalized <- ScaleData(object = expRaw, vars.to.regress = "nUMI")
###############################################################################

###############################################################################
## Dimensionality Reduction, Clustering and TSNE in Seurat
## Dimensionality Reduction using PCA
expRawNormalized <- RunPCA(object = expRawNormalized, 
                 pc.genes = expRawNormalized@var.genes, 
                 do.print = FALSE)

## Project PCA scores for all genes in the dataset based on PCA analysis 
expRawNormalized <- ProjectPCA(object = expRawNormalized, do.print = FALSE)

## Generate clusters of cells/cell types based on PCA
expRawNormalized <- FindClusters(object = expRawNormalized, 
                         reduction.type = "pca", 
                         dims.use = 1:20, 
                         resolution = 0.6, 
                         print.output = 0, 
                         save.SNN = TRUE)
## Perform TSNE with Seurat (for visualization)
expRawNormalized <- RunTSNE(object = expRawNormalized, 
                            dims.use = 1:20, 
                            do.fast = TRUE)

## Generate TSNE plot for comparison
TSNEPlot(object = expRawNormalized, do.label = TRUE, pt.size = 0.5)
dev.off()

##############################################################################
# Retrieving raw UMI counts, cluster and PCA info from a Seurat Object
##############################################################################


getClusterInfo <- function(ASeuratObject) {
        ## A function to retrieve raw UMI counts from a Seurat Object
        ## Requires the user to pass a valid Seurat Object
        
        ## Retrieve the raw count matrix
        rawCountMatrix <- as.matrix(ASeuratObject@raw.data)
        
        ## Calculate the total number of UMIs found in each cell
        UMIs <- as.data.frame(colSums(rawCountMatrix))
        colnames(UMIs) <- "UMI_Count"
        
        ## Retrieve the Cluster information from Seurat Info
        clustInfo <- ASeuratObject@ident
        
        ## Retrieve the distance measurement for each cluster from Seurat Object
        distances <- ASeuratObject@dr$pca@cell.embeddings
        
        ## Map the UMI counts to Clusters
        umiClust <- cbind(UMIs,clustInfo)
                     
        ## Add PCA distances to umiClust
        umiClust <- cbind(umiClust,distances)
        ## Return a data frame that has UMI counts, Clustering information
        ## and distance measures for each clusterfor each cell type
        umiClust
}

getClusterDistance <- function(df) {
        
        all.clust <- unique(df$RawClust.ident)
        
        for(cluster in all.clust){
                temp <- df[which(df$RawClust.ident == cluster),3:ncol(df)]
                centroids <- apply(temp, 2, mean)
                dists<-c()
                for(i in 1:ncol(temp)){
                        dists<-cbind(dists,apply(temp[,i, drop = FALSE], 1, function(x) (x - centroids[i]) ^ 2))
                }
                
                x1<-apply(dists, 1, function(x) sqrt(sum(x)))
                x1<-cbind(df[which(df$RawClust.ident == cluster),1:2], x1)
                
        } 
}

umiClust <- getClusterInfo(expRawNormalized)


