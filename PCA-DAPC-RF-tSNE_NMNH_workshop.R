###############################################
###############################################
## script modified for the NMNH workshop, derived from: 
## Derkarabetian S., Castillo S., Peter K.K., Ovchinnikov S., Hedin M. 
## "A Demonstration of Unsupervised Machine Learning in Species Delimitation"
## https://doi.org/10.1016/j.ympev.2019.106562
## https://github.com/shahanderkarabetian/uml_species_delim
###############################################
###############################################

#If you need to install any of the packages below, type: install.packages("package name"), or search for the package in the Packages/Install tab in RStudio.

#required packages and versions used
library("adegenet")
# adegenet 2.1.1
library("randomForest")
# randomForest 4.6-14
library("tsne")
# tsne 0.1-3
library("PCDimension")
# PCDimension 1.1.9
library("mclust")
# mclust 5.4.1
library("cluster")
# cluster 2.0.7-1
library("factoextra")
# factoextra 1.0.5


###############################################
###############################################
# PCA and DAPC
###############################################
###############################################

# import str file. Adjust input file name, n.ind, and n.loc for specific file/dataset.
# example dataset used in the above study and for the workshop
data <- read.structure("Metano_UCE_SNPs_forRscript.str", n.ind=30, n.loc=316, onerowperind=FALSE, col.lab=1, col.pop=3, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)

# scale genetic data for PCA
data_scaled <- scaleGen(data, center=FALSE, scale=FALSE, NA.method=c("zero"), nf)

# PCA, can adjust nf to include more components
pca1 <- dudi.pca(data_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA", pch=16)
s.label(pca1$li, clabel=0.5, grid=0)

# DAPC, determine optimal number of genetic clusters
# set max.n.clust to maximum number of pops
clusters <- find.clusters(data, max.n.clust=20, n.iter=1e6, n.start=10)
grp_k <- nlevels(clusters$grp)
results <- dapc(data, grp=clusters$grp, perc.pca=NULL)
assignplot(results)
compoplot(results)

# PCA with DAPC groupings
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA with DAPC clusters", col=clusters$grp, pch=16)

# in case you need 3D PCA
# requires XQuartz for Mac, not sure for PC
require("rgl")
pca3D <- dudi.pca(data_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=3)
plot3d(pca3D$l1, col=clusters$grp)


###############################################
###############################################
# into the Random Forest, unsupervised
###############################################
###############################################

# convert genind scaled data to factors for randomForest
data_conv <- as.data.frame(data_scaled)
data_conv[is.na(data_conv)] <- ""
data_conv[sapply(data_conv, is.integer)] <- lapply(data_conv[sapply(data_conv, is.integer)], as.factor)
data_conv[sapply(data_conv, is.character)] <- lapply(data_conv[sapply(data_conv, is.character)], as.factor)
nsamp <- nrow(data_conv)

# unsupervised random forest
rftest <- randomForest(data_conv, ntree=5000)

############
# classic multi-dimensional scaling (cMDS)
############

# cMDS with optimal number of components to retain using broken-stick
# may need to adjust number of dimensions if given error
cmdsplot1 <- MDSplot(rftest, clusters$grp, nsamp-1)
cmdsplot_bstick <- bsDimension(cmdsplot1$eig)
cmdsplot2 <- MDSplot(rftest, clusters$grp, cmdsplot_bstick)

# cMDS with optimal DAPC k and DAPC clusters
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters", col=clusters$grp, pch=16)
s.label(cmdsplot2$points, clabel=0.5, grid=0)

############
# alternative clustering algorithms
############

####
# pam clustering
####

# pam clustering on cMDS output with optimal k from DAPC
DAPC_pam_clust_cMDS <- pam(cmdsplot2$points, grp_k)
# cMDS with optimal k of DAPC and clusters via PAM
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters (PAM clustering, cMDS output)", col=DAPC_pam_clust_cMDS$clustering, pch=16)

# pam clustering: test range and select optimal k with highest silhouette value
pop<-c()
for (i in 2:10){
  pop[i]<-mean(silhouette(pam(cmdsplot2$points, i))[, "sil_width"])
}
plot(pop,type = "o", xlab = "K", ylab = "PAM silhouette", main="random forest PAM")

# optional: can also do pam clustering on proximity scores output with optimal k from DAPC
DAPC_pam_clust_prox <- pam(rftest$proximity, grp_k)
# cMDS with optimal k of DAPC and clusters via PAM
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters (PAM clustering, proximity scores)", col=DAPC_pam_clust_prox$clustering, pch=16)

####
# gap statistic
####

# determine optimal k from cMDS using gap statistic with PAM clusters, from proximity scores
# can adjust k.max
cmds_nbclust <- fviz_nbclust(cmdsplot1$points, kmeans, nstart = 25,  method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
cmds_nbclust
cmds_nbclust_k <- cmds_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
cmds_nbclust_clust <- pam(cmdsplot1$points, cmds_nbclust_k)
# cMDS with optimal k of RF via gap statistic and clusters via PAM (euc)
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="proximity RF gap statistic optimal K and clusters (PAM clustering)", col=cmds_nbclust_clust$clustering, pch=16)

####
# hierarchical clustering
####

# determine optimal k from cMDS via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
cmdsplot_clust <- Mclust(cmdsplot2$points, G=1:10)
mclust_grps_cmdsplot <- as.numeric(cmdsplot_clust$classification)
max(mclust_grps_cmdsplot)
# cMDS with optimal k and clusters of RF via hierarchical clustering
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_cmdsplot, pch=16)
mclust_grps_cmdsplot


###############################################
###############################################
# t-SNE
###############################################
###############################################

# prepare plot labels and such
# two ways to set up, depending on input

# 1) this makes it so it is grouped by DAPC clusters
colors = rainbow(length(unique(clusters$grp)))
names(colors) = unique(clusters$grp)
ecb = function(x,y){plot(x,t='n'); text(x, labels=clusters$grp, col=colors[clusters$grp])}

# 2) or, if present, this uses population assignment (a priori species) column in the .str file.
#colors = rainbow(length(unique(data$pop)))
#names(colors) = unique(data$pop)
#ecb = function(x,y){plot(x,t='n'); text(x, labels=data$pop, col=colors[data$pop])}

# t-SNE, on principal components of scaled data
# adjust perplexity, initial_dims
# can do k=3 for 3D plot
# can use data_conv matrix too
tsne_p5 = tsne(pca1$tab, epoch_callback=ecb, max_iter=5000, perplexity=5, initial_dims=5, k=2)

# tSNE plot with DAPC groupings
plot(tsne_p5, main="t-SNE perplexity=5 with DAPC optimal k and clusters", col=clusters$grp, pch=16)
sample_names <- data_scaled[,0]
tsne_plot <- cbind(sample_names, tsne_p5)
s.label(tsne_plot, clabel=0.5, grid=0)

############
# alternative clustering algorithms
############

####
# gap statistic
####

# determine optimal k using gap statistic with pam clustering
# can adjust k.max
tsne_p5_nbclust <- fviz_nbclust(tsne_p5, kmeans, nstart = 25,  method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
tsne_p5_nbclust
tsne_p5_nbclust_k <- tsne_p5_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
tsne_p5_nbclust_clust <- pam(tsne_p5, tsne_p5_nbclust_k)
# t-SNE with optimal k of RF via gap statistic and clusters via PAM
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF gap statistic optimal K and clusters (PAM clustering)", col=tsne_p5_nbclust_clust$clustering, pch=16)

####
# hierarchical clustering
####

# determine optimal k of RF via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
tsne_p5_clust <- Mclust(tsne_p5, G=1:20)
mclust_grps_tsne_p5 <- as.numeric(tsne_p5_clust$classification)
max(mclust_grps_tsne_p5)
# t-SNE p5 with optimal k and clusters of RF via hierarchical clustering
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_tsne_p5, pch=16)
mclust_grps_tsne_p5


###############################################
###############################################
# Output files
###############################################
###############################################

write.csv(pca1$li, "file_PCA.csv")

pdf("PCA_all_output.pdf")
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA", pch=16)
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA with a priori species", col=clusters$grp, pch=16)
s.label(pca1$li, clabel=0.5, grid=0)
dev.off()

write.csv(cmdsplot2$points, "file_cMDS.csv")
write.csv(mclust_grps_cmdsplot, "file_mclust_cmds_groups.csv")

pdf("RF_all_output.pdf")
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS with a priori species", col=clusters$grp, pch=16)
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS a priori species and clusters (PAM clustering)", col=DAPC_pam_clust_prox$clustering, pch=16)
cmds_nbclust
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS RF gap statistic optimal K and clusters (PAM clustering)", col=cmds_nbclust_clust$clustering, pch=16)
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_cmdsplot, pch=16)
s.label(cmdsplot2$points, clabel=0.5, grid=0)
dev.off()

pdf("tSNE_all_output.pdf")
plot(tsne_p5, main="t-SNE perplexity=5 with a priori species", col=clusters$grp, pch=16)
tsne_p5_nbclust
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF gap statistic optimal K and clusters (PAM clustering)", col=tsne_p5_nbclust_clust$clustering, pch=16)
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_tsne_p5, pch=16)
s.label(tsne_plot, clabel=0.5, grid=0)
dev.off()

