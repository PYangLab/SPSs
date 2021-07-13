
## identify SPSs framework
identifySPSs <- function(recurrence, maxFC) {
    library(ggplot2)
    library(PhosR)
    library(MASS)
    library(preprocessCore)
    # fitting the recurrence data with a gamma distribution
    recur <- recurrence[,2]
    names(recur) <- recurrence[,1]
    f1.params <- fitdistr(recur, densfun="gamma")
    x <- recur
    den <- density(x)
    dat <- data.frame(x=den$x, y = den$y)
    
    p1 <- pgamma(recur, shape=f1.params$estimate[1], rate=f1.params$estimate[2], lower.tail = FALSE)
    names(p1) <- names(recur)
    
    ## fitting the fold change data with a gamma distribution
    # filter phosphosites that have very low quantification rates across datasets
    maxFC <- maxFC[which(rowSums(!is.na(maxFC)) > round(ncol(maxFC) / 5)),]
    
    mn <- normalize.quantiles(as.matrix(maxFC))
    rownames(mn) <- rownames(maxFC)
    colnames(mn) <- colnames(maxFC)
    
    x <- apply(mn, 1, mean, na.rm=TRUE)
    
    f2.params <- fitdistr(x, densfun="gamma")
    den <- density(x)
    dat <- data.frame(x=den$x, y = den$y)
    
    p2 <- pgamma(x, shape=f2.params$estimate[1], rate=f2.params$estimate[2], lower.tail = TRUE)
    names(p2) <- rownames(mn)
    
    
    # Fisher's method for combining the two elements
    o <- intersect(names(p1), names(p2))
    
    ps <- cbind(p1[o], p2[o])
    fisher.p <- apply(ps, 1, function(x){
        pchisq(-2*sum(log(x)), 2*length(x), lower.tail = FALSE)
    })
    
    fisher.adj.p <- p.adjust(fisher.p, method = "BH")
    set1 <- fisher.adj.p[which(fisher.adj.p < 0.01)]
    set2 <- fisher.adj.p[which(fisher.adj.p < 0.05)]
    
    return (list(set1, set2, fisher.adj.p))
    
}

## color

scPalette <- function(n) {
    colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                    '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69',
                    '#8DD3C7','#999999')
    if (n <= length(colorSpace)) {
        colors <- colorSpace[1:n]
    } else {
        colors <- grDevices::colorRampPalette(colorSpace)(n)
    }
    return(colors)
}

## compare log2 FC of different SPS sets in different datsets
compareFC <- function(ppe, assay = "normalised") {
    mat <- ppe@assays@data[[assay]]
    rownames(mat) <- paste(toupper(ppe@GeneSymbol), ";", ppe@Residue, ppe@Site, ";", sep = "")
    mat.max <- apply(mat, 1, function(x) max(abs(x)))
    boxplot(mat.max[which(rownames(mat) %in% SPSs)],
            mat.max[which(rownames(mat) %in% mid.rank)],
            mat.max[which(rownames(mat) %in% bot.rank)],
            mat.max[-which(rownames(mat) %in% SPSs)],
            notch = TRUE, col = scPalette(5)[c(1,2,3,5)], names = c("SPSs", "Middle", "Bottom", "non-SPSs"), las = 2)
    paste(wilcox.test(mat.max[which(rownames(mat) %in% SPSs)], 
                      mat.max[which(rownames(mat) %in% mid.rank)], 
                      alternative = "two.sided")$p.value, 
          wilcox.test(mat.max[which(rownames(mat) %in% SPSs)], 
                      mat.max[which(rownames(mat) %in% bot.rank)], 
                      alternative = "two.sided")$p.value,
          wilcox.test(mat.max[which(rownames(mat) %in% SPSs)], 
                      mat.max[-which(rownames(mat) %in% SPSs)], 
                      alternative = "two.sided")$p.value,
          sep = " ")
}

## reform datsets from qPhos dataset to phosphoExperiments Object
mat <- function(idx) {
    load("/Users/dxiao/Dropbox (Sydney Uni)/Mining of publicly available phosphoproteomics data identify common signalome/Datasets/qPhos/datasets.qPhos.Rda")
    l <- length(idx)
    o <- c()
    for (i in idx) {
        rnames <- rownames(as.data.frame(combine.datasets[[i]]))
        o <- union(o, rnames)
    }
    o <- unique(o)
    mat <- matrix(NA, nrow = length(o), ncol = 0)
    for (i in idx) {
        mat <- cbind(mat, as.data.frame(combine.datasets[[i]])[o, 1])
    }
    rownames(mat) <- o
    colnames(mat) <- names(combine.datasets[idx])
    
    ppe <- PhosR::PhosphoExperiment(assays = list(normalised = mat))
    ppe@GeneSymbol <- sapply(strsplit(rownames(ppe), ";"), function(x)x[1])
    ppe@Residue <- sapply(strsplit(rownames(ppe), ";"), function(x)gsub("[0-9]", "", x[2]))
    ppe@Site <- sapply(strsplit(rownames(ppe), ";"), function(x)as.numeric(gsub("[STY]", "", x[2])))
    
    return(ppe)
}

## extract sites from ppe
Sites <- function(ppe) {
    sites <- paste(toupper(ppe@GeneSymbol), ";", paste(ppe@Residue, ppe@Site, sep = ""), ";", sep = "")
    return(sites)
}

getARI <- function(obj.list, method = "kmeans", k) {
    if (method == "kmeans") {
        truth <- obj.list[["truth"]]
        obj.list["truth"] <- NULL
        clusterings <- lapply(lapply(obj.list, t), kmeans, centers = k, iter.max=10)
        results <- sapply(clusterings, FUN = function(i) {
            adjustedRandIndex(truth, i$cluster)
        }, USE.NAMES = T)
    }
    
    return (results)
}

getFMindex <- function(obj.list, method = "kmeans", k) {
    library(dendextend)
    if (method == "kmeans") {
        truth <- obj.list[["truth"]]
        obj.list["truth"] <- NULL
        clusterings <- lapply(lapply(obj.list, t), kmeans, centers = k, iter.max=10)
        results <- sapply(clusterings, FUN = function(i) {
            FM_index(truth, i$cluster)
        }, USE.NAMES = T)
    }
    return (results)
}

getNMI <- function(obj.list, method = "kmeans", k) {
    truth <- obj.list[["truth"]]
    obj.list["truth"] <- NULL
    if (method == "kmeans") {
        clusterings <- lapply(lapply(obj.list, t), kmeans, centers = k, iter.max=10)
        results <- sapply(clusterings, FUN = function(i) {
            igraph::compare(as.numeric(factor(i$cluster)), as.numeric(factor(truth)), method = "nmi")
        }, USE.NAMES = T)
    }
    
    return (results)
}

getPurity <- function(obj.list, method = "kmeans", k) {
    truth <- obj.list[["truth"]]
    obj.list["truth"] <- NULL
    if (method == "kmeans") {
        clusterings <- lapply(lapply(obj.list, t), kmeans, centers = k, iter.max=10)
        results <- sapply(clusterings, FUN = function(i) {
            ClusterPurity(i$cluster, truth)
        }, USE.NAMES = T)
    }
    return (results)
}

ClusterPurity <- function(clusters, labels) {
    sum(apply(table(labels, clusters), 2, max)) / length(clusters)
}

getJaccard <- function(obj.list, method = "kmeans", k) {
    library(clusteval)
    truth <- obj.list[["truth"]]
    obj.list["truth"] <- NULL
    if (method == "kmeans") {
        clusterings <- lapply(lapply(obj.list, t), kmeans, centers = k, iter.max=10)
        results <- sapply(clusterings, FUN = function(i) {
            cluster_similarity(as.numeric(factor(truth)), as.numeric(factor(i$cluster)), similarity = "jaccard", method = "independence")
        }, USE.NAMES = T)
    }
    return (results)
}

# A function to create the color palette 
cite_colorPal <- function(n) {
    
    cols_10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
                 "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")
    
    
    cols_20 <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F",
                 "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6",
                 "#E15759", "#FF9D9A", "#79706E", "#BAB0AC", "#D37295",
                 "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")
    
    if (n <= length(cols_10)) {
        return(cols_10[seq_len(n)])
    } else if (n <= length(cols_20)) {
        return(cols_20[seq_len(n)])
    } else {
        return(scales::gradient_n_pal(cols_20)(seq(0, 1, length.out = n)))
    }
    
}


cite_shapePal <- function(n) {
    return(c(16, 17, 15, 3, 6, 8, 1, 0, 5)[seq_len(n)])
}

sub.corheatmap <- function(mat, low=0.3, midpoint=0.6, high = 1, breaks = 0.1) {
    cormat <- cor(mat)
    melted_cormat <- melt(cormat)
    ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile(color = "white")+
        scale_fill_gradient2(low = "white", high = "#14123D", mid = "#85CFF3", 
                             midpoint = midpoint, limit = c(low,high), space = "Lab",
                             breaks = seq(low,high, breaks),
                             name="Pearson\nCorrelation") +
        theme_minimal()+
        theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                         size = 12, hjust = 1))+
        xlab(grps) +
        coord_fixed()
}




