# Create on 12/27/2017
# By: Donglei Yin
# Purpose: 

# Define function plotPCA.dy() based on plotPCA() from DESeq2:
# 1. Add sample label to the plot
# 2. Changing PCA Axis from PC1 vs PC2 to any PCi vs PCj(i,j<=4).

# Note:
# Adjust: scale_x_continuous(limits = c(, )) for specific case


plotPCA.dy <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE, i,j,xlimits=c(-10,10)) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  } else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data = d, aes_string(x = paste0("PC",i), y = paste0("PC",j), color = "group",label="name")) + 
    geom_point(size = 2) + xlab(paste0("PC",i, ":",round(percentVar[i] *100), "% variance")) +
    geom_text(hjust=-0.2, vjust=0.2)+
    ylab(paste0("PC",j, ":", round(percentVar[j] * 100), "% variance")) + coord_fixed()+scale_x_continuous(limits = xlimits)
}

