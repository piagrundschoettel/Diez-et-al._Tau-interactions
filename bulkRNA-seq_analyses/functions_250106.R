

## Heatmap


plotHeatmap <- function(input = norm_anno,
                        geneset = "all",
                        title = "",
                        keyType = "Ensembl",
                        gene_type = "all",
                        cutree_cols = 1,
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        sample_annotation = sample_table,
                        plot_mean = FALSE,
                        method = "complete",
                        distance = "euclidean")
{
  if (geneset[1] != "all") {
    if (keyType == "Ensembl") {
      input <- input[input$GENEID %in% geneset, ]
    }
    else if (keyType == "Symbol") {
      input <- input[input$SYMBOL %in% geneset, ]
    }
    else {
      print("Wrong keyType. Choose Ensembl or Symbol!")
    }
  }
  
  if (gene_type != "all") {
    input <- input[input$GENETYPE %in% gene_type, ]
  }
  
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep = "") 
  
  if(plot_mean == FALSE ){
    input <- input[ , colnames(input) %in% sample_annotation$ID]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[, order(sample_annotation[[plot_order]], decreasing = FALSE)]
    
    column_annotation <- plot_annotation
    
  } else {
    input <- input[ , colnames(input) %in% plot_annotation_mean$condition]
    input_scale <- t(scale(t(input)))
    
    column_annotation <- plot_annotation_mean
  }
  
  breaks <- scaleColors(data = input_scale, maxvalue = 2)[["breaks"]]
  
  pheatmap(input_scale, main = as.character(title),
           show_rownames = show_rownames,
           show_colnames = TRUE,
           cluster_cols = cluster_cols,
           clustering_method = method,
           clustering_distance_rows = distance,
           clustering_distance_cols = distance,
           fontsize = 7,
           cutree_cols = cutree_cols,
           border_color = "black",
           annotation_col = column_annotation,
           annotation_colors = ann_colors,
           breaks = breaks,
           color = scaleColors(data = input_scale, maxvalue = 2)[["color"]])
}


##PCA plot

plotPCA <- function(pca_input = dds_vst,
                    pca_sample_table = sample_table,
                    ntop=500,
                    xPC=1,
                    yPC=2,
                    color,
                    anno_colour,
                    add_density = T,
                    add_silhouette = T,
                    split_by = "NULL",
                    shape="NULL",
                    point_size=3,
                    title="PCA",
                    label = NULL,
                    label_subset = NULL){
  
  if(is.character(pca_input)){
    vst_matrix <- as.matrix(removedbatch_dds_vst)
  }else if(!is.data.frame(pca_input)){
    vst_matrix <- as.matrix(assay(pca_input))
  }else{
    vst_matrix <- pca_input
  }
  
  if(length(ntop)>1){
    pca <- prcomp(t(vst_matrix[dimnames(vst_matrix)[[1]] %in% ntop,]))
  }else if (ntop == "all"){
    pca <- prcomp(t(vst_matrix))
    
  }else{
    # select the ntop genes by variance
    select <- order(rowVars(vst_matrix), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(vst_matrix[select,]))
  }
  
  #calculate explained variance per PC
  explVar <- pca$sdev^2/sum(pca$sdev^2)
  # transform variance to percent
  percentVar <- round(100 * explVar[c(xPC,yPC)], digits=1)
  
  # Define data for plotting
  pcaData <- data.frame(xPC=pca$x[,xPC],
                        yPC=pca$x[,yPC],
                        color = pca_sample_table[[color]],
                        name= as.character(pca_sample_table$ID),
                        stringsAsFactors = F)
  
  #plot PCA
  if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size) 
    }else{
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, color = color, shape=shape)) +
        geom_point(size =point_size) +
        scale_shape_discrete(name=shape)
      
    }
    
    if(add_silhouette){
      pca_plot <- pca_plot + stat_ellipse(geom = "polygon", aes(fill = color), alpha = 0.05)
    }
    
    if(anno_colour[1] == "NULL"){
      pca_plot <- pca_plot + scale_color_discrete(name=color)
    }else{
      pca_plot <- pca_plot + scale_color_manual(values=anno_colour, name=color, aesthetics = c("colour", "fill"))
    }
    
  }else if(is.numeric(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color) +
        scale_fill_gradientn(colours = bluered(100),name=color)
    }else{
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color) +
        scale_fill_gradientn(colours = bluered(100),name=color) +
        scale_shape_discrete(name=shape)
    }
  }
  
  # adds a label to the plot. To label only specific points, put them in the arument label_subset
  if (!is.null(label) == TRUE){
    pcaData$label <- pca_sample_table[[label]]
    if(!is.null(label_subset) == TRUE){
      pcaData_labeled <- pcaData[pcaData$label %in% label_subset,]
    } else {
      pcaData_labeled <- pcaData
    }
    pca_plot <- pca_plot +
      geom_text_repel(data = pcaData_labeled, aes(label = label), nudge_x = 2, nudge_y = 2, colour = "black")
  }
  
  pca_plot <- pca_plot+
    xlab(paste0("PC ",xPC, ": ", percentVar[1], "% variance")) +
    ylab(paste0("PC ",yPC,": ", percentVar[2], "% variance")) +
    coord_fixed()+
    theme_bw()+
    theme(aspect.ratio = 1,
          legend.position = "bottom")+
    ggtitle(title)
  
  if(add_density){
    ggMarginal(pca_plot, groupFill = T, type = "density")
  } else if (split_by != "NULL"){
    pca_plot + facet_wrap(~pca_sample_table[[split_by]])
  } else {
    pca_plot
  }
}


## Heatmap colors


scaleColors <- function(data = input_scale, # data to use
                        maxvalue = NULL # value at which the color is fully red / blue
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  
  ret <- list(breaks = myBreaks, color = myColor)
  return(ret)
}


# Wrapper Function to perform DESeq2 differential testing


DEAnalysis <- function(input = dds,
                       condition,
                       comparison_table = comparison_table,
                       alpha = 0.05,
                       lfcThreshold = 0,
                       sigFC = 2,
                       top = NA,
                       multiple_testing = "IHW",
                       pAdjustMethod = "BH",
                       independentFiltering= TRUE,
                       shrinkage = TRUE,
                       shrinkType = "normal",
                       gene_col = c("GENEID",
                                    "SYMBOL",
                                    "GENETYPE")){
  
  setClass(Class = "DESeq2_analysis_object",
           slots = c(results="data.frame", DE_genes="list", Number_DE_genes="list", rep ="data.frame", rep_sig="data.frame"))
  
  # create results_list
  results_list <- list()
  # print parameters
  results_list$parameters <-list(multiple_testing = multiple_testing,
                                 p_value_threshold = alpha,
                                 log2_FC_threshold = lfcThreshold,
                                 sigFC = sigFC,
                                 shrinkage = shrinkage,
                                 shrinkage_type = shrinkType)
  # Run results() function on comparisons defined in comparison table
  for (x in 1:nrow(comparison_table)){
    # create DE_object
    DE_object <- new(Class = "DESeq2_analysis_object")
    # IHW
    if (multiple_testing=="IHW") {
      res_deseq_lfc <- results(input,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[x]),
                                            paste(comparison_table$control[x])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               filterFun = ihw,
                               pAdjustMethod = pAdjustMethod,
                               altHypothesis = "greaterAbs")
      # Independent Filtering
    }else {
      res_deseq_lfc <- results(input,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[x]),
                                            paste(comparison_table$control[x])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               independentFiltering = independentFiltering,
                               altHypothesis = "greaterAbs",
                               pAdjustMethod= pAdjustMethod)
    }
    if(shrinkage == TRUE){
      if(shrinkType %in% c("normal", "ashr")){
        
        res_deseq_lfc <- lfcShrink(input, 
                                   contrast = c(condition,
                                                paste(comparison_table$comparison[x]),
                                                paste(comparison_table$control[x])),
                                   res=res_deseq_lfc,
                                   type = shrinkType)
        
      }else if(shrinkType == "apeglm"){
        
        res_deseq_lfc <- lfcShrink(input, 
                                   coef = paste0(condition, "_",
                                                 comparison_table$comparison[x], "_vs_",
                                                 comparison_table$control[x]),
                                   res=res_deseq_lfc,
                                   type = shrinkType,
                                   returnList = F)
      }
    }
    res_deseq_lfc <- as.data.frame(res_deseq_lfc)
    # indicate significant DE genes
    if(is.na(top)){
      res_deseq_lfc$regulation <- ifelse(!is.na(res_deseq_lfc$padj)&
                                           res_deseq_lfc$padj <= alpha&
                                           res_deseq_lfc$log2FoldChange > log(sigFC,2),
                                         "up",
                                         ifelse(!is.na(res_deseq_lfc$padj)&
                                                  res_deseq_lfc$padj <= alpha&
                                                  res_deseq_lfc$log2FoldChange < -log(sigFC,2),
                                                "down",
                                                "n.s."))
    }else{
      res_deseq_lfc$regulation <- ifelse(!is.na(res_deseq_lfc$padj)&
                                           res_deseq_lfc$padj <= alpha&
                                           res_deseq_lfc$log2FoldChange > log(sigFC,2),
                                         "up",
                                         ifelse(!is.na(res_deseq_lfc$padj)&
                                                  res_deseq_lfc$padj <= alpha&
                                                  res_deseq_lfc$log2FoldChange < -log(sigFC,2),
                                                "down",
                                                "n.s."))
    }
    # add gene annotation to results table
    res_deseq_lfc$GENEID <- row.names(res_deseq_lfc) 
    res_deseq_lfc <- merge(res_deseq_lfc,
                           norm_anno[,gene_col],
                           by = "GENEID")
    row.names(res_deseq_lfc) <- res_deseq_lfc$GENEID
    res_deseq_lfc$comparison<-paste(comparison_table$comparison[x]," vs ",comparison_table$control[x],
                                    sep="")
    # re-order results table
    if (multiple_testing=="IHW") {
      res_deseq_lfc<-res_deseq_lfc[,c(gene_col,
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "pvalue",
                                      "padj"
      )]
    }else{
      res_deseq_lfc<-res_deseq_lfc[,c(gene_col,
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "pvalue",
                                      "padj")]
    }
    # print result table
    DE_object@results <- res_deseq_lfc
    # print DE genes in separate tables
    DE_object@DE_genes <- list(up_regulated_Genes = res_deseq_lfc[res_deseq_lfc$regulation =="up",],
                               down_regulated_Genes= res_deseq_lfc[res_deseq_lfc$regulation =="down",])
    # print the numbers of DE genes
    DE_object@Number_DE_genes <- list(up_regulated_Genes = nrow(DE_object@DE_genes$up_regulated_Genes),
                                      down_regulated_Genes= nrow(DE_object@DE_genes$down_regulated_Genes))
    # write DE_object into results_list
    results_list[[paste(comparison_table$comparison[x], "vs", comparison_table$control[x], sep="_")]] <- DE_object
  }
  return(results_list)
}



### Union of DEgenes

uDEG <- function(input = DEresults, keyType = "Ensembl", comparisons){
  uDEGs <- NULL
  tmp <- input[names(input) %in% comparisons]
  for(i in 1:length(comparisons)){
    DEGs <- as.data.frame(tmp[[i]]@results[tmp[[i]]@results$regulation %in% c("up","down"),])
    if(!keyType %in% c("Ensembl", "Symbol")){stop("keyType should be one of Ensembl or Symbol")}
    if(keyType == "Ensembl"){
      uDEGs <- unique(c(uDEGs, DEGs$GENEID))
    } else if (keyType == "Symbol"){
      uDEGs <- unique(c(uDEGs, DEGs$SYMBOL))
    }
  }
  uDEGs
}

### Heatmaps of DE genes based on pheatmap


plotDEHeatmap <- function(input=norm_anno,
                          sample_annotation=sample_table,
                          column_annotation=plot_annotation,
                          comparison,
                          factor,
                          gene_anno=tx_annotation,
                          conditions="all",
                          gene_type="all",
                          show_rownames = FALSE,
                          cluster_cols = FALSE,
                          plot_mean = FALSE){
  
  geneset <- DEresults[[comparison]]@results[DEresults[[comparison]]@results$regulation %in% c("up","down"),"GENEID"]
  
  input <- input[input$GENEID %in% geneset,]
  
  if(plot_mean){
    
    if(conditions[1] == "all"){
      input <- input[ ,colnames(input) %in% as.vector(sample_annotation[[factor]])]
      input_scale <- t(scale(t(input)))
    } else {
      input <- input[ ,colnames(input) %in% sample_annotation[as.vector(sample_annotation[[factor]]) %in% conditions,][[factor]]]
      input_scale <- t(scale(t(input)))
      sample_annotation<-subset(sample_annotation,sample_annotation[[factor]] %in% conditions)
    }
    
    
  } else {
    if(conditions[1] == "all"){
      input <- input[,colnames(input) %in% sample_annotation$ID]
      input_scale <- t(scale(t(input)))
    } else {
      input <- input[,colnames(input) %in% sample_annotation[as.vector(sample_annotation[[factor]]) %in% conditions,]$ID]
      input_scale <- t(scale(t(input)))
      sample_annotation<-subset(sample_annotation,sample_annotation[[factor]] %in% conditions)
    }
  }
  
  input_scale<-as.data.frame(input_scale)
  input_scale$GENEID <- rownames(input_scale)
  gene_anno <- gene_anno[match(rownames(input_scale), gene_anno$GENEID),]
  input_scale <- merge(input_scale,
                       gene_anno,
                       by = "GENEID")
  rownames(input_scale) <- input_scale$GENEID
  
  title=paste("Heatmap of significant DE genes in: ",comparison,sep="")
  
  plotHeatmap(input=input_scale,
              sample_annotation=sample_annotation,
              geneset = geneset,
              title = title,
              keyType = "Ensembl",
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              gene_type=gene_type,
              plot_mean = plot_mean)
}

##GSEA enrichment

enrichGSEA <- function(comparisons,
                       DE_results = DEresults,
                       universe = universe,
                       GeneSets = c("GO","KEGG", "HALLMARK", "REACTOME"), 
                       pCorrection = "bonferroni", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05,
                       col_gene= "SYMBOL"){
  
  
  ## Results lists
  GOresults.list <- list()
  KEGGresults.list <- list()
  HALLMARKresults.list <- list()
  REACTOMEresults.list <- list()
  positionalresults.list <- list()
  
  ## Enrichment
  for(i in comparisons){
    
    print(paste0("Performing enrichment for ", i))
    
    ## DEGs
    genes <- DE_results[[i]]@results
    DE_up <- unlist(c(genes[genes$regulation=="up",][col_gene]))
    DE_down <- unlist(c(genes[genes$regulation=="down",][col_gene]))
    
    #Compare to positional geneset
    if("positional" %in% GeneSets){
      print("Performing positional enrichment")
      respositionalup <- as.data.frame(enricher(gene = DE_up,
                                                universe = universe,
                                                TERM2GENE= positional_terms,
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff  = qvalueCutoff))
      if(nrow(respositionalup)>0){
        respositionalup$comparison <- paste0(as.character(i), "_up")
        respositionalup$regulation <- "up"
        positionalresults.list[[paste(i, "_up")]] <- respositionalup
      }
      respositionaldown <- as.data.frame(enricher(gene = DE_down,
                                                  universe = universe,
                                                  TERM2GENE= positional_terms,
                                                  pAdjustMethod = pCorrection,
                                                  pvalueCutoff  = pvalueCutoff,
                                                  qvalueCutoff  = qvalueCutoff))
      if(nrow(respositionaldown)>0){
        respositionaldown$comparison <- paste0(as.character(i), "_down")
        respositionaldown$regulation <- "down"
        positionalresults.list[[paste(i, "_down")]] <- respositionaldown
      }
    }
    
    # Compare the Clusters regarding their GO enrichment
    if("GO" %in% GeneSets){
      print("Performing GO enrichment")
      resGOup <- as.data.frame(enricher(gene = DE_up,
                                        universe = universe,
                                        TERM2GENE= GO_terms,
                                        pAdjustMethod = pCorrection,
                                        pvalueCutoff  = pvalueCutoff,
                                        qvalueCutoff  = qvalueCutoff))
      if(nrow(resGOup)>0){
        resGOup$comparison <- paste0(as.character(i), "_up")
        resGOup$regulation <- "up"
        GOresults.list[[paste(i, "_up")]] <- resGOup
      }
      resGOdown <- as.data.frame(enricher(gene = DE_down,
                                          universe = universe,
                                          TERM2GENE= GO_terms,
                                          pAdjustMethod = pCorrection,
                                          pvalueCutoff  = pvalueCutoff,
                                          qvalueCutoff  = qvalueCutoff))
      if(nrow(resGOdown)>0){
        resGOdown$comparison <- paste0(as.character(i), "_down")
        resGOdown$regulation <- "down"
        GOresults.list[[paste(i, "_down")]] <- resGOdown
      }
    }
    # Compare the Clusters regarding their HALLMARK enrichment
    if("HALLMARK" %in% GeneSets){
      print("Performing HALLMARK enrichment")
      resHMup <- as.data.frame(enricher(gene = DE_up,
                                        universe = universe,
                                        TERM2GENE = HALLMARK_terms,
                                        pAdjustMethod = pCorrection,
                                        pvalueCutoff  = pvalueCutoff,
                                        qvalueCutoff  = qvalueCutoff))
      if(nrow(resHMup)>0){
        resHMup$comparison <- paste0(as.character(i), "_up")
        resHMup$regulation <- "up"
        HALLMARKresults.list[[paste(i, "_up")]] <- resHMup
      }
      resHMdown <- as.data.frame(enricher(gene = DE_down,
                                          universe = universe,
                                          TERM2GENE = HALLMARK_terms,
                                          pAdjustMethod = pCorrection,
                                          pvalueCutoff  = pvalueCutoff,
                                          qvalueCutoff  = qvalueCutoff))
      if(nrow(resHMdown)>0){
        resHMdown$comparison <- paste0(as.character(i), "_down")
        resHMdown$regulation <- "down"
        HALLMARKresults.list[[paste(i, "_down")]] <- resHMdown
      }
    }
    # Compare the Clusters regarding their KEGG enrichment
    if("KEGG" %in% GeneSets){
      print("Performing KEGG enrichment")
      resKEGGup <- as.data.frame(enricher(gene = DE_up,
                                          universe = universe,
                                          TERM2GENE = KEGG_terms,
                                          pAdjustMethod = pCorrection,
                                          pvalueCutoff  = pvalueCutoff,
                                          qvalueCutoff  = qvalueCutoff))
      if(nrow(resKEGGup)>0){
        resKEGGup$comparison <- paste0(as.character(i), "_up")
        resKEGGup$regulation <- "up"
        KEGGresults.list[[paste(i, "_up")]] <- resKEGGup
      }
      resKEGGdown <- as.data.frame(enricher(gene = DE_down,
                                            universe = universe,
                                            TERM2GENE = KEGG_terms,
                                            pAdjustMethod = pCorrection,
                                            pvalueCutoff  = pvalueCutoff,
                                            qvalueCutoff  = qvalueCutoff))
      if(nrow(resKEGGdown)>0){
        resKEGGdown$comparison <- paste0(as.character(i), "_down")
        resKEGGdown$regulation <- "down"
        KEGGresults.list[[paste(i, "_down")]] <- resKEGGdown
      }
    }
    # Compare the Clusters regarding their REACTOME enrichment
    if("REACTOME" %in% GeneSets){
      print("Performing REACTOME enrichment")
      resREACTOMEup <- as.data.frame(enricher(gene = DE_up,
                                              universe = universe,
                                              TERM2GENE = REACTOME_terms,
                                              pAdjustMethod = pCorrection,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff  = qvalueCutoff))
      if(nrow(resREACTOMEup)>0){
        resREACTOMEup$comparison <- paste0(as.character(i), "_up")
        resREACTOMEup$regulation <- "up"
        REACTOMEresults.list[[paste(i, "_up")]] <- resREACTOMEup
      }
      resREACTOMEdown <- as.data.frame(enricher(gene = DE_down,
                                                universe = universe,
                                                TERM2GENE = REACTOME_terms,
                                                pAdjustMethod = pCorrection,
                                                pvalueCutoff  = pvalueCutoff,
                                                qvalueCutoff  = qvalueCutoff))
      if(nrow(resREACTOMEdown)>0){
        resREACTOMEdown$comparison <- paste0(as.character(i), "_down")
        resREACTOMEdown$regulation <- "down"
        REACTOMEresults.list[[paste(i, "_down")]] <- resREACTOMEdown
      }
    }
    
  }
  positionalresults <- do.call("rbind", positionalresults.list)
  GOresults <- do.call("rbind", GOresults.list)
  HALLMARKresults <- do.call("rbind", HALLMARKresults.list)
  KEGGresults <- do.call("rbind", KEGGresults.list)
  REACTOMEresults <- do.call("rbind", REACTOMEresults.list)
  
  results.list <- list("GO" = GOresults,
                       "HALLMARK" = HALLMARKresults,
                       "KEGG" = KEGGresults,
                       "REACTOME" = REACTOMEresults,
                       "positional" = positionalresults)
  return(results.list)
}




dotplotGSEA <- function(enrich.df,
                        show =10,
                        orderBy = c("padj", "count", "GeneRatio"),
                        colorBy = c("padj", "regulation"),
                        scaleBy = c("count", "GeneRatio")
){
  if(nrow(enrich.df)<1){ print("No enrichment found.")
  } else{
    # ordering
    
    enrich.df$gene.ratio <- sapply(strsplit(enrich.df$GeneRatio, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    
    if(orderBy=="padj"){
      enrich.df %>% group_by(comparison) %>% arrange(desc(Count), .by_group = TRUE) -> x
      x %>% group_by(comparison) %>% arrange(desc(Count), .by_group = TRUE) %>% dplyr::top_n(n = -show, wt = p.adjust) %>% dplyr::pull(Description) -> terms
      x <- x[x$Description %in% terms,]
      x$Description <- ifelse(nchar(x$Description)>80, paste(substr(x$Description, 1, 80),"[...]",sep=""), x$Description)
      x$Description <- factor(x$Description, levels = rev(unique(x$Description)))
      
    } else if (orderBy=="count"){
      enrich.df %>% group_by(comparison) %>% arrange(p.adjust, .by_group = TRUE) -> x
      x %>% group_by(comparison) %>% arrange(p.adjust, .by_group = TRUE) %>% dplyr::top_n(n = show, wt = Count) %>% dplyr::pull(Description) -> terms
      x <- x[x$Description %in% terms,]
      x$Description <- ifelse(nchar(x$Description) > 80, paste(substr(x$Description, 1, 80),"[...]",sep=""), x$Description)
      x$Description <- factor(x$Description, levels = rev(unique(x$Description)))
      
    } else if (orderBy=="GeneRatio"){
      enrich.df %>% group_by(comparison) %>% arrange(p.adjust, .by_group = TRUE) -> x
      x %>% group_by(comparison) %>% arrange(p.adjust, .by_group = TRUE) %>% dplyr::top_n(n = show, wt = gene.ratio) %>% dplyr::pull(Description) -> terms
      x <- x[x$Description %in% terms,]
      x$Description <- ifelse(nchar(x$Description) > 80, paste(substr(x$Description, 1, 80),"[...]",sep=""), x$Description)
      x$Description <- factor(x$Description, levels = rev(unique(x$Description)))
    }else {
      "Please select orderBy='padj', orderBy='count'or orderBy='GeneRatio'."
    }
    
    ## coloring
    if(colorBy=="regulation" & scaleBy == "count"){
      p <- ggplot(data = x ,aes(x=comparison, y=Description, size=Count , fill=regulation)) +
        geom_point(pch=21) +
        scale_fill_manual(values = c(up = "firebrick4", down = "dodgerblue3")) +
        xlab("") +
        ylab("") +
        scale_radius() +
        theme_linedraw()+
        theme(axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(angle=90, vjust=0.5,hjust=1, color = "black"),
              text = element_text(size = 12),
              panel.grid.major = element_line(colour = "grey70"))
      plot(p)
      
    } else if (colorBy=="regulation" & scaleBy == "GeneRatio"){
      p <- ggplot(data = x ,aes(x=comparison, y=Description, size=gene.ratio , fill=regulation)) +
        geom_point(pch=21) +
        scale_fill_manual(values = c(up = "firebrick4", down = "dodgerblue3")) +
        xlab("") +
        ylab("") +
        scale_radius() +
        theme_linedraw()+
        theme(axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(angle=90, vjust=0.5,hjust=1, color = "black"),
              text = element_text(size = 12),
              panel.grid.major = element_line(colour = "grey70"))
      plot(p)
      
    } else if (colorBy=="padj" & scaleBy=="count"){
      p <- ggplot(data = x, aes(x=comparison, y=Description, size=Count, color=p.adjust)) +
        geom_point(pch=21) +
        scale_colour_gradientn(colours = c('red', 'orange', 'darkblue', 'darkblue'),
                               limits = c(0, 1),
                               values = c(0, 0.05, 0.2, 0.5, 1),
                               breaks = c(0.05, 0.2, 1),
                               labels = format(c(0.05, 0.2, 1))) +
        xlab("") +
        ylab("") +
        scale_radius() +
        theme_linedraw()+
        theme(axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, color = "black"),
              text = element_text(size = 12),
              panel.grid.major = element_line(colour = "grey70"))
      plot(p)
      
    } else if (colorBy=="padj" & scaleBy=="GeneRatio"){
      p <- ggplot(data = x, aes(x=comparison, y=Description, size=gene.ratio, color=p.adjust)) +
        geom_point(pch=21) +
        scale_colour_gradientn(colours = c('red', 'orange', 'darkblue', 'darkblue'),
                               limits = c(0, 1),
                               values = c(0, 0.05, 0.2, 0.5, 1),
                               breaks = c(0.05, 0.2, 1),
                               labels = format(c(0.05, 0.2, 1))) +
        xlab("") +
        ylab("") +
        scale_radius() +
        theme_linedraw()+
        theme(axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(angle=90, vjust=0.5,hjust=1, color = "black"),
              text = element_text(size = 12),
              panel.grid.major = element_line(colour = "grey70"))
      plot(p)
      
    }else {
      "Please select colorBy='regulation' or colorBy='padj' and scaleBy='count' or scaleBy='GeneRatio'."
    }
  }
}


##PC regression

plot_pc_regression <- function(input = removedbatch_dds_vst,
                               meta = sample_table,
                               ntop = "all",
                               meta_variables, 
                               nPCs = 10,
                               title = "PC contribution"){
  
  #compute the PCA embedding outside the original function (optionally select the most variable features of the data)
  
  if(class(input) == "DESeqTransform"){
    if (length(ntop)>1){
      select <- dimnames(as.matrix(assay(input)))[[1]] %in% ntop
    } else if(ntop == "all") {
      select <- order(rowVars(as.matrix(assay(input))), decreasing=TRUE)[1:nrow(as.matrix(assay(input)))]
    } else {
      select <- order(rowVars(as.matrix(assay(input))), decreasing=TRUE)[1:ntop]
    }
    pca <- prcomp(t(as.matrix(assay(input))[select,]))
    
  } else if(class(input) == "data.frame") {
    if(length(ntop)>1){
      select <- rownames(input) %in% ntop
    } else if(ntop == "all") {
      select <- order(rowVars(as.matrix(input)), decreasing=TRUE)[1:nrow(input)]
    } else {
      select <- order(rowVars(input), decreasing=TRUE)[1:ntop]
    }
    pca <- prcomp(t(as.matrix(input)[select, ]))
    
  } else { print("unknown input format")}
  df_pca <- as.data.frame(pca[["x"]])
  
  
  # Calculate variance attributed to metadata
  M <- meta[ , which(colnames(meta) %in% meta_variables)]
  
  for(i in colnames(M)){
    if(length(unique(M[ ,i])) <2 ){
      print(paste("exclude", i, sep = " "))
    }
  }
  
  pc_adj_r_squared <- matrix(NA, ncol = dim(df_pca)[2], nrow = dim(M)[2])
  for(i in 1:dim(df_pca)[2]){
    for(j in 1:dim(M)[2]){
      pc_adj_r_squared[j,i] <- summary(lm(df_pca[,i] ~ M[,j], na.action = na.exclude))$adj.r.squared
    }
  }
  
  pc_adj_r_squared <- as.data.frame(pc_adj_r_squared)
  colnames(pc_adj_r_squared) <- colnames(df_pca)
  rownames(pc_adj_r_squared) <- colnames(M)
  
  df <- pc_adj_r_squared[, 1:nPCs]
  
  paletteLength<-50
  my_palette <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  breakList <- c(seq(min(df), 0, length.out = ceiling(paletteLength/2) + 1), 
                 seq(max(df)/paletteLength, max(df), length.out = floor(paletteLength/2)))
  
  hm <- pheatmap(as.matrix(df),
                 cluster_rows = F,
                 cluster_cols = F,
                 show_rownames = T,
                 show_colnames = T,
                 scale = "none",
                 main = title,
                 display_numbers = T,
                 color= my_palette,
                 breaks = breakList)
}

#Boxplot of single Genes
plotSingleGene <-function(data,
                          symbol  ,
                          condition="cond_set" ,
                          anno_colour=col_cond_set,
                          shape = NULL,
                          ncol = 3,
                          ylab = "Normalized counts"){
  between
  
  input<-as.data.frame(data)
  rownames(input)<- input$GENEID
  
  if(length(symbol[symbol %in% input$GENEID]) == 0){
    stop("None of the genes are present in your input.")
  }else{
    symbol.all <- symbol
    symbol <- symbol[symbol %in% input$GENEID]
    
    if(length(symbol.all) > length(symbol)) {
      print(paste0(paste(symbol.all[!symbol.all %in% symbol], collapse = ", "), " is/are not present"))
    }
    
    plots<-list()
    for (i in 1:length(symbol)) {
      geneCounts <- as.data.frame(t(input[input$GENEID == symbol[i], colnames(input) %in% sample_table$ID]))
      geneCounts$condition <- factor(sample_table[[condition]],
                                     levels = levels(sample_table[[condition]]))
      
      
      for (k in 1){
        GENEID<-colnames(geneCounts)[k]
        colnames(geneCounts)[k]<-"y"
        
        if(!is.null(anno_colour)){
          if (is.null(shape)){
            plot<-ggplot(geneCounts, aes(x = condition, y = y, fill=condition)) +
              geom_boxplot(width=.75,alpha=1, outlier.shape = NA) + 
              stat_boxplot(geom ='errorbar',width=.25) +
              geom_jitter(shape = 21, width = 0.25, na.rm=T, alpha = 0.75) +
              scale_fill_manual(values=anno_colour) 
          }else{
            geneCounts$shape <- sample_table[[shape]]
            legend_shape<-paste0(shape)
            plot<-ggplot(geneCounts, aes_string(x = condition, y = "y", fill=condition)) +
              scale_fill_manual(values=anno_colour)+
              geom_boxplot(width=.75,alpha=1, outlier.shape = NA) + 
              stat_boxplot(geom ='errorbar',width=.25) +
              geom_jitter(shape = 21, width = 0.25, na.rm=T, alpha = 0.75) +
              scale_shape(name=legend_shape)
          }
        }else{
          if (is.null(shape_opt)){
            plot<-ggplot(geneCounts, aes_string(x = condition, y = "y", fill=condition)) +
              geom_boxplot(width=.75,alpha=1, outlier.shape = NA) + 
              stat_boxplot(geom ='errorbar',width=.25) +
              geom_jitter(shape = 21, width = 0.25, na.rm=T, alpha = 0.75) +
              scale_color_brewer(palette = "Spectral")
          }else{
            geneCounts$shape <- sample_table[[shape]]
            legend_shape<-paste0(shape)
            plot<-ggplot(geneCounts, aes_string(x = condition, y = "y", fill=condition)) +
              geom_boxplot(width=.75,alpha=1, outlier.shape = NA) + 
              stat_boxplot(geom ='errorbar',width=.25) +
              geom_jitter(shape = 21, width = 0.25, na.rm=T, alpha = 0.75) +
              scale_color_brewer(palette = "Spectral")+
              scale_shape(name=legend_shape)
          }
        }
        
        
        plots[[length(plots)+1]]<-plot+
          ylab(ylab) +
          scale_y_continuous(expand=c(0.05,0.25)) +
          expand_limits(y=0) +
          labs(title=paste(symbol[i], GENEID, sep=": "))+
          theme_classic()+
          theme(plot.title = element_text(hjust=0.5),
                axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
        
        colnames(geneCounts)[k] <- GENEID
      }
    }
    if(length(plots)>1){
      patchwork::wrap_plots(plots, guides = "collect", ncol = ncol)
    }else{
      plots
    }
  }}

### Volcano Plot
plotVolcano <-  function(DEresults_obj = DEresults,
                         comparison,
                         labelnum=20,
                         y.max = NULL,
                         signature = NULL){
  # specify labeling
  upDE <-  as.data.frame(DEresults_obj[[comparison]]@results[DEresults_obj[[comparison]]@results$regulation =="up",])
  if(!is.null(signature))
    upDE <- upDE %>% dplyr::filter(.,SYMBOL %in% signature)
  FClabel_up <- upDE[order(abs(upDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_up)>labelnum){
    FClabel_up <- as.character(FClabel_up[c(1:labelnum),"GENEID"])
  } else {
    FClabel_up <- as.character(FClabel_up$GENEID)}
  plabel_up <- upDE[order(upDE$padj, decreasing = FALSE),]
  if(nrow(plabel_up)>labelnum){
    plabel_up <- as.character(plabel_up[c(1:labelnum),"GENEID"])
  } else {
    plabel_up <- as.character(plabel_up$GENEID)}
  
  downDE <-  as.data.frame(DEresults_obj[[comparison]]@results[DEresults_obj[[comparison]]@results$regulation =="down",])
  if(!is.null(signature))
    downDE <- downDE %>% dplyr::filter(SYMBOL %in% signature)
  FClabel_down <- downDE[order(abs(downDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_down)>labelnum){
    FClabel_down <- as.character(FClabel_down[c(1:labelnum),"GENEID"])
  } else {
    FClabel_down <- as.character(FClabel_down$GENEID)}
  plabel_down <- downDE[order(downDE$padj, decreasing = FALSE),]
  if(nrow(plabel_down)>labelnum){
    plabel_down <- as.character(plabel_down[c(1:labelnum),"GENEID"])
  } else {
    plabel_down <- as.character(plabel_down$GENEID)}
  
  
  label<- unique(c(FClabel_up, plabel_up, FClabel_down, plabel_down))
  
  data <- DEresults_obj[[comparison]]@results
  data$label<- ifelse(data$GENEID %in% label == "TRUE",as.character(data$GENEID), "")
  data <- data[,colnames(data) %in% c("label", "log2FoldChange", "padj", "regulation")]
  data$color <- apply(data, 1, function(x){
    if(x["label"] != "") "black"
    else x["regulation"]
  }) %>% unlist()
  limits <- ceiling(max(DEresults_obj[[comparison]]@results$log2FoldChange))
  
  # Volcano Plot
  volcano <- ggplot(data=na.omit(data), aes(x=log2FoldChange, y=-log10(padj), fill=regulation, color = color)) +
    geom_point(shape = 21, alpha=0.75, size=1.75) +
    scale_color_manual(values = c(c("down" = "dodgerblue3","n.s." = "grey","up" = "firebrick", "black" = "black"))) +
    scale_fill_manual(values=c("down" = "dodgerblue3", "n.s." = "grey", "up" = "firebrick"))+
    xlab("log2(FoldChange)") +
    ylab("-log10(padj)") +
    geom_vline(xintercept = c(-log(DEresults_obj$parameters$sigFC,2),log(DEresults_obj$parameters$sigFC,2)), colour="darkgrey", linetype = "dashed")+
    geom_hline(yintercept=-log(0.05,10),colour="darkgrey", linetype = "dashed")+
    geom_text_repel(data=na.omit(data[!data$label =="",]),aes(label=label), size=3)+
 #   scale_x_continuous(limits = c(-(limits), limits), breaks = seq(-(limits), limits, by=1))+
    guides(color="none") +
    ggtitle(paste("Volcano Plot of: ",comparison,sep="")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if(!is.null(y.max))
    volcano <- volcano +
    scale_y_continuous(limits = c(0, y.max), breaks = seq(0, y.max, by=1))
  else 
    volcano <- volcano +
    scale_y_continuous()
  
  volcano
}

#DEGs for each chromosome
DEGs_chrom <- function(DEres = DEresults, i){
  results_df <- DEres[[i]]@results[DEres[[i]]@results$regulation != "n.s.",]
  DEGs_chrom <- rbind(as.data.frame(table(results_df$chromosome_name, results_df$regulation)),as.data.frame(table(results_df$chromosome_name,results_df$comparison)))
  DEGs_chrom <- cbind(as.matrix(table(results_df$chromosome_name, results_df$regulation)),as.matrix(table(results_df$chromosome_name)))
  DEGs_chrom <- DEGs_chrom[match(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19", "X"), rownames(DEGs_chrom)),]
  colnames(DEGs_chrom) <- c("down", "up", "all")
  print(Heatmap(DEGs_chrom, cluster_rows = F, cluster_columns = F, row_names_side = "left",column_title_side ="top", col=viridis(100), column_title =i, cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", DEGs_chrom[i, j]), x, y, gp = gpar(fontsize = 10))}))
}
