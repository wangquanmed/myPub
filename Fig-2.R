# General Preset
library(magrittr)
library(reshape2)

expr = round(na.omit(expr))
N_value0 = apply(expr,1,function(x) table(x)["0"])
expr = expr[c(N_value0 < (ncol(expr) - 1)/20 | is.na(N_value0) == T),]#g exp>1/20 sp
colnames(expr) = substr(colnames(expr),1,16)
expr %<>% 
  .[,(
    (substr(colnames(.),14,16) %in% c("01A","11A")) &
      (substr(colnames(.),1,12) %in% p_info$PID)) | 
      colnames(.) == "entrez"]

cons = colnames(expr)
exps = cons[substr(cons,14,16) == "01A"]
cons = cons[substr(cons,14,16) == "11A"]
expr_tumor = expr[,c("entrez",exps)]

getCount2Norm(r = expr_tumor,
              r_gname = "entrez",
              symbolORgeneid = "gene_id",
              mmORhs = "hs",
              to = "tpm") ->
  expr_tumor_tpm
data.frame(symbol = expr_tumor_tpm$symbol,
           N = round(
             rowMeans(
               expr_tumor_tpm[,
                              intersect(colnames(expr_tumor_tpm),
                                        paste0(c(n0,n1a,n1b),".01A"))]),3)) ->
  expr_tumor_tpm_mean
merge(expr_tumor_tpm_mean,
      data.frame(symbol = expr_tumor_tpm$symbol,
                 n0 = round(rowMeans(expr_tumor_tpm[,intersect(colnames(expr_tumor_tpm),
                                                               paste0(n0,".01A"))]),3),
                 n1a = round(rowMeans(expr_tumor_tpm[,intersect(colnames(expr_tumor_tpm),
                                                                paste0(n1a,".01A"))]),3),
                 n1b = round(rowMeans(expr_tumor_tpm[,intersect(colnames(expr_tumor_tpm),
                                                                paste0(n1b,".01A"))]),3)),
      by = "symbol",
      all.x = T) ->
  expr_tumor_tpm_mean
expr_tumor_tpm_mean = getDF_sorted(DF = expr_tumor_tpm_mean,by = "N")
#write.csv(expr_tumor_tpm_mean,file = "./data/expr_tumor_tpm_mean.csv",row.names = F)

load("./data/x");getsave <- function(){
  getX(expr = expr_tumor,
       gname = "entrez",
       cons = intersect(paste0(n0,".01A"),colnames(expr_tumor)),
       exps = intersect(paste0(n1a,".01A"),colnames(expr_tumor)),
       method = "rc-deseq2") %>% 
    merge(g_info,.,by.x = "gene_id",by.y = "entrez") %>% 
    getDF_sorted(.,by = "pValue") ->
    x_n0_n1a
  colnames(x_n0_n1a)
  x_n0_n1a = x_n0_n1a[,c("gene_id","pValue","log2FoldChange")]
  colnames(x_n0_n1a) = c("gene_id",paste0(colnames(x_n0_n1a),"_n0_n1a")[2:3])
  getX(expr = expr_tumor,
       gname = "entrez",
       cons = intersect(paste0(n0,".01A"),colnames(expr_tumor)),
       exps = intersect(paste0(n1b,".01A"),colnames(expr_tumor)),
       method = "rc-deseq2") %>% 
    merge(g_info,.,by.x = "gene_id",by.y = "entrez") %>% 
    getDF_sorted(.,by = "pValue") ->
    x_n0_n1b
  colnames(x_n0_n1b)
  x_n0_n1b = x_n0_n1b[,c("gene_id","pValue","log2FoldChange")]
  colnames(x_n0_n1b) = c("gene_id",paste0(colnames(x_n0_n1b),"_n0_n1b")[2:3])
  getX(expr = expr_tumor,
       gname = "entrez",
       cons = intersect(paste0(n1a,".01A"),colnames(expr_tumor)),
       exps = intersect(paste0(n1b,".01A"),colnames(expr_tumor)),
       method = "rc-deseq2") %>% 
    merge(g_info,.,by.x = "gene_id",by.y = "entrez") %>% 
    getDF_sorted(.,by = "pValue") ->
    x_n1a_n1b
  x_n1a_n1b = x_n1a_n1b[,c("gene_id","symbol","Chromosome","pValue","log2FoldChange")]
  x_n1a_n1b %>% 
    merge(.,x_n0_n1b,by = "gene_id") %>% 
    merge(.,x_n0_n1a,by = "gene_id") ->
    x;rm(x_n0_n1a,x_n0_n1b,x_n1a_n1b)
  save(x,file = "./data/x")
}
subset(x,
       (pValue_n0_n1a <pset & abs(log2FoldChange_n0_n1a) > logfcset)) ->
  x_s_1
subset(x,
       (pValue_n0_n1b <pset & abs(log2FoldChange_n0_n1b) > logfcset)) ->
  x_s_2
subset(x,
       (pValue <pset & abs(log2FoldChange) > logfcset)) ->
  x_s_3
subset(x,
       (pValue_n0_n1a <pset & abs(log2FoldChange_n0_n1a) > logfcset) | 
         (pValue_n0_n1b < pset & abs(log2FoldChange_n0_n1b) > logfcset) | 
         (pValue < pset & abs(log2FoldChange) > logfcset)) ->
  x_s
#write.csv(x,file = "./data/x.csv",row.names = F)
#write.csv(x_s,file = "./data/x_s.csv",row.names = F)

getFig2 <- function(){
  pdf(paste0("./data/Fig2c_n1a2n1b.v2",".pdf"),6,5)
  getMP(x = x_s,SNP = "symbol",CHR = "Chromosome",pV = "pValue",Nautosome = 22)
  graphics.off()
  pdf(paste0("./data/Fig2c_n1b",".pdf"),6,5)
  getMP(x = x_s,SNP = "symbol",CHR = "Chromosome",pV = "pValue_n0_n1b",Nautosome = 22)
  graphics.off()
  pdf(paste0("./data/Fig2c_n1a",".pdf"),6,5)
  getMP(x = x_s,SNP = "symbol",CHR = "Chromosome",pV = "pValue_n0_n1a",Nautosome = 22)
  graphics.off()
  
  
  
  fo = merge(expr_tumor_tpm_mean,x_s,by = "symbol")
  ifelse(fo$pValue_n0_n1a < pset & fo$log2FoldChange_n0_n1a > logfcset,1,
         ifelse(fo$pValue_n0_n1a < pset & fo$log2FoldChange_n0_n1a < -logfcset,1,0)) ->
    fo$tag1
  ifelse(fo$pValue_n0_n1b < pset & fo$log2FoldChange_n0_n1b > logfcset,1,
         ifelse(fo$pValue_n0_n1b < pset & fo$log2FoldChange_n0_n1b < -logfcset,1,0)) ->
    fo$tag2
  ifelse(fo$pValue < pset & fo$log2FoldChange > logfcset,1,
         ifelse(fo$pValue < pset & fo$log2FoldChange < -logfcset,1,0)) ->
    fo$tag3
  fo$tag = paste0(fo$tag1,fo$tag2,fo$tag3)
  table(fo$tag1);table(fo$tag2);table(fo$tag3)
  fo = subset(fo,tag != "000")
  table(fo$tag1)
  subset(fo,tag != "111" & pValue < 0.01 & pValue_n0_n1a < 0.01 & pValue_n0_n1b < 0.01) ->
    fo_p001
  subset(expr_tumor_tpm,symbol %in% fo_p001$symbol) -> expr_tumor_tpm_p001
  melt(expr_tumor_tpm_p001) -> expr_tumor_tpm_p001
  expr_tumor_tpm_p001$variable = substr(expr_tumor_tpm_p001$variable,1,12)
  merge(expr_tumor_tpm_p001,p_info,by.x = "variable",by.y = "PID") -> 
    expr_tumor_tpm_p001
  expr_tumor_tpm_p001 = expr_tumor_tpm_p001[,c(2,3,12)]
  
  library(ggridges)
  library(scales)
  data = expr_tumor_tpm_p001
  pdf(paste0("./data/Figs2",".pdf"),6,6)
  ggplot(data, aes(
    x = value,
    y = symbol,
    fill = stat(x),                  
    group = interaction(symbol, path_n) 
  )) +
    geom_density_ridges_gradient(
      alpha = 0.8,
      scale = 0.9,
      color = "white",
      quantile_lines = TRUE,
      quantiles = 2,
      size = 0.5
    ) +
    scale_fill_gradientn(
      name = "Value Level",
      colors = viridis::viridis(3, option = "C"),
      limits = c(-10, 100)
    ) +
    theme_minimal() + 
    theme(
      strip.background = element_rect(fill = "#34495E"),
      strip.text = element_text(color = "white"),
      legend.position = "none"
    ) +
    labs(
      x = "Gene Expression Value",
      y = "Symbol"
    )+
    xlim(c(-10,100))
  graphics.off()
  
  #pdf(paste0(,".pdf"),6,4)
  library(UpSetR)
  library(plyr)
  library(gridExtra)
  library(grid)
  library(scales)
  upset(fo,order.by = "freq",
        sets = c("tag1","tag2","tag3"),
        keep.order = T,
        matrix.color = "purple",
        main.bar.color = hue_pal()(7),
        mainbar.y.label = "DEGs overlapping counts",
        sets.bar.color = "gray",
        sets.x.label = "DEGs size",
        att.color = "black",
        shade.color = "gray")
  #graphics.off()
  #pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,4)
  library(ggtern)
  library(dplyr)
  softmax <- function(x){exp(x)/sum(exp(x))}
  fo2 = fo[,c("symbol","n0","n1a","n1b","tag")]
  for (i in 1:nrow(fo2)) {#i = 1
    if(mean(as.numeric(fo2[i,2:4])) < 30){fo2[i,2:4] = round(softmax(fo2[i,2:4]),5)}
  }
  fo %>%
    ggtern(.,
           aes(x = n0,
               y = n1b,
               z = n1a,
               color=tag,
               fill = tag))+
    #geom_confidence_tern()+
    geom_hex_tern(binwidth=0.01,alpha=0.4)+
    theme_rgbw()+
    geom_point()
  #graphics.off()
  
  expr_tumor_tpm[,colnames(expr_tumor_tpm) %in% 
                   c("symbol",paste0(c(n0,n1a,n1b),".01A"))] ->
    expr_hm
  selnames = subset(fo,tag == "111")$symbol
  expr_hm = expr_hm[expr_hm$symbol %in% selnames,]
  getHM2 <- function(x,genename){
    library(pheatmap)
    library(ggplot2)
    rownames(x) = x[,1]
    x[,-1] -> x
    pheatmap(mat = x,show_colnames = F,
             color=colorRampPalette(c("purple","gray90","red"))(100),
             fontsize_row=6,
             cluster_cols = T,
             scale="row",cellheight=6.2,cellwidth=1.1) ->
      p
    return(p)
  }
  pdf(paste0("./data/Fig2c.v100",".pdf"),6,5)
  getHM2(x = expr_hm,genename = "symbol")
  graphics.off()
  
  
  se <- function(x) sqrt(var(x)/length(x)) 
  
  net2c = getPPInet2c(scoreSET = scoreSET,org = "hs")
  UPs = subset(x,pValue_n0_n1a < pset & log2FoldChange_n0_n1a > logfcset)$symbol
  DNs = subset(x,pValue_n0_n1a < pset & log2FoldChange_n0_n1a < -logfcset)$symbol
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),9,9)
  getPPInet2c2V(UPs = UPs,DNs = DNs,net2c = net2c,showLegend = F)
  graphics.off()
  UPs = subset(x,pValue_n0_n1b < pset & log2FoldChange_n0_n1b > logfcset)$symbol
  DNs = subset(x,pValue_n0_n1b < pset & log2FoldChange_n0_n1b < -logfcset)$symbol
  getPPInet2c2V(UPs = UPs,DNs = DNs,net2c = net2c,showLegend = T)
  UPs = subset(x_n1a_n1b,pValue < pset & log2FoldChange > logfcset)$symbol
  DNs = subset(x_n1a_n1b,pValue < pset & log2FoldChange < -logfcset)$symbol
  getPPInet2c2V(UPs = UPs,DNs = DNs,net2c = net2c,showLegend = T)
  
  library(clusterProfiler)
  load("gsea_n1a");load("gsea_n1b");load("gsea_n1a_n1b");getsave <- function(){
    library(clusterProfiler)
    getT2G(org = "hs") -> t2g
    getLST(dfexp = x,log2fc = "log2FoldChange_n0_n1a",gene_id = "gene_id") -> lst
    getGSEA(LST = lst,t2g = t2g,org = "hs",getRS = F) -> gsea_n1a
    getLST(dfexp = x,log2fc = "log2FoldChange_n0_n1b",gene_id = "gene_id") -> lst
    getGSEA(LST = lst,t2g = t2g,org = "hs",getRS = F) -> gsea_n1b
    getLST(dfexp = x,log2fc = "pValue",gene_id = "gene_id") -> lst
    getGSEA(LST = lst,t2g = t2g,org = "hs",getRS = F) -> gsea_n1a_n1b
    save(gsea_n1a,file = "gsea_n1a");rm(gsea_n1a)
    save(gsea_n1b,file = "gsea_n1b");rm(gsea_n1b)
    save(gsea_n1a_n1b,file = "gsea_n1a_n1b");rm(gsea_n1a_n1b)
  }
  
  TOPN = 5
  gs = gsea_n1a
  gs = subset(gs,EnrichType %in% c("GOBP","GOCC","GOMF","KEGG"))
  subset(gs,status == "UP" & EnrichType == "GOBP")[1:TOPN,1:16] %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "GOCC")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "GOMF")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOBP")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOCC")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOMF")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "KEGG")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "KEGG")[1:TOPN,1:16]) ->
    TOPN
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,6)
  getClev(enrichDF = TOPN,split_by = c("EnrichType","status"))
  graphics.off()
  
  TOPN = 5
  gs = gsea_n1b
  gs = subset(gs,EnrichType %in% c("GOBP","GOCC","GOMF","KEGG"))
  subset(gs,status == "UP" & EnrichType == "GOBP")[1:TOPN,1:16] %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "GOCC")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "GOMF")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOBP")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOCC")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOMF")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "KEGG")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "KEGG")[1:TOPN,1:16]) ->
    TOPN
  getClev(enrichDF = TOPN,split_by = c("EnrichType","status"))
  
  TOPN = 5
  gs = gsea_n1a_n1b
  gs = subset(gs,EnrichType %in% c("GOBP","GOCC","GOMF","KEGG"))
  subset(gs,status == "UP" & EnrichType == "GOBP")[1:TOPN,1:16] %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "GOCC")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "GOMF")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOBP")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOCC")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "DN" & EnrichType == "GOMF")[1:TOPN,1:16]) %>% 
    rbind(.,subset(gs,status == "UP" & EnrichType == "KEGG")[1:TOPN,1:16])  ->
    TOPN
  getClev(enrichDF = TOPN,split_by = c("EnrichType","status"))
  
  clipr::write_clip(TOPN[1:15,2])
  clipr::write_clip(TOPN[16:30,2])
}




load("g_wgcna");getWGCNA <- function(expr,SP01grp,TopPerc,Nsiggrp){
  TopPerc = 0.2
  #Nsiggrp = 3
  expr = expr_tumor
  cls_wg_info = data.frame(name = colnames(expr))
  p_info[,c("PID","path_n")] -> cls_wg_info2
  cls_wg_info2$PID = paste0(cls_wg_info2$PID,".01A")
  merge(cls_wg_info,cls_wg_info2,by.x = "name",by.y = "PID") ->
    mm;SP01grp = mm$path_n
  rownames(expr) = expr[,1]
  expr = expr[,-1]
  library(WGCNA)
  library(dplyr)
  options(stringsAsFactors = FALSE)
  allowWGCNAThreads()
  SP01color = numbers2colors(SP01grp, signed = FALSE)
  nS = ncol(expr)
  
  #Screen genes with top variance over samples
  id = apply(X = expr,MARGIN = 1,FUN = var)
  expr = expr[which(percent_rank(id) > 1-TopPerc),]#TopPerc = 0.2
  expr = as.data.frame(t(expr))
  
  ##Cluster Samples
  gsg = goodSamplesGenes(expr, verbose = 8)
  if(gsg$allOK == T){
    powers = c(c(1:10), seq(from = 12, to=30, by=2))
    sft = pickSoftThreshold(data = expr,
                            powerVector=powers,
                            networkType="unsigned",
                            verbose=5)
    power = sft$powerEstimate
    
    
    net = blockwiseModules(datExpr = expr,
                           power = power,
                           miSPnameoduleSize = 30,
                           maxBlockSize = 10000,
                           TOMType = "unsigned",
                           corType = "pearson",
                           maxPOutliers=1,
                           reassignThreshold = 0,
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           pamRespectsDendro = FALSE,
                           saveTOMs=F,
                           verbose = 3)
    moduleColors = labels2colors(net$colors)
    plotDendroAndColors(dendro = net$dendrograms[[1]],
                        colors = moduleColors[net$blockGenes[[1]]],
                        groupLabels = "Module colors",
                        dendroLabels = FALSE,hang = 0.03,
                        appiGuide = TRUE,guideHang = 0.05)
    #Cluster Genes Over
    # Module-Trait Cor within a heatmap plot
    MEs = moduleEigengenes(expr, moduleColors)$eigengenes
    MEs = orderMEs(MEs)
    moduleTraitCor = cor(x = MEs,y = SP01grp,use = "p")
    moduleTraitCorP = corPvalueStudent(cor = moduleTraitCor,nSamples = nS)
    textMatrix = paste('P = ',signif(moduleTraitCorP, 1), "\n(",
                       signif(moduleTraitCor, 2), ")")
    dim(textMatrix) = dim(moduleTraitCor)
    
    pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),9,6)
    par(mar = c(8, 8, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = "X",
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = colorRampPalette(c("blue","white","#FF6369"))(100),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    graphics.off()
    
    
    apply(data.frame(moduleTraitCorP < 0.05),2,FUN = "+") %>%
      apply(., 1, sum) ->
      kp
    kp = substr(names(kp[kp >= 1]),3,99)
    
    fo = c()
    for (i in 1:length(kp)){
      addgen = names(expr)[moduleColors == kp[i]]
      addcor = rep(kp[i],length(addgen))
      names(addcor) = addgen
      fo = c(fo,addcor)
    }
  }
  
  data.frame(mols = names(fo),colors = fo) -> fo
  fo -> g_wgcna
  rm(g_wgcna)
}
merge(x,g_wgcna,by.x = "gene_id",by.y = "mols") ->
  g_wgcna
getMP2 <- function(x,SNP,CHR,pV,Nautosome){
  #x = g_wgcna,SNP = "symbol",CHR = "Chromosome",pV = "pValue",Nautosome = 22
  x = x[,c(SNP,CHR,pV,c("colors"))]
  colnames(x) = c("SNP","CHR","PValue","colors")
  x = x[is.na(x$SNP)==F,]
  x %<>%
    subset(.,CHR %in% c(1:Nautosome))
  
  HighLights = c()
  for (i in c(1:Nautosome)) {
    x %>%
      subset(.,CHR == i) %>%
      subset(.,PValue == min(PValue)) %>%
      .$SNP ->
      y
    HighLights = c(HighLights,y)
  }
  
  x$CHR = factor(as.numeric(x$CHR))
  x$label <- x$SNP
  x$label[!(x$label %in% HighLights)] <- ""
  
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,4)
  ggplot(x,aes(x = CHR,y=-log10(PValue),color=colors,
               label = label)) +
    geom_jitter() +
    ggrepel::geom_label_repel(max.overlaps = 100) +
    theme_bw() +
    theme(legend.position = "none")
  graphics.off()
}

library(clusterProfiler);t2g = getT2G(org = "hs")
lst1 = getLST(dfexp = x,log2fc = "log2FoldChange_n0_n1a",gene_id = "gene_id")
lst2 = getLST(dfexp = x,log2fc = "log2FoldChange_n0_n1b",gene_id = "gene_id")
lst3 = getLST(dfexp = x,log2fc = "log2FoldChange",gene_id = "gene_id")
gs1 = getGSEA(LST = lst1,org = "hs",t2g = t2g,getRS = F)
gs1$prank_1 = dplyr::percent_rank(gs1$enrichmentScore)
gs2 = getGSEA(LST = lst2,org = "hs",t2g = t2g,getRS = F)
gs2$prank_2 = dplyr::percent_rank(gs2$enrichmentScore)
gs3 = getGSEA(LST = lst3,org = "hs",t2g = t2g,getRS = F)
gs3$prank_3 = dplyr::percent_rank(gs3$enrichmentScore)
getFig3c <- function(){
  merge(merge(x = gs1[,c("Description","prank_1")],
              y = gs2[,c("Description","prank_2")],
              all = T,
              by = "Description"),
        gs3[,c("Description","prank_3","EnrichType")],
        all = T,
        by = "Description") ->
    fo
  fo = fo[abs(fo$prank_1 - 0.5) > 0.45 | 
            abs(fo$prank_2 - 0.5) > 0.45 | 
            abs(fo$prank_3 - 0.5) > 0.45,]
  fo$prank_1  = dplyr::percent_rank(fo$prank_1)
  fo$prank_2  = dplyr::percent_rank(fo$prank_2)
  fo$prank_3  = dplyr::percent_rank(fo$prank_3)
  fo$Description = stringr::str_to_sentence(substr(fo$Description,6,100))
  fo$alpha = abs(apply(fo[,2:4], 1, max) - apply(fo[,2:4], 1, min) - 0.5) + 0.5
  
  library(ggtern)
  fo$Description[fo$prank_1 < 0.008 | fo$prank_2 < 0.008 |fo$prank_3 < 0.008 |
                   fo$prank_1 > 0.992 | fo$prank_2 > 0.992 |fo$prank_3 > 0.992] ->
    fo_name
  subset(fo,Description %in% fo_name) -> fod
  subset(fo,!(Description %in% fo_name)) -> fo
  fo$Description = NA
  rbind(fo,fod) -> fo
  #pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),12,8)
  fo %>%
    ggtern(.,
           aes(x = prank_1,
               y = prank_3,
               z = prank_2,
               color = Description,
               alpha = alpha,
               shape = EnrichType))+
    stat_density_tern() +
    #geom_confidence_tern()+
    #geom_hex_tern(binwidth=0.05,alpha=0.2)+
    theme_rgbw()+
    geom_point(size = 2)+
    labs(title = "Ternary Plot", x = "N0 ->N1a", y = "N1a ->N1b", z = "N0 ->N1b") +
    theme(legend.position = "none")
  graphics.off()
}


TOPN = 4
subset(gs1,status == "UP" & EnrichType == "KEGG")[1:TOPN,2] %>% 
  c(.,subset(gs1,status == "DN" & EnrichType == "KEGG")[1:TOPN,2])%>% 
  c(.,subset(gs1,status == "UP" & EnrichType == "GOBP")[1:TOPN,2])%>% 
  c(.,subset(gs1,status == "DN" & EnrichType == "GOBP")[1:TOPN,2])%>% 
  c(.,subset(gs1,status == "UP" & EnrichType == "GOCC")[1:TOPN,2])%>% 
  c(.,subset(gs1,status == "DN" & EnrichType == "GOCC")[1:TOPN,2])%>% 
  c(.,subset(gs1,status == "UP" & EnrichType == "GOMF")[1:TOPN,2])%>% 
  c(.,subset(gs1,status == "DN" & EnrichType == "GOMF")[1:TOPN,2]) ->
  TOPN
TOPN_df1 = subset(gs1,Description %in% TOPN)
TOPN_df1_up = subset(TOPN_df1,status %in% c("UP"))
TOPN_df1_dn = subset(TOPN_df1,status %in% c("DN"))
TOPN = 4
subset(gs2,status == "UP" & EnrichType == "KEGG")[1:TOPN,2] %>% 
  c(.,subset(gs2,status == "DN" & EnrichType == "KEGG")[1:TOPN,2])%>% 
  c(.,subset(gs2,status == "UP" & EnrichType == "GOBP")[1:TOPN,2])%>% 
  c(.,subset(gs2,status == "DN" & EnrichType == "GOBP")[1:TOPN,2])%>% 
  c(.,subset(gs2,status == "UP" & EnrichType == "GOCC")[1:TOPN,2])%>% 
  c(.,subset(gs2,status == "DN" & EnrichType == "GOCC")[1:TOPN,2])%>% 
  c(.,subset(gs2,status == "UP" & EnrichType == "GOMF")[1:TOPN,2])%>% 
  c(.,subset(gs2,status == "DN" & EnrichType == "GOMF")[1:TOPN,2]) ->
  TOPN
TOPN_df2 = subset(gs2,Description %in% TOPN)
TOPN_df2_up = subset(TOPN_df2,status %in% c("UP"))
TOPN_df2_dn = subset(TOPN_df2,status %in% c("DN"))
TOPN = 4
subset(gs3,status == "UP" & EnrichType == "KEGG")[1:TOPN,2] %>% 
  c(.,subset(gs3,status == "DN" & EnrichType == "KEGG")[1:TOPN,2])%>% 
  c(.,subset(gs3,status == "UP" & EnrichType == "GOBP")[1:TOPN,2])%>% 
  c(.,subset(gs3,status == "DN" & EnrichType == "GOBP")[1:TOPN,2])%>% 
  c(.,subset(gs3,status == "UP" & EnrichType == "GOCC")[1:TOPN,2])%>% 
  c(.,subset(gs3,status == "DN" & EnrichType == "GOCC")[1:TOPN,2])%>% 
  c(.,subset(gs3,status == "UP" & EnrichType == "GOMF")[1:TOPN,2])%>% 
  c(.,subset(gs3,status == "DN" & EnrichType == "GOMF")[1:TOPN,2]) ->
  TOPN
TOPN_df3 = subset(gs3,Description %in% TOPN)

getClev(enrichDF = TOPN_df1,split_by = c("EnrichType","status"),DESwithPREF = T) -> p1
p1+theme(axis.text.y = element_text(size=8),legend.position = "none") -> p1
getClev(enrichDF = TOPN_df2,split_by = c("EnrichType","status"),F) -> p2
p2+theme(axis.text.y = element_text(size=8),legend.position = "none") -> p2
getClev(enrichDF = TOPN_df3,split_by = c("EnrichType","status"),T) -> p3
p3+theme(axis.text.y = element_text(size=8)) -> p3
library(patchwork)
pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),15,6.5)
p1 | p2 | p3
graphics.off()

selTerms = TOPN[,"Description"]
t2g_s = subset(t2g,TERM %in% selTerms)
getGSEA(LST = lst,org = "hs",t2g = t2g_s,getRS = T) ->
  gs;
shDF = gs[[2]];gs = gs[[1]]
getRS2V(shDF = shDF,gs = gs)
library(ggplot2)
library(ggthemes)
library(reshape2)
N = nrow(shDF)/2
s = shDF[1:N,]
s = melt(s,id.vars = c("gene_id","log2FoldChange","rank"))
s$rank = dplyr::percent_rank(s$rank)
s$variable = selTerms[as.numeric(substr(s$variable,3,10))]
s$value = abs(s$value)
ggplot(s, aes(x = rank,y = value,colour = variable,fill = variable)) +
  #geom_area() +
  #geom_polygon()+
  geom_line() +
  scale_size_area() +
  #coord_radial(start = -0.6 * pi, end = 0.6 * pi,, inner.radius = 0.3)+
  #scale_colour_brewer(palette = "Set1") +
  coord_polar(theta = "x") +
  #geom_segment(aes(x = 0,y = 1,xend = 0.5,yend = 1),colour = "#E41A1C") +
  #geom_segment(aes(x = 0.5,y = 1,xend = 1,yend = 1),colour = "#377EB8")+
  ylim(0.3,0.92) + 
  ylab("Running Score")+
  xlab("Position in the Ranked List of Genes") ->
  p;p
clipr::write_clip(subset(TOPN,status == "DN")$Description)
p+guides(fill = guide_legend(title = NULL))+
  theme_minimal()+
  theme(legend.text = element_text(size=4),legend.position = "bottom")

