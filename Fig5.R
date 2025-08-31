
library(magrittr)
library(reshape2)

read.csv("./data/expr.csv") -> expr
round(na.omit(expr)) -> expr
N_value0 = apply(expr,1,function(x) table(x)["0"])
expr = expr[c(N_value0 < (ncol(expr) - 1)/20 | is.na(N_value0) == T),]#g exp>1/20 sp
colnames(expr) = substr(colnames(expr),1,16)
expr %<>% 
  .[,(
    (substr(colnames(.),14,16) %in% c("01A","11A")) &
      (substr(colnames(.),1,12) %in% p_info$PID)) | 
      colnames(.) == "entrez"]
cons = colnames(expr)
cons[substr(cons,14,16) == "01A"] -> exps
cons[substr(cons,14,16) == "11A"] -> cons
expr_tumor = expr[,c("entrez",exps)]

library(Seurat)
df1 = expr_tumor
df1$entrez <- g_info$symbol[match(df1$entrez,g_info$gene_id)]
df1 = na.omit(df1)
getUmap <- function(mat){
  mat %>% CreateSeuratObject(.) -> 
    se0
  se0 %<>% #Normalization
    NormalizeData(., normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(., selection.method = "vst", nfeatures = 2000)
  se0 %<>%#scaling to 0,1 distribution to PCA
    ScaleData(., features = rownames(.)) %>%
    RunPCA(., features = VariableFeatures(.))
  pct <- (se0[["pca"]]@stdev)/(sum(se0[["pca"]]@stdev))
  ProperN4PC = which(cumsum(pct) > 0.9 & pct < 0.05)[1]
  
  #top10 <- head(VariableFeatures(se0), 10)
  #plot1 <- VariableFeaturePlot(se0)
  #LabelPoints(plot = plot1, points = top10, repel = TRUE)
  
  se0 %>% 
    FindNeighbors(., dims = 1:ProperN4PC) %>%
    FindClusters(.) %>%
    RunUMAP(.,dims = 1:ProperN4PC) ->
    se.umap
  #DimHeatmap(se.umap, dims = 1:6, cells = 500, balanced = TRUE)
  #DimPlot(se.umap, reduction = "umap")
  
  return(se.umap)
}
getUmap(mat = expr_tumor) -> 
  seumap
DimHeatmap(seumap, dims = 1:6, balanced = TRUE)
DimPlot(seumap, reduction = "umap")
seumap@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% 
  cbind(tx = seumap@meta.data$seurat_clusters) %>% 
  mutate(.,lb = substr(rownames(.),9,12)) %>% 
  mutate(.,PID = substr(rownames(.),1,12)) %>% 
  merge(.,p_info,by = "PID")->
  seump_xy;data.frame(SID = names(Idents(seumap)),
                      cls = Idents(seumap),
                      row.names = NULL) %>% 
  subset(.,
         substr(SID,1,12) %in% p_info$PID) -> 
  cls_info
seump_xy <- seump_xy %>% mutate(path_n = as.factor(path_n))
library(ggrepel)
pdf(paste0("./data/Fig5a.v5",".pdf"),6,5)
ggplot(seump_xy, aes(x = umap_1, y = umap_2, color = tx, shape = path_n)) + 
  geom_point(size = 3) +  
  geom_text_repel(
    aes(label = lb), 
    force = 3, 
    size = 2, 
    max.overlaps = 20,
    fontface = "italic", 
    arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"), 
    nudge_x = 0.1, 
    nudge_y = 0.06
  ) +
  theme_pubclean() +
  labs(color = "Cluster", shape = "N stage")  # 优化图例标题

graphics.off()
cls_info$SID = substr(cls_info$SID,1,12)
merge(p_info,cls_info,by.x = "PID",by.y = "SID") -> 
  pat_info_new
tab_c2n = data.frame(table(pat_info_new$cls,pat_info_new$path_n))
reshape2::dcast(tab_c2n,formula = Var1~Var2) ->
  tab_c2n;tab_c2n_plot = tab_c2n
colnames(tab_c2n_plot) = c("Cluster","N0","N1a","N1b")
tab_c2n_plot = melt(tab_c2n_plot,id.vars = "Cluster")
colnames(tab_c2n_plot) = c("Cluster","N","count")
tab_c2n_plot$Cluster = paste0("cls",tab_c2n_plot$Cluster)
tab_c2n_plot$lb = c(1,1,0,0,1,1,0,0,1,1,0,1)
pdf(paste0("./data/Fig5b",".pdf"),6,5)
ggplot(tab_c2n_plot, aes(x=Cluster, y=N)) +
  geom_point(aes(size=count), shape=21, colour="black", fill="purple") +
  scale_size_area(max_size=20) +
  geom_text(aes(y=as.numeric(N)-0.25, label= paste0(lb,"/",count)), vjust=1,
            colour="grey60", size=4)+
  theme(legend.position = "none")+
  ylab("N stage")
graphics.off()
tab_c2n_pct = apply(tab_c2n[,!(colnames(tab_c2n) == "Var1")],1,sum)
tab_c2n_pct = tab_c2n[,!(colnames(tab_c2n) == "Var1")]/tab_c2n_pct
tab_c2n_pct = round(tab_c2n_pct,3)
tab_c2n_pct$cls = paste0("cls--",as.numeric(rownames(tab_c2n_pct)) -1)
pat_info_new[,c("times","patient.vital_status","path_m","path_n","path_t","cls")] -> df
library(survival)
sfit <- survfit(Surv(times, patient.vital_status)~cls+path_n, data=df)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)

#N0 low risk -> high risk
cons = intersect(n0,cls_info$SID[cls_info$cls == 2])
exps = intersect(n0,cls_info$SID[cls_info$cls == 0 | cls_info$cls == 1])

load("drug_lr_hr_n0");getsave <- function(){
  load("save_x_lr_hr_n0");getsave <- function(){
    getX(expr = expr_tumor,
         gname = "entrez",
         cons = intersect(paste0(cons,".01A"),colnames(expr_tumor)),
         exps = intersect(paste0(exps,".01A"),colnames(expr_tumor)),
         method = "rc-deseq2") %>% 
      merge(g_info,.,by.x = "gene_id",by.y = "entrez") %>% 
      getDF_sorted(.,by = "pValue") ->
      x_lr_hr_n0
    x_lr_hr_n0 |> colnames()
    x_lr_hr_n0 = x_lr_hr_n0[,c("gene_id","pValue","log2FoldChange")]
    save(x_lr_hr_n0,file = "save_x_lr_hr_n0");rm(x_lr_hr_n0)
  }
  PSET2 = 0.001;LOGFCSET2 = 2
  merge(x_lr_hr_n0,g_info,by = "gene_id") -> x_lr_hr
  UPs = subset(x_lr_hr,pValue < PSET2 & log2FoldChange > LOGFCSET2)$symbol
  UPs_hid = subset(g_info,symbol %in% stringr::str_to_upper(UPs))$gene_id
  DNs = subset(x_lr_hr,pValue < PSET2 & log2FoldChange < -LOGFCSET2)$symbol
  DNs_hid = subset(g_info,symbol %in% stringr::str_to_upper(DNs))$gene_id
  
  
  
  if("getDistancePre" != 0){
    drug = read.csv("drugs.csv")
    getPPInet2c(scoreSET = scoreSET,org = "hs") -> ppi2c
    BG = unique((c(ppi2c[,1],ppi2c[,2])))
    graph_from_data_frame(ppi2c,directed = F) -> ppi2c_graph
    ppis = distances(ppi2c_graph, mode="out")
    ppid = degree(ppi2c_graph)}
  if("GetDistance" != 0){
    dd = rep("NA",nrow(drug))
    for (tag in 1:nrow(drug)) {
      #tag = 11
      drug[tag,2] %>%
        stringr::str_split(.,",") %>%
        unlist(.) ->
        tg
      getDT2DD(tg = tg,dg = c(UPs,DNs),ppi2c = ppi2c) -> dd[tag]
    }
    drug$distance = as.numeric(dd)
    drug = drug[,c("Drugs","distance","Target")]
    drug$percDistance = dplyr::percent_rank(drug$distance)
  }
  drug -> drug_lr_hr_n0_distance_v250518
  save(drug_lr_hr_n0_distance_v250518,file = "drug_lr_hr_n0_distance_v250518")
  rm(drug_lr_hr_n0_distance)
  ??percent_rank
  
  if("getL1000pre" != 0){
    library(cmapR)
    LI = "~/Downloads/BIGDATA//GSE92742_Broad_LINCS_sig_info.txt"
    col_meta <- read.delim(LI, sep="\t", stringsAsFactors=F)
    cl = names(table(col_meta$cell_id))[1:76]
  }
  if(0 != "get L1000"){
    t1 = Sys.time()
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    getL1000(UPs_hid,DNs_hid,cellline = cl) -> l1000
    print(Sys.time() - t1)
  }
  l1000 -> drug_lr_hr_n0
  save(drug_lr_hr_n0,file = "drug_lr_hr_n0")
  rm(drug_lr_hr_n0)
}



#N1a low risk -> high risk
cons = intersect(n1a,cls_info$SID[cls_info$cls == 2])
exps = intersect(n1a,cls_info$SID[cls_info$cls == 0 | cls_info$cls == 1])

load("drug_lr_hr_n1a");getsave <- function(){
  load("save_x_lr_hr_n1a");getsave <- function(){
    getX(expr = expr_tumor,
         gname = "entrez",
         cons = intersect(paste0(cons,".01A"),colnames(expr_tumor)),
         exps = intersect(paste0(exps,".01A"),colnames(expr_tumor)),
         method = "rc-deseq2") %>% 
      merge(g_info,.,by.x = "gene_id",by.y = "entrez") %>% 
      getDF_sorted(.,by = "pValue") ->
      x_lr_hr_n1a
    x_lr_hr_n1a |> colnames()
    x_lr_hr_n1a = x_lr_hr_n1a[,c("gene_id","pValue","log2FoldChange")]
    save(x_lr_hr_n1a,file = "save_x_lr_hr_n1a");rm(x_lr_hr_n1a)
  }
  PSET2 = 0.001;LOGFCSET2 = 2
  merge(x_lr_hr_n1a,g_info,by = "gene_id") -> x_lr_hr
  UPs = subset(x_lr_hr,pValue < PSET2 & log2FoldChange > LOGFCSET2)$symbol
  UPs_hid = subset(g_info,symbol %in% stringr::str_to_upper(UPs))$gene_id
  DNs = subset(x_lr_hr,pValue < PSET2 & log2FoldChange < -LOGFCSET2)$symbol
  DNs_hid = subset(g_info,symbol %in% stringr::str_to_upper(DNs))$gene_id
  
  
  
  if("getDistancePre" != 0){
    drug = read.csv("drugs.csv")
    getPPInet2c(scoreSET = scoreSET,org = "hs") -> ppi2c
    BG = unique((c(ppi2c[,1],ppi2c[,2])))
    graph_from_data_frame(ppi2c,directed = F) -> ppi2c_graph
    ppis = distances(ppi2c_graph, mode="out")
    ppid = degree(ppi2c_graph)}
  if("GetDistance" != 0){
    dd = rep("NA",nrow(drug))
    for (tag in 1:nrow(drug)) {
      #tag = 11
      drug[tag,2] %>%
        stringr::str_split(.,",") %>%
        unlist(.) ->
        tg
      getDT2DD(tg = tg,dg = c(UPs,DNs),ppi2c = ppi2c) -> dd[tag]
    }
    drug$distance = as.numeric(dd)
    drug = drug[,c("Drugs","distance","Target")]
    drug$percDistance = dplyr::percent_rank(drug$distance)
  }
  drug -> drug_lr_hr_n1a_distance_v250518
  save(drug_lr_hr_n1a_distance_v250518,file = "drug_lr_hr_n1a_distance_v250518")
  rm(drug_lr_hr_n1a_distance)
  
  
  
  if("getL1000pre" != 0){
    library(cmapR)
    LI = "~/Downloads/BIGDATA//GSE92742_Broad_LINCS_sig_info.txt"
    col_meta <- read.delim(LI, sep="\t", stringsAsFactors=F)
    cl = names(table(col_meta$cell_id))[1:76]
  }
  if(0 != "get L1000"){
    t1 = Sys.time()
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    getL1000(UPs_hid,DNs_hid,cellline = cl) -> l1000
    print(Sys.time() - t1)
  }
  l1000 -> drug_lr_hr_n1a
  save(drug_lr_hr_n1a,file = "drug_lr_hr_n1a")
  rm(drug_lr_hr_n1a)
}



#N1b low risk -> high risk
cons = intersect(n1b,cls_info$SID[cls_info$cls == 2])
exps = intersect(n1b,cls_info$SID[cls_info$cls == 0 | cls_info$cls == 1])

load("drug_lr_hr_n1b");getsave <- function(){
  load("save_x_lr_hr_n1b");getsave <- function(){
    getX(expr = expr_tumor,
         gname = "entrez",
         cons = intersect(paste0(cons,".01A"),colnames(expr_tumor)),
         exps = intersect(paste0(exps,".01A"),colnames(expr_tumor)),
         method = "rc-deseq2") %>% 
      merge(g_info,.,by.x = "gene_id",by.y = "entrez") %>% 
      getDF_sorted(.,by = "pValue") ->
      x_lr_hr_n1b
    x_lr_hr_n1b |> colnames()
    x_lr_hr_n1b = x_lr_hr_n1b[,c("gene_id","pValue","log2FoldChange")]
    save(x_lr_hr_n1b,file = "save_x_lr_hr_n1b");rm(x_lr_hr_n1a)
  }
  PSET2 = 0.001;LOGFCSET2 = 2
  merge(x_lr_hr_n1b,g_info,by = "gene_id") -> x_lr_hr
  UPs = subset(x_lr_hr,pValue < PSET2 & log2FoldChange > LOGFCSET2)$symbol
  UPs_hid = subset(g_info,symbol %in% stringr::str_to_upper(UPs))$gene_id
  DNs = subset(x_lr_hr,pValue < PSET2 & log2FoldChange < -LOGFCSET2)$symbol
  DNs_hid = subset(g_info,symbol %in% stringr::str_to_upper(DNs))$gene_id
  
  
  
  if("getDistancePre" != 0){
    drug = read.csv("drugs.csv")
    getPPInet2c(scoreSET = scoreSET,org = "hs") -> ppi2c
    BG = unique((c(ppi2c[,1],ppi2c[,2])))
    graph_from_data_frame(ppi2c,directed = F) -> ppi2c_graph
    ppis = distances(ppi2c_graph, mode="out")
    ppid = degree(ppi2c_graph)}
  if("GetDistance" != 0){
    dd = rep("NA",nrow(drug))
    for (tag in 1:nrow(drug)) {
      #tag = 11
      drug[tag,2] %>%
        stringr::str_split(.,",") %>%
        unlist(.) ->
        tg
      getDT2DD(tg = tg,dg = c(UPs,DNs),ppi2c = ppi2c) -> dd[tag]
    }
    drug$distance = as.numeric(dd)
    drug = drug[,c("Drugs","distance","Target")]
    drug$percDistance = dplyr::percent_rank(drug$distance)
  }
  drug -> drug_lr_hr_n1b_distance_v250518
  save(drug_lr_hr_n1b_distance_v250518,file = "drug_lr_hr_n1b_distance_v250518")
  rm(drug_lr_hr_n1b_distance)
  
  
  
  
  if("getL1000pre" != 0){
    library(cmapR)
    LI = "~/Downloads/BIGDATA//GSE92742_Broad_LINCS_sig_info.txt"
    col_meta <- read.delim(LI, sep="\t", stringsAsFactors=F)
    cl = names(table(col_meta$cell_id))[1:76]
  }
  if(0 != "get L1000"){
    t1 = Sys.time()
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    getL1000(UPs_hid,DNs_hid,cellline = cl) -> l1000
    print(Sys.time() - t1)
  }
  l1000 -> drug_lr_hr_n1b
  save(drug_lr_hr_n1b,file = "drug_lr_hr_n1b")
  rm(drug_lr_hr_n1b)
}

write.csv(drug_lr_hr_n0,file = "drug_lr_hr_n0.csv",row.names = F)
write.csv(drug_lr_hr_n1a,file = "drug_lr_hr_n1a.csv",row.names = F)
write.csv(drug_lr_hr_n1b,file = "drug_lr_hr_n1b.csv",row.names = F)

getmain4drugplot <- function(){
  #rm(list = ls());
  t=Sys.time();
  c4c="~/Library/CloudStorage/OneDrive-Personal/code4code/"
  #install.packages(paste0(c4c,"U/xbox/"),repos = NULL,type = "source")
  library(xbox)
  library(magrittr)
  library(reshape2)
  
  pset = 0.05;nesset = -1.64
  
  #N0 low risk -> high risk
  load("drug_lr_hr_n0_distance_v250518")
  d_dis = drug_lr_hr_n0_distance_v250518;rm(drug_lr_hr_n0_distance_v250518)
  d_dis_sig = subset(d_dis,percDistance < 0.01)
  
  load("drug_lr_hr_n0");
  d_l1000_all = drug_lr_hr_n0;rm(drug_lr_hr_n0)
  data.frame(Drugs = names(table(d_l1000_all$pert_iname)),
             Ncond = as.numeric(table(d_l1000_all$pert_iname))) -> 
    Ncond
  d_l1000_sig = subset(d_l1000_all,pV < pset & NES < nesset)
  data.frame(Drugs = names(table(d_l1000_sig$pert_iname)),
             Nsig = as.numeric(table(d_l1000_sig$pert_iname))) -> 
    Nsig
  
  merge(Nsig,Ncond,by = "Drugs") -> N;rm(Nsig);rm(Ncond)
  d =merge(d_dis_sig,N,by = "Drugs")
  d$Category = "Arresting in N0";
  d-> d1;rm(d);rm(d_dis);rm(d_dis_sig);rm(d_l1000_all);rm(d_l1000_sig);rm(N)
  DT::datatable(d1)
  
  
  #N1a low risk -> high risk
  load("drug_lr_hr_n1a_distance_v250518")
  d_dis = drug_lr_hr_n1a_distance_v250518;rm(drug_lr_hr_n1a_distance_v250518)
  d_dis_sig = subset(d_dis,percDistance < 0.01)
  
  load("drug_lr_hr_n1a")
  d_l1000_all = drug_lr_hr_n1a;rm(drug_lr_hr_n1a)
  data.frame(Drugs = names(table(d_l1000_all$pert_iname)),
             Ncond = as.numeric(table(d_l1000_all$pert_iname))) -> 
    Ncond
  d_l1000_sig = subset(d_l1000_all,pV < pset & NES < nesset)
  data.frame(Drugs = names(table(d_l1000_sig$pert_iname)),
             Nsig = as.numeric(table(d_l1000_sig$pert_iname))) -> 
    Nsig
  
  merge(Nsig,Ncond,by = "Drugs") -> N;rm(Nsig);rm(Ncond)
  d =merge(d_dis_sig,N,by = "Drugs")
  d$Category = "Arresting in n1a";
  d-> d2;rm(d);rm(d_dis);rm(d_dis_sig);rm(d_l1000_all);rm(d_l1000_sig);rm(N)
  DT::datatable(d2)
  
  #N1b low risk -> high risk
  load("drug_lr_hr_n1b_distance_v250518")
  d_dis = drug_lr_hr_n1b_distance_v250518;rm(drug_lr_hr_n1b_distance_v250518)
  d_dis_sig = subset(d_dis,percDistance < 0.01)
  
  load("drug_lr_hr_n1b")
  d_l1000_all = drug_lr_hr_n1b;rm(drug_lr_hr_n1b)
  data.frame(Drugs = names(table(d_l1000_all$pert_iname)),
             Ncond = as.numeric(table(d_l1000_all$pert_iname))) -> 
    Ncond
  d_l1000_sig = subset(d_l1000_all,pV < pset & NES < nesset)
  data.frame(Drugs = names(table(d_l1000_sig$pert_iname)),
             Nsig = as.numeric(table(d_l1000_sig$pert_iname))) -> 
    Nsig
  
  merge(Nsig,Ncond,by = "Drugs") -> N;rm(Nsig);rm(Ncond)
  d =merge(d_dis_sig,N,by = "Drugs")
  d$Category = "Arresting in N1b";
  d-> d3;rm(d);rm(d_dis);rm(d_dis_sig);rm(d_l1000_all);rm(d_l1000_sig);rm(N)
  DT::datatable(d3)
  
  rbind(d1,rbind(d2,d3)) -> d
  d = d[,c(1,2,3,7)]
  DT::datatable(d)
  write.csv(d,file = "Table.drug.csv",row.names = F)
  
}