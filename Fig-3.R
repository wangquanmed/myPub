
getLogisticM <- function(){
  pset = 0.05;logfcset = 1;
  set_g_exp_tag = 1
  
  library(ggplot2)
  library(pROC)
  c(expr_tumor_tpm_mean$n1a > set_g_exp_tag |expr_tumor_tpm_mean$n0 > set_g_exp_tag) %>% 
    expr_tumor_tpm_mean$symbol[.] ->
    g_exp_tag
  
  
  degs = x_s$symbol[x_s$pValue_n0_n1a < pset & abs(x_s$log2FoldChange_n0_n1a) > logfcset]
  degs = intersect(degs,g_exp_tag)
  merge(g_info,expr_tumor_tpm,by = "symbol") ->
    df
  df = df[df$symbol %in% degs,
          colnames(df) %in% c(paste0(c(n0,n1a),".01A"),"symbol")]
  getExpr2Mat <- function(expr,by){
    rn = expr[,by]
    rownames(expr) = rn
    expr = expr[,!(colnames(expr) %in% by)]
    return(expr)
  }
  getExpr2Mat(expr = df,by = "symbol") ->
    df
  df = data.frame(t(df));
  shapiro.test(df[,5])$p.value
  df = df + 0.01
  library(car)
  rn = rownames(df)
  for (i in 1:ncol(df)) {
    df[rn,i] %<>%  bcPower(.,powerTransform(.)$lambda)
  };densityPlot(df[,5])
  shapiro.test(df[,5]);
  
  cbind(data.frame(Outcome = ifelse(merge(data.frame(PID = substr(rownames(df),1,12)),
                                          p_info,by = "PID")$path_n == "1",0,1)),
        df) ->
    df
  
  df[,c(1:3)] %>% 
    reshape2::melt(.,id.vars = "Outcome") %>% 
    ggplot(data = .,mapping = aes(x = variable,y = value,color=factor(Outcome)))+
    geom_jitter()+
    xlab("")
  
  for (i in 2:(ncol(df))) {
    print(roc(df$Outcome,df[,i])$auc)
  }
  auc_set = 0.7
  g_sig = c()
  for (i in 2:(ncol(df))) {
    if(roc(df$Outcome,df[,i])$auc > auc_set){g_sig = c(g_sig,i)}
  }
  
  df = df[,colnames(df)[c(1,g_sig)]]
  
  ggroc(list(Exap1=roc(df$Outcome,as.numeric(df[,2])),
             Exap2=roc(df$Outcome,as.numeric(df[,3]))))
  library(rms)
  outcome = "Outcome ~ "
  lm.sol=lm(formula = as.formula(paste(outcome, 
                                       paste(colnames(df)[2:ncol(df)], collapse= "+"))),
            data=df)
  summary(lm.sol)
  lm.step = step(lm.sol)
  a = summary(lm.step)$coefficients
  b = rownames(a)[a[,4] < 0.1];b = b[b != "(Intercept)"]
  a = as.formula(paste(outcome, paste(b, collapse= "+")))
  
  ddist <- datadist(df)
  options(datadist = 'ddist')
  nom <- nomogram(lrm(formula = a,data = df),
                  lp=T,
                  lp.at = seq(-3,4,by=0.5),
                  fun=function(x) 1/(1+exp(-x)),
                  funlabel = 'Risk of Malignancy',
                  fun.at = c(0.05,seq(0.1,0.9,by=0.1),0.95),
                  conf.int = c(0.1,0.7))
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),6,4)
  plot(nom,
       lplabel = 'Linear Predictor',
       fun.side = c(1,1,1,1,3,1,3,1,1,1,1),
       label.every = 3,
       col.conf = c('purple','grey'),
       conf.space = c(0.1,0.2),
       col.grid = gray(c(0.8,0.95)),
       which='shock')
  graphics.off()
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),4,4)
  ggroc(list(C16orf89=roc(df$Outcome,as.numeric(df[,"C16orf89"])),
             CLEC3B=roc(df$Outcome,as.numeric(df[,"CLEC3B"])),
             USP2=roc(df$Outcome,as.numeric(df[,"USP2"])),
             WFIKKN1=roc(df$Outcome,as.numeric(df[,"WFIKKN1"]))
  ))
  graphics.off()
  
  
  library(ggplot2)
  library(pROC)
  c(expr_tumor_tpm_mean$n1b > set_g_exp_tag |expr_tumor_tpm_mean$n0 > set_g_exp_tag) %>% 
    expr_tumor_tpm_mean$symbol[.] ->
    g_exp_tag
  
  degs = x_s$symbol[x_s$pValue_n0_n1b < pset & abs(x_s$log2FoldChange_n0_n1b) > logfcset]
  degs = intersect(degs,g_exp_tag)
  
  
  merge(g_info,expr_tumor_tpm,by = "symbol") ->
    df
  df = df[df$symbol %in% degs,
          colnames(df) %in% c(paste0(c(n0,n1b),".01A"),"symbol")]
  getExpr2Mat <- function(expr,by){
    rn = expr[,by]
    rownames(expr) = rn
    expr = expr[,!(colnames(expr) %in% by)]
    return(expr)
  }
  getExpr2Mat(expr = df,by = "symbol") ->
    df
  df = data.frame(t(df));
  shapiro.test(df[,5])$p.value
  df = df + 0.01
  library(car)
  rn = rownames(df)
  for (i in 1:ncol(df)) {
    df[rn,i] %<>%  bcPower(.,powerTransform(.)$lambda)
  };densityPlot(df[,5])
  shapiro.test(df[,5]);
  
  cbind(data.frame(Outcome = ifelse(merge(data.frame(PID = substr(rownames(df),1,12)),
                                          p_info,by = "PID")$path_n == "0",0,1)),
        df) ->
    df
  
  df[,c(1:3)] %>% 
    reshape2::melt(.,id.vars = "Outcome") %>% 
    ggplot(data = .,mapping = aes(x = variable,y = value,color=factor(Outcome)))+
    geom_jitter()+
    xlab("")
  
  for (i in 2:(ncol(df))) {
    print(roc(df$Outcome,df[,i])$auc)
  }
  auc_set = 0.7
  g_sig = c()
  for (i in 2:(ncol(df))) {
    if(roc(df$Outcome,df[,i])$auc > auc_set){g_sig = c(g_sig,i)}
  }
  
  df = df[,colnames(df)[c(1,g_sig)]]
  
  ggroc(list(Exap1=roc(df$Outcome,as.numeric(df[,2])),
             Exap2=roc(df$Outcome,as.numeric(df[,3]))))
  library(rms)
  outcome = "Outcome ~ "
  lm.sol=lm(formula = as.formula(paste(outcome, 
                                       paste(colnames(df)[2:ncol(df)], collapse= "+"))),
            data=df)
  summary(lm.sol)
  lm.step = step(lm.sol)
  a = summary(lm.step)$coefficients
  b = rownames(a)[a[,4] < 0.1/10];b = b[b != "(Intercept)"]
  a = as.formula(paste(outcome, paste(b, collapse= "+")))
  
  ddist <- datadist(df)
  options(datadist = 'ddist')
  nom <- nomogram(lrm(formula = a,data = df),
                  lp=T,
                  lp.at = seq(-3,4,by=0.5),
                  fun=function(x) 1/(1+exp(-x)),
                  funlabel = 'Risk of Malignancy',
                  fun.at = c(0.05,seq(0.1,0.9,by=0.1),0.95),
                  conf.int = c(0.1,0.7))
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),6,4)
  plot(nom,
       lplabel = 'Linear Predictor',
       fun.side = c(1,1,1,1,3,1,3,1,1,1,1),
       label.every = 3,
       col.conf = c('purple','grey'),
       conf.space = c(0.1,0.2),
       col.grid = gray(c(0.8,0.95)),
       which='shock')
  graphics.off()
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),4,4)
  ggroc(list(FKBP5=roc(df$Outcome,as.numeric(df[,"FKBP5"])),
             TREM1=roc(df$Outcome,as.numeric(df[,"TREM1"]))
  ))
  graphics.off()
  
  
  library(ggplot2)
  library(pROC)
  c(expr_tumor_tpm_mean$n1a > set_g_exp_tag |expr_tumor_tpm_mean$n1b > set_g_exp_tag) %>% 
    expr_tumor_tpm_mean$symbol[.] ->
    g_exp_tag
  
  degs = x_s$symbol[x_s$pValue < pset & abs(x_s$log2FoldChange) > logfcset]
  degs = intersect(degs,g_exp_tag)
  merge(g_info,expr_tumor_tpm,by = "symbol") ->
    df
  df = df[df$symbol %in% degs,
          colnames(df) %in% c(paste0(c(n1a,n1b),".01A"),"symbol")]
  getExpr2Mat <- function(expr,by){
    rn = expr[,by]
    rownames(expr) = rn
    expr = expr[,!(colnames(expr) %in% by)]
    return(expr)
  }
  getExpr2Mat(expr = df,by = "symbol") ->
    df
  df = data.frame(t(df));
  shapiro.test(df[,5])$p.value
  df = df + 0.01
  library(car)
  rn = rownames(df)
  for (i in 1:ncol(df)) {
    df[rn,i] %<>%  bcPower(.,powerTransform(.)$lambda)
  };densityPlot(df[,5])
  shapiro.test(df[,5]);
  
  cbind(data.frame(Outcome = ifelse(merge(data.frame(PID = substr(rownames(df),1,12)),
                                          p_info,by = "PID")$path_n == "1",0,1)),
        df) ->
    df
  
  df[,c(1:3)] %>% 
    reshape2::melt(.,id.vars = "Outcome") %>% 
    ggplot(data = .,mapping = aes(x = variable,y = value,color=factor(Outcome)))+
    geom_jitter()+
    xlab("")
  
  for (i in 2:(ncol(df))) {
    print(roc(df$Outcome,df[,i])$auc)
  }
  auc_set = 0.7
  g_sig = c()
  for (i in 2:(ncol(df))) {
    if(roc(df$Outcome,df[,i])$auc > auc_set){g_sig = c(g_sig,i)}
  }
  
  df = df[,colnames(df)[c(1,g_sig)]]
  
  ggroc(list(Exap1=roc(df$Outcome,as.numeric(df[,2])),
             Exap2=roc(df$Outcome,as.numeric(df[,3]))))
  library(rms)
  outcome = "Outcome ~ "
  lm.sol=lm(formula = as.formula(paste(outcome, 
                                       paste(colnames(df)[2:ncol(df)], collapse= "+"))),
            data=df)
  summary(lm.sol)
  lm.step = step(lm.sol)
  a = summary(lm.step)$coefficients
  b = rownames(a)[a[,4] < 0.1];b = b[b != "(Intercept)"]
  a = as.formula(paste(outcome, paste(b, collapse= "+")))
  
  ddist <- datadist(df)
  options(datadist = 'ddist')
  nom <- nomogram(lrm(formula = a,data = df),
                  lp=T,
                  lp.at = seq(-3,4,by=0.5),
                  fun=function(x) 1/(1+exp(-x)),
                  funlabel = 'Risk of Malignancy',
                  fun.at = c(0.05,seq(0.1,0.9,by=0.1),0.95),
                  conf.int = c(0.1,0.7))
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),6,4)
  plot(nom,
       lplabel = 'Linear Predictor',
       fun.side = c(1,1,1,1,3,1,3,1,1,1,1),
       label.every = 3,
       col.conf = c('purple','grey'),
       conf.space = c(0.1,0.2),
       col.grid = gray(c(0.8,0.95)),
       which='shock')
  graphics.off()
  pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),4,4)
  ggroc(list(HBA1=roc(df$Outcome,as.numeric(df[,"HBA1"])),
             SERPINE1=roc(df$Outcome,as.numeric(df[,"SERPINE1"]))
  ));graphics.off()
  
  #se.markers=FindAllMarkers(seumap, only.pos = TRUE,min.pct = 0.4,logfc.threshold = 0.4)
  merge(pat_info,sam_info_cls_tumor,by.x = "PID",by.y = "SID") -> 
    pat_info_new
  tab_c2n = data.frame(table(pat_info_new$cls,pat_info_new$path_n))
  reshape2::dcast(tab_c2n,formula = Var1~Var2) ->
    tab_c2n
  tab_c2n_pct = apply(tab_c2n[,!(colnames(tab_c2n) == "Var1")],1,sum)
  tab_c2n_pct = tab_c2n[,!(colnames(tab_c2n) == "Var1")]/tab_c2n_pct
  tab_c2n_pct = round(tab_c2n_pct,3)
  tab_c2n_pct$cls = paste0("cls--",as.numeric(rownames(tab_c2n_pct)) -1)
  
  expr %>% 
    .[,(substr(colnames(.),14,16) %in% c("01A")) & 
        (substr(colnames(.),1,12) %in% pat_info$PID)] ->
    expr_tumor
  library(dplyr)
  se.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) ->
    sm
  Sys.time() - t
  
  library(survival)
  sfit <- survfit(Surv(times, vital_status)~cls, data=pat_info_new)
  print(sfit)
  ggsurvplot(sfit, conf.int=F, pval=TRUE)
  
  if("wgcna"){
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
    library(DoubletFinder)
    library(WGCNA)
    Ud = "~/Library/CloudStorage/OneDrive-Personal/code4code/U/get/"
    Us = dir(Ud)
    e = getE("mm")
    read.csv("~/Library/CloudStorage/OneDrive-Personal/DrugBank/drugs-220421.csv") ->
      drug
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    rm(i);rm(Ud);rm(Us)
  };print(Sys.time() - t1)
  
  
  
  
  
  
  pp$t4 = ifelse(pp$path_t %in% c("t4","t4a"),1,0)
  pp$t3 = ifelse(pp$path_t %in% c("t3"),1,0)
  pp$t2 = ifelse(pp$path_t %in% c("t2"),1,0)
  pp$t1 = ifelse(pp$path_t %in% c("t1","t1a","t1b"),1,0)
  pp$n1 = ifelse(pp$path_n %in% c("n1","n1a","n1b"),1,0)
  pp$n1a = ifelse(pp$path_n %in% c("n1a"),1,0)
  pp$n1b = ifelse(pp$path_n %in% c("n1b"),1,0)
  pp$n0 = ifelse(pp$path_n %in% c("n0"),1,0)
  pp$m1 = ifelse(pp$path_m %in% c("m1"),1,0)
  pp$m0 = ifelse(pp$path_m %in% c("m0"),1,0)
  n0 = pp$sid[pp$n0 == 1]
  n1a = pp$sid[pp$n1a == 1]
  n1b = pp$sid[pp$n1b == 1]
  rowMeans(expr_t[,paste0(n0,".01A")]) -> n0mean
  rowMeans(expr_t[,paste0(n1a,".01A")]) -> n1amean
  rowMeans(expr_t[,paste0(n1b,".01A")]) -> n1bmean;nsum = n0mean+n1amean+n1bmean
  
  
  data.frame(entrez = expr_tumor$entrez,
             n0mean = round(n0mean,3),
             n1amean = round(n1amean,3),
             n1bmean = round(n1bmean,3),
             n0mean_pct = round(n0mean/nsum,3),
             n1amean_pct = round(n1amean/nsum,3),
             n1bmean_pct = round(n1bmean/nsum,3)) ->
    expr_tumor_ln_mean
  
  pp_ln_exact = subset(pp,n1a == 1 | n1b == 1 | n0 == 1)
  expr_ln_exact = expr_t[,paste0(pp_ln_exact$sid,".01A")]
  pp_ln_exact = pp_ln_exact[,c("n0","n1a","n1b")]
  getWGCNA(expr = expr_ln_exact,SP01grp = pp_ln_exact,TopPerc = 0.2) ->
    fo_ln
  merge(fo_ln,g_info,by.x = "mols",by.y = "gene_id") ->
    fo_ln
  merge(fo_ln,expr_tumor_ln_mean,
        by.x = "mols",by.y = "entrez",all.x = T) ->
    fo_ln_mean
  
  
  library(Seurat)
  library(SeuratObject)
  library(monocle3)
  
  sobj = expr_tumor[,intersect(paste0(c(n0,n1a,n1b),".01A"),
                               colnames(expr_tumor))]
  CreateSeuratObject(counts = sobj,
                     project = "mysc",
                     min.cells = 3,
                     min.features = 200) ->
    s0
  
  cell_metadata = colnames(sobj)
  metadata = merge(data.frame(entrez = expr_rn),
                   g_info,by.x = "entrez",by.y = "gene_id",all.x = T,sort = F)
  
  count <- GetAssayData(s0, assay = "RNA", slot = "counts")
  
}


cls1 = cls_info$SID[cls_info$cls == 1]
cls02 = cls_info$SID[cls_info$cls %in% c(0,2)]
intersect(cls1,p_info$PID[p_info$path_n == 0]) ->
  n0_classic
intersect(cls02,p_info$PID[p_info$path_n %in% c(1,2)]) ->
  n1_classic
p_info_new = p_info[p_info$PID %in% c(n0_classic,n1_classic),]
p_info_new[,c("times","patient.vital_status","path_m","path_n","path_t")] -> df
library(survival)
#pdf(paste0("P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,14),".pdf"),9,6)
sfit <- survfit(Surv(times, patient.vital_status)~path_t, data=df)#print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
sfit <- survfit(Surv(times, patient.vital_status)~path_n, data=df)#print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
sfit <- survfit(Surv(times, patient.vital_status)~path_m, data=df)#print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
#dev.off()



















rm(expr);rm(expr_tumor);rm(expr_tumor_tpm);rm(expr_tumor_tpm_mean);
if("Dist&L1000" != 0){
  library(network)
  library(igraph)
  
  library(cmapR)
  L1000INFO = "~/Downloads/BIGDATA//GSE92742_Broad_LINCS_sig_info.txt"
  col_meta <- read.delim(L1000INFO, sep="\t", stringsAsFactors=F)
  db_lincs_drugs = read.csv("~/Downloads/BIGDATA/db_lincs92742_drugs.csv")$drugs
  
  load("~/Library/CloudStorage/OneDrive-Personal/STRING/MM/ppi500.Rdata");
  pps = distances(ppi[[3]], mode="out")
  ppid = degree(ppi[[3]])
  ggplot() + geom_density(aes(x = log2(V(ppi[[3]])$degree)))+
    xlab("Degree")+
    ylab("")
}


PSET2 = 0.01;LOGFCSET2 = 1
UPs = subset(x,pValue_n0_n1a < PSET2 & log2FoldChange_n0_n1a > LOGFCSET2)$symbol
UPs_hid = subset(g_info,symbol %in% stringr::str_to_upper(UPs))$gene_id
DNs = subset(x,pValue_n0_n1a < PSET2 & log2FoldChange_n0_n1a < -LOGFCSET2)$symbol
DNs_hid = subset(g_info,symbol %in% stringr::str_to_upper(DNs))$gene_id
drug = read.csv("drugs.csv");names(table(col_meta$cell_id))[1:76] -> cl
if("Get Distance" != 0){
  dd = rep("NA",nrow(drug))
  for (tag in 1:nrow(drug)) {
    #tag = 11
    drug[tag,2] %>%
      stringr::str_split(.,",") %>%
      unlist(.) ->
      tg
    getDT2DD(tg = tg,dg = c(UPs,DNs)) -> dd[tag]
  }
  drug$distance = as.numeric(dd)
  drug = drug[,c("Drugs","distance","Target")]
}

load("drug_n1a_n1b");getsave <- function(){
  if(0 != "get L1000"){
    t1 = Sys.time()
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    getL1000(UPs_hid,DNs_hid,cellline = cl) -> l1000
    #l1000 -> l1000a#"11954yes"Time difference of 20.24243 mins
    #l1000 -> l1000b#"10853yes"Time difference of 18.22806 mins
    print(Sys.time() - t1)
  }
  l1000-> drug_n1a_n1b
  save(drug_n1a_n1b,file = "drug_n1a_n1b")
  rm(drug_n1a_n1b)
}
load("drug_n0_n1b");getsave <- function(){
  if(0 != "get L1000"){
    t1 = Sys.time()
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    getL1000(UPs_hid,DNs_hid,cellline = cl) -> l1000
    print(Sys.time() - t1)
  }
  l1000 -> drug_n0_n1b
  save(drug_n0_n1b,file = "drug_n0_n1b")
  rm(drug_n0_n1b)
}
load("drug_n0_n1a");getsave <- function(){
  if(0 != "get L1000"){
    t1 = Sys.time()
    wg = read.csv("~/Downloads/BIGDATA/wellgenes.csv")$id
    getL1000(UPs_hid,DNs_hid,cellline = cl) -> l1000
    print(Sys.time() - t1)
  }
  l1000 -> drug_n0_n1a
  save(drug_n0_n1a,file = "drug_n0_n1a")
  rm(drug_n0_n1a)
}


