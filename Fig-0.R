# General Preset
rm(list = ls())
library(magrittr)
library(reshape2)
library(RTCGA)
library(RTCGA.rnaseq)
library(survminer)
library(survival)
library(org.Hs.eg.db)
c4c="~/Library/CloudStorage/OneDrive-Personal/code4code/"
library(xbox)#install.packages(paste0(c4c,"U/xbox/"),repos = NULL,type = "source")

# Specific Preset
t = Sys.time()
pset = 0.05
logfcset = log2(1.5)
scoreSET = 600
merge(unique(toTable(org.Hs.egSYMBOL)),
      unique(toTable(org.Hs.egENSEMBL)),
      by = 'gene_id') %>% 
  merge(.,
        toTable(org.Hs.egCHRLOC),
        by = "gene_id",
        all.x = T) ->
  g_info
g_info = g_info[!duplicated(g_info$gene_id),]
  


getCD <- function(opt){
  RTCGA.clinical::THCA.clinical %>%
    data.frame(.) %>%
    survivalTCGA(.) %>%
    dplyr::mutate(PID = stringr::str_replace_all(bcr_patient_barcode,"-",".")) ->
    p_info
  cd_all = data.frame(RTCGA.clinical::THCA.clinical)
  fr=c("patient.bcr_patient_barcode",
       "patient.gender",
       "patient.age_at_initial_pathologic_diagnosis",
       "patient.biospecimen_cqcf.histological_type",
       "patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
       "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n",
       "patient.stage_event.tnm_categories.pathologic_categories.pathologic_m",
       "patient.stage_event.pathologic_stage",
       "patient.vital_status")
  to=c("patient.bcr_patient_barcode",
       "gender",
       "age",
       "his",
       "path_t",
       "path_n",
       "path_m",
       "pstage",
       "vital_status")
  cd = fr %>% cd_all[,.] # cd means clinical data
  colnames(cd) = to
  cd %<>%  
    dplyr::mutate(
      PID = stringr::str_to_upper(stringr::str_replace_all(patient.bcr_patient_barcode,
                                                           "-",
                                                           ".")))
  if(opt == "all"){#getVarAndTrans
    #remove opt0 Vars except patinent barcode
    cd_all = cd_all[,!(colnames(cd_all) %in% fr[2:length(fr)])] # 1732-->1724
    #get cd with na < 0.1. 
    cd_all %<>% .[,colSums(is.na(.)) < ncol(.)*0.1] #get 1724-->350
    
    #str to factor c(1,2,--)
    cd_all %>% 
      apply(., 2, function(x) as.numeric((as.character(x)))) %>% 
      .[,colSums(is.na(.)) == nrow(.)] %>% 
      colnames(.) -> 
      cd_factor;cd_all[,cd_factor] -> 
      cd_factor #189
    b = apply(cd_factor, 2, function(x) as.numeric(as.factor(x)))
    #remove const and all differents
    b = b[,apply(b, 2, function(x) !max(x) %in% c(1))] # 189 --> 146
    b = b[,apply(b, 2, function(x) !max(x) %in% c(nrow(b)))] # 146 --> 110
    
    cd_all %>% 
      apply(., 2, function(x) as.numeric((as.character(x)))) %>% 
      .[,colSums(is.na(.)) != nrow(.)] %>% 
      colnames(.) -> 
      cd_contin;cd_all[,cd_contin] -> 
      cd_contin # 161
    cd_contin = apply(cd_contin, 2, function(x) as.numeric(x))
    #remove const
    cd_contin[,
              colnames(cd_contin)[apply(cd_contin, 
                                        2, 
                                        function(x) sd(x,na.rm = T)) != 0]] ->
      cd_contin # 141
    
    cd_factorPLUSconti = cbind(data.frame(b),data.frame(cd_contin))
    
    cd = cbind(cd,cd_factorPLUSconti)
    
  }
  merge(p_info,
        cd,
        by = "PID") ->
    p_info
  return(p_info)
}
p_info = getCD(opt = "all")
p_info = p_info[p_info$times >=0,]
subset(p_info,
       path_t %in% c("t1","t1a","t1b","t2","t3","t4","t4a") &
         path_n %in% c("n0","n1a","n1b") &
         p_info$path_m %in% c("m0","m1")) ->
  p_info
ifelse(p_info$path_t %in% c("t1","t1a","t1b"),1,
       ifelse(p_info$path_t %in% c("t2"),2,
              ifelse(p_info$path_t == "t3",3,4))) ->
  p_info$path_t
table(p_info$path_t)
ifelse(p_info$path_n == "n0",0,
       ifelse(p_info$path_n == "n1a",1,
              ifelse(p_info$path_n == "n1b",2,NA))) ->
  p_info$path_n
table(p_info$path_n)

table(p_info$path_n[p_info$path_m == "m0"],p_info$path_t[p_info$path_m == "m0"])
table(p_info$path_n[p_info$path_m == "m1"],p_info$path_t[p_info$path_m == "m1"])

#write.csv(g_info,file = "./data/g_info.csv",row.names = F)
#write.csv(p_info,file = "./data/p_info.csv",row.names = F)
g_info = read.csv("./data/g_info.csv")
p_info = read.csv("./data/p_info_m0.csv")
n0 = p_info$PID[p_info$path_n == "0"]
n1a = p_info$PID[p_info$path_n == "1"]
n1b = p_info$PID[p_info$path_n == "2"]

expr = read.csv("./data/expr.csv")
Sys.time() -t
