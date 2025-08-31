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
#install.packages(paste0(c4c,"U/xbox/"),repos = NULL,type = "source")
library(xbox)
pset = 0.05
logfcset = log2(1.5)
g_info = read.csv("./data/g_info.csv")
p_info = read.csv("./data/p_info_m0.csv")
n0 = p_info$PID[p_info$path_n == "0"]
n1a = p_info$PID[p_info$path_n == "1"]
n1b = p_info$PID[p_info$path_n == "2"]

expr = read.csv("./data/expr.csv")
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

library(clusterProfiler)
t2g = getT2G(org = "hs")
lst1 = getLST(dfexp = x,log2fc = "log2FoldChange_n0_n1a",gene_id = "gene_id")
lst2 = getLST(dfexp = x,log2fc = "log2FoldChange_n0_n1b",gene_id = "gene_id")
lst3 = getLST(dfexp = x,log2fc = "log2FoldChange",gene_id = "gene_id")
gs1 = getGSEA(LST = lst1,org = "hs",t2g = t2g,getRS = F)
gs1$prank_1 = dplyr::percent_rank(gs1$enrichmentScore)
gs2 = getGSEA(LST = lst2,org = "hs",t2g = t2g,getRS = F)
gs2$prank_2 = dplyr::percent_rank(gs2$enrichmentScore)
gs3 = getGSEA(LST = lst3,org = "hs",t2g = t2g,getRS = F)
gs3$prank_3 = dplyr::percent_rank(gs3$enrichmentScore)

merge(merge(x = gs1[,c("Description","prank_1")],
            y = gs2[,c("Description","prank_2")],
            all = T,
            by = "Description"),
      gs3[,c("Description","prank_3","EnrichType")],
      all = T,
      by = "Description") ->
  fo
fo = subset(fo,EnrichType %in% c("GOBP","GOCC","GOMF"))
fo = subset(fo,EnrichType == "KEGG")



library(energy)
n_points <- nrow(fo)
M_simulations <- 1000
generate_uniform_cube_points <- function(n) {matrix(runif(n * 3), ncol = 3)}
set.seed(42)
uniform_cube_obs <- generate_uniform_cube_points(n_points)
non_uniform_cube_obs <- matrix(NA, nrow = n_points, ncol = 3)
for (i in 1:n_points) {
  non_uniform_cube_obs[i, 1] <- fo[i,2]
  non_uniform_cube_obs[i, 2] <- fo[i,3]
  non_uniform_cube_obs[i, 3] <- fo[i,4]
}
# 生成用于与观测数据进行比较的均匀参照样本
uniform_cube_reference_sample <- generate_uniform_cube_points(n_points)
# 合并观测样本和参照样本
combined_uniform_obs <- rbind(uniform_cube_obs, uniform_cube_reference_sample)
combined_non_uniform_obs <- rbind(non_uniform_cube_obs, uniform_cube_reference_sample)
# 计算观测统计量
T_obs_uniform <- eqdist.e(combined_uniform_obs, 
                          sizes = c(nrow(uniform_cube_obs), 
                                    nrow(uniform_cube_reference_sample)))
T_obs_non_uniform <- eqdist.e(combined_non_uniform_obs, 
                              sizes = c(nrow(non_uniform_cube_obs), 
                                        nrow(uniform_cube_reference_sample)))
message(paste("均匀数据观测 Energy 统计量:", T_obs_uniform))
message(paste("非均匀数据观测 Energy 统计量:", T_obs_non_uniform))
# 模拟在原假设（均匀分布）下的统计量分布
simulated_stats_null <- numeric(M_simulations)
for (m in 1:M_simulations) {
  # 在原假设下生成 N 个均匀分布的点 (作为“模拟观测数据”)
  sim_data_null <- generate_uniform_cube_points(n_points)
  # 生成新的均匀参照样本（每次模拟都独立）
  sim_uniform_reference_sample <- generate_uniform_cube_points(n_points)
  # 合并模拟数据和参照样本
  combined_sim_data_null <- rbind(sim_data_null, sim_uniform_reference_sample)
  # 计算模拟统计量
  simulated_stats_null[m] <- eqdist.e(combined_sim_data_null, 
                                      sizes = c(nrow(sim_data_null), 
                                                nrow(sim_uniform_reference_sample)))}
# --- 6. 计算并打印 p 值 ---
# 均匀数据的 p 值
p_value_uniform <- (sum(simulated_stats_null >= T_obs_uniform) + 1) / (M_simulations + 1)
message("\n--- P 值 ---")
message(paste("均匀数据 P 值 (Energy Test 结合蒙特卡洛):", p_value_uniform))
# 非均匀数据的 p 值
p_value_non_uniform <- (sum(simulated_stats_null >= T_obs_non_uniform) + 1) / (M_simulations + 1)
message(paste("非均匀数据 P 值 (Energy Test 结合蒙特卡洛):", p_value_non_uniform))




fo$prank_1  = dplyr::percent_rank(fo$prank_1)
fo$prank_2  = dplyr::percent_rank(fo$prank_2)
fo$prank_3  = dplyr::percent_rank(fo$prank_3)
#fo$Description = stringr::str_to_sentence(substr(fo$Description,6,100))
fo$alpha = abs(apply(fo[,2:4], 1, max) - apply(fo[,2:4], 1, min) - 0.5) + 0.5

library(ggtern)
cate <- dplyr::case_when(
  abs(fo$prank_1 - 0.5) > 0.25 &abs(fo$prank_2 - 0.5) > 0.25 & abs(fo$prank_3 - 0.5) > 0.25 ~ "Effective in All 3 Transition",
  abs(fo$prank_1 - 0.5) > 0.25 &abs(fo$prank_2 - 0.5) > 0.25 ~ "N0 -> N1a & N0 -> N1b",
  abs(fo$prank_2 - 0.5) > 0.25 &abs(fo$prank_3 - 0.5) > 0.25 ~ "N0 -> N1b & N1a -> N1b",
  abs(fo$prank_1 - 0.5) > 0.25 &abs(fo$prank_3 - 0.5) > 0.25 ~ "N0 -> N1a & N1a -> N1b",
  TRUE ~ "Other"
)
table(cate)
fo$cate = cate
pdf(paste0("./data/","Fig4b.pdf"),width = 2.9,height = 2.6)
fo %>%
  ggtern(.,
         aes(x = prank_1,
             y = prank_3,
             z = prank_2,
             color = cate,
             #shape = EnrichType,
             alpha = alpha))+
  stat_density_tern() +
  #geom_confidence_tern()+
  #geom_hex_tern(binwidth=0.05,alpha=0.2)+
  theme_rgbw()+
  #geom_point(size = 2)+
  labs(x = "N0 ->N1a", y = "N1a ->N1b", z = "N0 ->N1b") +
  theme(legend.position = "right")
graphics.off()


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
TOPN_df3 = subset(gs3,)

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

