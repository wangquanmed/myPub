# General Preset
library(magrittr)
library(reshape2)

ps = p_info[,c("PID","times","patient.vital_status","path_t","path_n","path_m","pstage")]
ps = ps[ps$path_n %in% c(0,1,2) &
          ps$path_t %in% c(1,2,3,4) &
          ps$path_m %in% c("m0","m1"),]

# Fig A
sfit <- survfit(Surv(times, patient.vital_status)~path_t, data=ps)
print(sfit)
#pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,4)
ggsurvplot(sfit, conf.int=F, pval=T)
#graphics.off()

# Fig B
#pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,4)
ggsurvplot(survfit(Surv(times, patient.vital_status)~path_n, data=ps), 
           conf.int=F, pval=T)
#graphics.off()

# Fig C
#pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,4)
ggsurvplot(survfit(Surv(times, patient.vital_status)~path_m, data=ps), 
           conf.int=F, pval=T)
#graphics.off()

# Fig D
sfit <- survfit(Surv(times, patient.vital_status)~pstage, data=ps)
print(sfit)
#pdf(paste0("./img/P",substr(gsub('[^[:alnum:]]', '', Sys.time()),3,13),".pdf"),6,4)
ggsurvplot(survfit(Surv(times, patient.vital_status)~pstage, data=ps), 
           conf.int=F, pval=T)
#graphics.off()

p_info$vital_status = ifelse(p_info$vital_status == "dead",1,0)
p_info_m0 = p_info[p_info$path_m == "m0",]
table(p_info_m0$path_n)
table(p_info_m0$his)
summary(as.numeric(p_info_m0$age))
#write.csv(p_info_m0,file = "./data/p_info_m0.csv",row.names = F)



getCorVar <- function(){
  ppt = p_info_m0[,c("path_n",colnames(p_info)[14:264])]
  library(Hmisc)
  ppt_p = data.frame(rcorr(as.matrix(ppt))$P)
  data.frame(a = rownames(ppt_p),b = (ppt_p$path_n < 0.0001)) -> bb
  bb = subset(bb,b == T)
  ppt = ppt[,c(bb$a)]
  
  library(PerformanceAnalytics)
  
  
  return(ppt)
};

getCorVar() ->
  p_info_cor
p_info_cor %>% colnames()
#my_data <- p_info_cor[,1:9]
#chart.Correlation(my_data, histogram=TRUE, pch=19)

getMv2BvDF <- function(df){
  uid = paste0("uid",rownames(df))
  df_cast = data.frame(uid = uid)
  for (i in 1:ncol(df)) {
    v = paste0(colnames(df)[i],"_")
    df[,i] = paste0(v,df[,i])
    cbind(data.frame(uid = uid),v = df[,i]) -> df2
    reshape2::dcast(df2,formula = as.formula("uid~v"),value.var = "v") -> df2
    N = ncol(df2)
    cbind(df_cast,df2[,2:N]) -> df_cast
  }
  df_cast[not(is.na(df_cast))] <- 1
  df_cast[(is.na(df_cast))] <- 0
  df_cast$uid = uid
  return(df_cast)
}
#getMv2BvDF(df = p_info_cor) -> p_info_cor
#p_info_cor = p_info_cor[,2:ncol(p_info_cor)]

p_info[,c(6,8:11)] ->
  p_info_ghtnm
p_info_s = cbind(p_info_ghtnm,p_info_cor);rm(p_info_ghtnm);rm(p_info_cor)
colnames(p_info_s) = stringr::str_replace_all(colnames(p_info_s), "[^[:alnum:]]", "_")
cbind(p_info[,c(1,2,4,7)],p_info_s) -> 
  p_info;rm(p_info_s)


