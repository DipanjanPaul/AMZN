amznprof<-function(file1,file2) {
  setwd("~/UPS General/Analytics/Amazon")
  
  r_df<-data.frame()
  r_df2<-data.frame()
  
  sh1<-read.csv(file1,sep=";",colClasses="character",strip.white=TRUE,header=F)
  names(sh1)<-names(sh1)<-c("AC_NR","WND_DT","INF_SRC","CUS_CLS_CD","RCH_TYP_CD","FRS_STS_CD",
                          "DAT_SRC_CD","SBB_IR","RET_SHP_CD","MVM_DRC","MIN_CHG","PKG_QY",
                          "SVC_FEA","CTNR_TYP","SVC_TYP","ACQ_MTH","BIL_TER",
                          "ZN_NR","SVC_F_OR","PKG_WT","GRS_RVN","NET_FST","NET_FNA")
  
  sh1$SVC_F_OR<-ifelse((sh1$RCH_TYP_CD=="02" | sh1$SVC_F_OR == ""), sh1$SVC_FEA, sh1$SVC_F_OR)
  sh1$GRS_RVN<-as.numeric(sh1$GRS_RVN)
  sh1$NET_FST<-as.numeric(sh1$NET_FST)
  sh1$NET_FNA<-as.numeric(sh1$NET_FNA)
  sh1$PKG_QY<-as.numeric(sh1$PKG_QY)
  sh1$PKG_WT<-as.numeric(sh1$PKG_WT)

  l<-list(sh1$AC_NR,sh1$INF_SRC,sh1$CUS_CLS_CD,sh1$FRS_STS_CD,sh1$DAT_SRC_CD,sh1$SBB_IR,sh1$RET_SHP_CD,
          sh1$MVM_DRC,sh1$MIN_CHG,sh1$SVC_F_OR,sh1$CTNR_TYP,sh1$SVC_TYP,sh1$ACQ_MTH,
          sh1$BIL_TER,sh1$ZN_NR)

  pitem<-c("AC_NR","INF_SRC","CUS_CLS_CD","FRS_STS_CD","DAT_SRC_CD","SBB_IR",
           "RET_SHP_CD","MVM_DRC","MIN_CHG","SVC_F_OR","CTNR_TYP","SVC_TYP","ACQ_MTH",
           "BIL_TER","ZN_NR")

  p_len<-length(pitem)
  
  pdata<-split(sh1,l,drop=T)
  
  rm(sh1)
  
  for (i in 1:length(pdata)) {

    if (nrow(r_df) == 0) {
      r_df<-(pdata[[i]][1,pitem][1:p_len])
    } else {
      r_df[i,pitem]<-(pdata[[i]][1,pitem][1:p_len])
    }

    f_wt_all<-factor(pdata[[i]]$PKG_WT)    

    r_df[i,p_len + 1]<-mean(tapply(((pdata[[i]]$GRS_RVN - pdata[[i]]$NET_FST)/pdata[[i]]$GRS_RVN),f_wt_all,mean))
    r_df[i,p_len + 2]<-sd(tapply(((pdata[[i]]$GRS_RVN - pdata[[i]]$NET_FST)/pdata[[i]]$GRS_RVN),f_wt_all,mean))

    if (is.na(r_df[i,p_len + 2])) {
      r_df[i,p_len + 2] <- 0.04 * r_df[i,p_len + 1]
    }
    
    r_df[i,p_len + 3]<-mean(tapply(((pdata[[i]]$GRS_RVN - pdata[[i]]$NET_FNA)/pdata[[i]]$GRS_RVN),f_wt_all,mean))
    r_df[i,p_len + 4]<-sd(tapply(((pdata[[i]]$GRS_RVN - pdata[[i]]$NET_FNA)/pdata[[i]]$GRS_RVN),f_wt_all,mean))

    if (is.na(r_df[i,p_len + 4])) {
      r_df[i,p_len + 4] <- 0.04 * r_df[i,p_len + 3]
    }
                           
    pwkly<-split(pdata[[i]],pdata[[i]]$WND_DT)

    r_dfw<-array(list(), dim=c(8,length(pwkly)))
    
    for (j in 1:length(pwkly)) {
      f_wt<-factor(pwkly[[j]]$PKG_WT)
      
      r_dfw[[1,j]]<-tapply(pwkly[[j]]$PKG_QY,f_wt,sum)
      r_dfw[[2,j]]<-tapply(pwkly[[j]]$GRS_RVN,f_wt,sum)
      r_dfw[[3,j]]<-tapply(pwkly[[j]]$NET_FST,f_wt,sum)
      r_dfw[[4,j]]<-tapply(pwkly[[j]]$NET_FNA,f_wt,sum)
      
      r_dfw[[5,j]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FST)/pwkly[[j]]$GRS_RVN),f_wt,mean)
      r_dfw[[6,j]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FST)/pwkly[[j]]$GRS_RVN),f_wt,sd)
      
      bad<-is.na(r_dfw[[6,j]][])
      r_dfw[[6,j]][bad]<-ifelse(sum(!bad) == 0, 0, mean(r_dfw[[6,j]][!bad]))
                           
      r_dfw[[7,j]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FNA)/pwkly[[j]]$GRS_RVN),f_wt,mean)
      r_dfw[[8,j]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FNA)/pwkly[[j]]$GRS_RVN),f_wt,sd)
    
      bad<-is.na(r_dfw[[8,j]][])
      r_dfw[[8,j]][bad]<-ifelse(sum(!bad) == 0, 0, mean(r_dfw[[8,j]][!bad]))
    }

    # The following will give the mean and SD for each profile for the weekly means   
    
    r_df[i,p_len + 5]<-mean(sapply(r_dfw[5,],mean))
    r_df[i,p_len + 6]<-sd(sapply(r_dfw[5,],mean))
    
    if (is.na(r_df[i,p_len + 6])) {
      r_df[i,p_len + 6] <- 0.04 * r_df[i,p_len + 5]
    }
    
    r_df[i,p_len + 7]<-mean(sapply(r_dfw[7,],mean))
    r_df[i,p_len + 8]<-sd(sapply(r_dfw[5,],mean))
    
    if (is.na(r_df[i,p_len + 8])) {
      r_df[i,p_len + 8] <- 0.04 * r_df[i,p_len + 7]
    }
    
    fac<-pdata[[i]]$WND_DT
    
    r_df[i,p_len + 9]<-mean(tapply(pdata[[i]]$GRS_RVN,fac,sum))
    r_df[i,p_len + 10]<-sd(tapply(pdata[[i]]$GRS_RVN,fac,sum))
    
    if (is.na(r_df[i,p_len + 10])) {
      r_df[i,p_len + 10] <- 0.04 * r_df[i,p_len + 9]
    }
    
    r_df[i,p_len + 11]<-mean(tapply(pdata[[i]]$NET_FST,fac,sum))
    r_df[i,p_len + 12]<-sd(tapply(pdata[[i]]$NET_FST,fac,sum))
    
    if (is.na(r_df[i,p_len + 12])) {
      r_df[i,p_len + 12] <- 0.04 * r_df[i,p_len + 11]
    }
    
    r_df[i,p_len + 13]<-mean(tapply(pdata[[i]]$NET_FNA,fac,sum))
    r_df[i,p_len + 14]<-sd(tapply(pdata[[i]]$NET_FNA,fac,sum))
    
    if (is.na(r_df[i,p_len + 14])) {
      r_df[i,p_len + 14] <- 0.04 * r_df[i,p_len + 13]
    }
    
    r_df[i,p_len + 15]<-mean(tapply(pdata[[i]]$PKG_QY,fac,sum))
    r_df[i,p_len + 16]<-sd(tapply(pdata[[i]]$PKG_QY,fac,sum))
    
    if (is.na(r_df[i,p_len + 16])) {
      r_df[i,p_len + 16] <- 0.04 * r_df[i,p_len + 15]
    }

    if (length(pwkly) > 1) {
        r_df[i,p_len + 17]<-mean(as.numeric(combn((r_dfw[1,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 18]<-mean(as.numeric(combn((r_dfw[2,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 19]<-mean(as.numeric(combn((r_dfw[3,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 20]<-mean(as.numeric(combn((r_dfw[4,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 21]<-mean(as.numeric(combn((r_dfw[5,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 22]<-mean(as.numeric(combn((r_dfw[6,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 23]<-mean(as.numeric(combn((r_dfw[7,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
        r_df[i,p_len + 24]<-mean(as.numeric(combn((r_dfw[8,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))))
    } else {
      r_df[i,p_len + 17]<-c('NA')
      r_df[i,p_len + 18]<-c('NA')
      r_df[i,p_len + 19]<-c('NA')
      r_df[i,p_len + 20]<-c('NA')
      r_df[i,p_len + 21]<-c('NA')
      r_df[i,p_len + 22]<-c('NA')
      r_df[i,p_len + 23]<-c('NA')
      r_df[i,p_len + 24]<-c('NA')
    }

    r_df[i,p_len + 25]<-length(pwkly)
    r_df[i,p_len + 26]<-i
    r_df[i,p_len + 27]<-(pdata[[i]][1,5])


  }
    
  names(r_df)<-c("AC_NR","INF_SRC","CUS_CLS_CD","FRS_STS_CD","DAT_SRC_CD","SBB_IR","RET_SHP_CD",
                 "MVM_DRC","MIN_CHG","SVC_F_OR","CTNR_TYP","SVC_TYP","ACQ_MTH","BIL_TER","ZN_NR",
                 "F_M_INC_FST","F_SD_INC_FST","F_M_INC_FNA","F_SD_INC_FNA","W_M_INC_FST",
                 "W_SD_INC_FST","W_M_INC_FNA","W_SD_INC_FNA","W_M_GRS",
                 "W_SD_GRS","W_M_FST","W_SD_FST","W_M_FNA","W_SD_FNA","W_M_VOL","W_SD_VOL",
                 "KS_M_W_V_WT","KS_M_W_GRS_WT","KS_M_W_FST_WT","KS_M_W_FNA_WT","KS_M_W_I1_WT",
                 "KS_SD_W_I1_WT","KS_M_W_I2_WT","KS_SD_W_I2_WT","WK_CNT","SPL_NR","RCH_TYP_CD")   

  
  sh2<-read.csv(file2,sep=";",colClasses="character",strip.white=TRUE,header=T)
  names(sh2)<-c("AC_NR","WND_DT","INF_SRC","CUS_CLS_CD","RCH_TYP_CD","FRS_STS_CD",
                "DAT_SRC_CD","SBB_IR","RET_SHP_CD","MVM_DRC","MIN_CHG","PKG_QY",
                "SVC_FEA","CTNR_TYP","SVC_TYP","ACQ_MTH","BIL_TER",
                "ZN_NR","SVC_F_OR","PKG_WT","GRS_RVN","NET_FST","NET_FNA")
  
  sh2$SVC_F_OR<-ifelse((sh2$RCH_TYP_CD=="02" | sh2$SVC_F_OR == ""), sh2$SVC_FEA, sh2$SVC_F_OR)
  sh2$GRS_RVN<-as.numeric(sh2$GRS_RVN)
  sh2$NET_FST<-as.numeric(sh2$NET_FST)
  sh2$NET_FNA<-as.numeric(sh2$NET_FNA)
  sh2$PKG_QY<-as.numeric(sh2$PKG_QY)
  sh2$PKG_WT<-as.numeric(sh2$PKG_WT)

  l2<-list(sh2$AC_NR,sh2$INF_SRC,sh2$CUS_CLS_CD,sh2$FRS_STS_CD,sh2$DAT_SRC_CD,sh2$SBB_IR,sh2$RET_SHP_CD,
        sh2$MVM_DRC,sh2$MIN_CHG,sh2$SVC_F_OR,sh2$CTNR_TYP,sh2$SVC_TYP,sh2$ACQ_MTH,
        sh2$BIL_TER,sh2$ZN_NR)

  pdata2<-split(sh2,l2,drop=T)
  rm(sh2)

  for (i in 1:length(pdata2)) {
    
    if (nrow(r_df2) == 0) {
      r_df2<-(pdata2[[i]][1,pitem][1:p_len])
    } else {
      r_df2[i,pitem]<-(pdata2[[i]][1,pitem][1:p_len])
    }
    
    prof_his<-r_df[(r_df$AC_NR == r_df2[i,1] & r_df$INF_SRC == r_df2[i,2] & r_df$CUS_CLS_CD == r_df2[i,3] &                    
                      r_df$FRS_STS_CD == r_df2[i,4] & r_df$DAT_SRC_CD == r_df2[i,5] & r_df$SBB_IR == r_df2[i,6] &
                      r_df$RET_SHP_CD == r_df2[i,7] & r_df$MVM_DRC == r_df2[i,8] & r_df$MIN_CHG == r_df2[i,9] &
                      r_df$SVC_F_OR == r_df2[i,10] & r_df$CTNR_TYP == r_df2[i,11] & r_df$SVC_TYP == r_df2[i,12] & 
                      r_df$ACQ_MTH == r_df2[i,13] & r_df$BIL_TER == r_df2[i,14] & r_df$ZN_NR == r_df2[i,15]),][,]
         
    if (nrow(prof_his) == 1) {

        prof_data<-pdata[[prof_his$SPL_NR]] 
      
        pwkly<-split(prof_data,prof_data$WND_DT)
        
        r_dft<-array(list(), dim=c(8,length(pwkly) + 1))
        
        f_wt2<-factor(pdata2[[i]]$PKG_WT)

        r_df2[i,p_len + 1]<-mean(tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FST)/pdata2[[i]]$GRS_RVN),f_wt2,mean))
        w_m_norm<-(tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FST)/pdata2[[i]]$GRS_RVN),f_wt2,mean))
        
        r_df2[i,p_len + 2]<-ifelse(nrow(pdata2[[i]]) == 1 , 0, sd(tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FST)/pdata2[[i]]$GRS_RVN),f_wt2,mean)))
        w_m_norm2<-(tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FNA)/pdata2[[i]]$GRS_RVN),f_wt2,mean))
        
        r_df2[i,p_len + 3]<-mean(tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FNA)/pdata2[[i]]$GRS_RVN),f_wt2,mean))
        r_df2[i,p_len + 4]<-ifelse(nrow(pdata2[[i]]) == 1 , 0, sd(tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FNA)/pdata2[[i]]$GRS_RVN),f_wt2,mean)))
        
        
        r_dft[[1,1]]<-tapply(pdata2[[i]]$PKG_QY,f_wt2,sum)
        r_dft[[2,1]]<-tapply(pdata2[[i]]$GRS_RVN,f_wt2,sum)
        r_dft[[3,1]]<-tapply(pdata2[[i]]$NET_FST,f_wt2,sum)
        r_dft[[4,1]]<-tapply(pdata2[[i]]$NET_FNA,f_wt2,sum)
        
        r_dft[[5,1]]<-tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FST)/pdata2[[i]]$GRS_RVN),f_wt2,mean)
        r_dft[[6,1]]<-tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FST)/pdata2[[i]]$GRS_RVN),f_wt2,sd)
        
        bad<-is.na(r_dft[[6,1]][])
        r_dft[[6,1]][bad]<-ifelse(sum(!bad) == 0, 0, mean(r_dft[[6,1]][!bad]))
        
        r_dft[[7,1]]<-tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FNA)/pdata2[[i]]$GRS_RVN),f_wt2,mean)
        r_dft[[8,1]]<-tapply(((pdata2[[i]]$GRS_RVN - pdata2[[i]]$NET_FNA)/pdata2[[i]]$GRS_RVN),f_wt2,sd)
        
        bad<-is.na(r_dft[[8,1]][])
        r_dft[[8,1]][bad]<-ifelse(sum(!bad) == 0, 0, mean(r_dft[[8,1]][!bad]))
      
        for (j in 1:length(pwkly)) {
          f_wt<-factor(pwkly[[j]]$PKG_WT)
                
          r_dft[[1,j+1]]<-tapply(pwkly[[j]]$PKG_QY,f_wt,sum)
          r_dft[[2,j+1]]<-tapply(pwkly[[j]]$GRS_RVN,f_wt,sum)
          r_dft[[3,j+1]]<-tapply(pwkly[[j]]$NET_FST,f_wt,sum)
          r_dft[[4,j+1]]<-tapply(pwkly[[j]]$NET_FNA,f_wt,sum)
          
          r_dft[[5,j+1]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FST)/pwkly[[j]]$GRS_RVN),f_wt,mean)
          r_dft[[6,j+1]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FST)/pwkly[[j]]$GRS_RVN),f_wt,sd)
  
          bad<-is.na(r_dft[[6,j+1]][])
          r_dft[[6,j+1]][bad]<-ifelse(sum(!bad) == 0, 0, mean(r_dft[[6,j+1]][!bad]))
          
          r_dft[[7,j+1]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FNA)/pwkly[[j]]$GRS_RVN),f_wt,mean)
          r_dft[[8,j+1]]<-tapply(((pwkly[[j]]$GRS_RVN - pwkly[[j]]$NET_FNA)/pwkly[[j]]$GRS_RVN),f_wt,sd)
          
          bad<-is.na(r_dft[[8,j+1]][])
          r_dft[[8,j+1]][bad]<-ifelse(sum(!bad) == 0, 0, mean(r_dft[[6,j+1]][!bad]))
        }
        

        r_df2[i,p_len + 5]<-sum(pdata2[[i]]$GRS_RVN)
        r_df2[i,p_len + 6]<-sum(pdata2[[i]]$NET_FST)
        r_df2[i,p_len + 7]<-sum(pdata2[[i]]$NET_FNA)
        r_df2[i,p_len + 8]<-sum(pdata2[[i]]$PKG_QY)
        r_df2[i,p_len + 9]<-sum(pdata2[[i]]$PKG_WT)
        
        r_df2[i,p_len + 10]<-mean(as.numeric(combn((r_dft[1,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 11]<-mean(as.numeric(combn((r_dft[2,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 12]<-mean(as.numeric(combn((r_dft[3,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 13]<-mean(as.numeric(combn((r_dft[4,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 14]<-mean(as.numeric(combn((r_dft[5,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 15]<-mean(as.numeric(combn((r_dft[6,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 16]<-mean(as.numeric(combn((r_dft[7,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
        r_df2[i,p_len + 17]<-mean(as.numeric(combn((r_dft[8,]),2,function(x) ifelse(((sum(is.na(x[[1]])) > 0) | (sum(is.na(x[[2]]) > 0))), 0, ks.test(x[[1]],x[[2]])[1]))[1:j]))
 
#        r_df2[i,p_len + 18] <- (table((abs(w_m_norm - prof_his$W_M_INC_FST)) > 2.8 * prof_his$F_SD_INC_FST)[2])/length(w_m_norm)
        r_df2[i,p_len + 18] <- round((abs(mean(w_m_norm) - prof_his$W_M_INC_FST))/(prof_his$F_SD_INC_FST/(sqrt(length(w_m_norm)))),4)
        addchk <- round((abs(mean(w_m_norm2) - prof_his$W_M_INC_FNA))/(prof_his$F_SD_INC_FNA/(sqrt(length(w_m_norm2)))),4)
        
        r_df2[i,p_len + 18] <- ifelse(is.na(r_df2[i,p_len + 18]), 0, r_df2[i,p_len + 18])
        addchk <- ifelse(is.na(addchk), 0, addchk)
        
        if (((round(abs(r_df2[i,p_len + 1] - prof_his$W_M_INC_FST),4) > 3.5 * prof_his$W_SD_INC_FST) |
            (round(r_df2[i,p_len + 2],4) > 1.75 * prof_his$F_SD_INC_FST) | 
            (round(abs(r_df2[i,p_len + 3] - prof_his$W_M_INC_FNA),4) > 3.5 * prof_his$W_SD_INC_FNA) |
            (round(r_df2[i,p_len + 4],4) > 1.75 * prof_his$F_SD_INC_FNA) |
#            (as.numeric(prof_his$KS_M_W_I1_WT) <= 0.1 & (as.numeric(r_df2[i,p_len + 14]) > 1.5 * (as.numeric(prof_his$KS_M_W_I1_WT)))) |
#            (as.numeric(prof_his$KS_M_W_I2_WT) <= 0.1 & (as.numeric(r_df2[i,p_len + 16]) > 1.5 * (as.numeric(prof_his$KS_M_W_I2_WT)))) |
            ((round(prof_his$F_SD_INC_FNA,4) > 0) & (addchk > 2.2)) |
            ((round(prof_his$F_SD_INC_FST,4) > 0) & (r_df2[i,p_len + 18] > 2.2))) 
#            (r_df2[i,p_len + 15] > prof_his$KS_SD_W_I1_WT) |
#            (r_df2[i,p_len + 17] > prof_his$KS_SD_W_I2_WT))
          & (prof_his$WK_CNT > 1)
          & (prof_his$W_M_VOL > 100)) {
          r_df2[i,p_len + 19] <- "ERROR"
        } else {
          if (prof_his$W_M_VOL > 100) {
            r_df2[i,p_len + 19] <- "OK"
          } else {
            r_df2[i,p_len + 19] <- "LH"
          }
        }
        
        r_df2[i,p_len + 20] <- prof_his$SPL_NR
        
        
    } else {
        r_df2[i,p_len + 10]<-c('NA')
        r_df2[i,p_len + 11]<-c('NA')
        r_df2[i,p_len + 12]<-c('NA')
        r_df2[i,p_len + 13]<-c('NA')
        r_df2[i,p_len + 14]<-c('NA')
        r_df2[i,p_len + 15]<-c('NA')
        r_df2[i,p_len + 16]<-c('NA')
        r_df2[i,p_len + 17]<-c('NA')
        r_df2[i,p_len + 18]<-c('NA')
        r_df2[i,p_len + 19] <- "NOHIS"
        r_df2[i,p_len + 19]<-c('NA')
        r_df2[i,p_len + 20]<-c('NA')
        
    }
     
    r_df2[i,p_len + 21]<-(pdata2[[i]][1,5])
    rm(prof_his)
   
  }

  names(r_df2)<-c("AC_NR","INF_SRC","CUS_CLS_CD","FRS_STS_CD","DAT_SRC_CD","SBB_IR","RET_SHP_CD",
                  "MVM_DRC","MIN_CHG","SVC_F_OR","CTNR_TYP","SVC_TYP","ACQ_MTH","BIL_TER","ZN_NR",
                  "F_M_INC_FST","F_SD_INC_FST","F_M_INC_FNA","F_SD_INC_FNA","T_GRS","T_FST",
                  "T_FNA","T_QTY","T_WT","KS_M_W_V_WT","KS_M_W_GRS_WT","KS_M_W_FST_WT",
                  "KS_M_W_FNA_WT","KS_M_W_I1_WT","KS_SD_W_I1_WT","KS_M_W_I2_WT",
                  "KS_SD_W_I2_WT","ERR_FAC","STATUS","SPL_NR","RCH_TYP_CD")  

  write.table(r_df,"amznprof.txt",sep=";",col.names=T,row.names=F)
  write.table(r_df2,"amzncurr.txt",sep=";",col.names=T,row.names=F)
  
  r_list<-list(r_df,r_df2)
  return(r_list)
}
