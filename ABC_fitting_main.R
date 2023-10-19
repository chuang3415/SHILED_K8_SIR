library("dplyr")
library("tidyverse")
library("stringr")
library(reshape2)
library(ggplot2)
library(lubridate)
library(readxl)
library(MMWRweek)
library(mc2d)

setwd("#set working directory")
folder="#set the file path to store output files"
result_m1=readRDS("#read output distributions from the first ABC fitting procedure")
low=c(0,0,0,0,0,0,0.44,0,0,0,1,1,1)
high=c(1.2,1.2,1,1,1.2,1,1,3,3,3,3.5,3.5,5)
n=10 #generation
m=1000 #particles
result=array(NA,dim=c(dim(result_m1$results)[1],m),dimnames=list(c("beta_ss1","beta_ss2","beta_se1","beta_se2","beta_ee","omega","mask_unmask","ext_s1","ext_s2","ext_e","alpha","delta","omicron")))
result=result_m1$results[,ncol(result_m1$results),sample(c(1:dim(result_m1$results)[3]),m)] #resample previous results as prior from first ABC fitting procedure

#initial vaccination rate
v=cbind(vacci_date_k8 %>% group_by(MMWRyear,MMWRweek) %>% group_keys(),"youth"=do.call(c,vacci_date_k8 %>% group_by(MMWRyear,MMWRweek) %>% group_map(~mean(.$youth))),"adult"=do.call(c,vacci_date_k8 %>% group_by(MMWRyear,MMWRweek) %>% group_map(~mean(.$adult))))
v$MMWRweek[v$MMWRyear==2020]=0;v$MMWRweek[v$MMWRyear==2022]=v$MMWRweek[v$MMWRyear==2022]+52;v=v[,-1]
prop=c(stud_prop_1,stud_prop_2,0.8,0.6,test_prop$stud_test_prop[1],test_prop$emply_test_prop[1],0.316,sum((49-18)/(100-18)*0.326,(64-50)/(100-18)*0.282,(100-65)/(100-18)*0.115)) #Initialize student proportions,mask_student,mask_adult,testing proportions in student and employees, and seroprevalence in students and adults.
name="k8_pop0_2nd_update1" #output file name

fit_param=data.frame(beta_ss1=0.03,beta_ss2=0.03,beta_se1=0.03,beta_se2=0.01,beta_ee=0.1,omega=0.01,mask_unmask=0.5,ext_s1=0.1,ext_s2=0.1,ext_e=1,alpha=1.1,delta=1.2,omicron=1.5)
obs_inc=obs_inc_0 #input
endT=444
#duration to date and month
times=seq(0,endT,by=1)
date=as.Date("2021-01-01")+times
times=MMWRweek(date)
#initialize school days
school_days=rep(1,nrow(times))
school_days[times$MMWRday %in% c(1,7)]=rep(0,sum(times$MMWRday %in% c(1,7)))
holidays=list(date=as.Date(c("2021-01-01","2021-01-04","2021-01-18","2021-02-12","2021-02-15","2021-03-22","2021-03-23","2021-03-24","2021-03-25","2021-03-26","2021-03-29","2021-03-30","2021-03-31","2021-04-01","2021-04-02","2021-05-31","2021-06-11","2021-08-23","2021-09-06","2021-10-11","2021-11-24","2021-11-25","2021-11-26","2021-12-23","2022-01-03","2022-02-18","2022-02-21","2022-03-04")),week=c(seq(24,33,by=1),47,51,52)) 
school_days[times$MMWRweek %in% holidays$week]=rep(0,sum(times$MMWRweek %in% holidays$week))
school_days[date %in% holidays$date]=rep(0,sum(date %in% holidays$date))
#vaccination efficacy and infectious period
mth=month(date);mth[366:length(mth)]=mth[366:length(mth)]+12
vacci_e=cbind(mth,rep(NA,length(mth)),rep(NA,length(mth)))
vacci_e[vacci_e[,1]<12,2]=rtriang(sum(vacci_e[,1]<12),min=0.79,mode=0.92,max=0.97)
vacci_e[vacci_e[,1]>11,2]=rtriang(sum(vacci_e[,1]>11),min=0.51,mode=0.81,max=0.93)
vacci_e[vacci_e[,1]<7,3]=rtriang(sum(vacci_e[,1]<7),min=0.81,mode=0.91,max=0.96)
vacci_e[vacci_e[,1]<9&vacci_e[,1]>6,3]=rtriang(sum(vacci_e[,1]<9&vacci_e[,1]>6),min=0.26,mode=0.66,max=0.84)
vacci_e[vacci_e[,1]>8,3]=rtriang(sum(vacci_e[,1]>8),min=0.48,mode=0.67,max=0.80)
vacci_e[vacci_e[,1]>12,2]=runif(sum(vacci_e[,1]>12),min=0.09,max=0.79)
vacci_e[vacci_e[,1]>12,3]=runif(sum(vacci_e[,1]>12),min=0.18,max=0.62)
gamma=cbind(mth,rep(NA,length(mth)))
gamma[gamma[,1]<7,2]=rtriang(sum(gamma[,1]<7),min=4.7,mode=5.5,max=6.5)
gamma[gamma[,1]>6,2]=rtriang(sum(gamma[,1]>6),min=4.1,mode=4.7,max=5.6)

comm_inc=cdc_22_wt$positive_rate[3:nrow(cdc_22_wt)] #CDC community positive case weighted by county population
vr=apply(vacci_date_k8[,4:5],1,function (x) rep(x,c(8,4))) #vaccination rate increments CDC
wk_test=data.frame(test_prop %>% mutate(MMWRweek(.$date)) %>% group_by(MMWRyear,MMWRweek) %>% group_keys(),"student_prop"=do.call(c,test_prop %>% mutate(MMWRweek(.$date)) %>% group_by(MMWRyear,MMWRweek) %>% group_map(~sum(.$stud_test_prop))),"employee_prop"=do.call(c,test_prop %>% mutate(MMWRweek(.$date)) %>% group_by(MMWRyear,MMWRweek) %>% group_map(~sum(.$emply_test_prop))))
wk_test$student_prop[which(wk_test$student_prop>1)]=1
wk_test$employee_prop[which(wk_test$employee_prop>1)]=1
wk_test$MMWRweek[wk_test$MMWRyear==2022]=wk_test$MMWRweek[wk_test$MMWRyear==2022]+52
wk_test$MMWRweek[wk_test$MMWRyear==2020]=0
wk_test=wk_test[,-1]
schl_test=schl_test %>% filter(MMWRyear!=2022,!(MMWRweek %in% c(seq(24,33,by=1),51,52)))

param=list("vr"=vr,"vacci_e"=vacci_e,"gamma"=gamma,"school_days"=school_days,"S0"=pop_school,"test_prop"=test_prop,"wk_test"=wk_test,"tint"=0,"comm_int"=cdc_in_21,"comm_age"=cdc_21_county_age,"mask"=cc_data_IL,"folder"=folder,"schl_test"=schl_test)
#remove summer break testing data, winter break 51,52
obs_inc=obs_inc %>% filter(!(MMWRweek %in% c(seq(24,33,by=1),51,52)))

#Run first generation
out=list(NULL)
out[[1]]=Fit.group1(fit_param,obs_inc,comm_inc,n,low,high,m,name,prop,v,endT,param,result)
#Run the rest of the generation
tint=2 #start from 2nd generation
for (t in tint:n) {
  out[[t]]=Fit.groupn(out[[t-1]],t,obs_inc,comm_inc,n,low,high,m,name,prop,v,endT,param)
}
output=list(out[[n]]$results,out[[n]]$w,out[[n]]$d,out[[n]]$out_I[[1]],out[[n]]$out_VI[[1]])