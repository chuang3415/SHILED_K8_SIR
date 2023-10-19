#ABC fitting procedure----
#First round initialization----
Fit.group1=function(parms,obs_inc,comm_inc,n,low,high,m,name,prop,v,endT,param,result_m1){
  time=seq(0,endT,by=1)
  day=as.Date("2021-01-01")+ddays(time)
  week=MMWRweek(day)$MMWRweek;week[week==53]=0;week[MMWRweek(day)$MMWRyear==2022]=week[MMWRweek(day)$MMWRyear==2022]+52
  S0=param$S0
  wk_test=param$wk_test;test_prop=param$test_prop
  comm_int=param$comm_int
  comm_age=param$comm_age
  mask=param$mask;schl_test=param$schl_test
  #Initialize
  results=array(NA,dim=c(length(low),n,m),dimnames=list(c("beta_ss1","beta_ss2","beta_se1","beta_se2","beta_ee","omega","mask_unmask","ext_s1","ext_s2","ext_e","alpha","delta","omicron")))
  d=rep(NA,m);w=matrix(NA,nrow=n,ncol=m)
  t=1;i=1;N_iter=0
  w_numer=(prod(high-low))^-1
  out_inc_I=list(NULL);out_inc_VI=list(NULL);I_before=list(NULL);VI_before=list(NULL)
  #First Group initialization
  for(l in 1:length(results[,1,1])){
    results[l,t,]=runif(m,low[l],high[l])
    # results[l,t,]=result_m1[l,]
  }
  
  for(i in 1:m){
    parms$beta_ss1=results[1,t,i];parms$beta_ss2=results[2,t,i]
    parms$beta_se1=results[3,t,i];parms$beta_se2=results[4,t,i]
    parms$beta_ee=results[5,t,i]
    parms$omega=results[6,t,i];parms$mask_unmask=results[7,t,i]
    parms$ext_s1=results[8,t,i];parms$ext_s2=results[9,t,i];parms$ext_e=results[10,t,i]
    parms$alpha=results[11,t,i];parms$delta=results[12,t,i];parms$omicron=results[13,t,i]
    
    v0=v[1,];k=sample(1:nrow(S0),1)
    prop=c(S0[k,2],S0[k,3],mask$simulated_adhr_student[1],mask$average_mask_adhr[1],test_prop$stud_test_prop[1],test_prop$emply_test_prop[1],0.316,sum((49-18)/(100-18)*0.326,(64-50)/(100-18)*0.282,(100-65)/(100-18)*0.115))
    comm_int_0=as.data.frame(comm_age$data[match(S0$school_fip[k],comm_age$county_fips_code)][[1]])
    comm_int_0=comm_int_0[1,6:7]
    pop0=pop_int(prop,v0,S0[k,1],comm_int_0)
    out_1=SIR_model(parms,time[which(week==24)[1]-1],comm_inc,pop0,param)
    out_1=out_1[-1,]
    param2=param;param2$tint=time[which(week==34)[1]]
    comm_int_0=as.data.frame(comm_age$data[match(S0$school_fip[k],comm_age$county_fips_code)][[1]])
    comm_int_0=comm_int_0[which(week==34)[1],6:7]
    prop[4]=mask$average_mask_adhr[which(week==34)[1]];prop[3]=mask$simulated_adhr_student[which(week==34)[1]];prop[7]=0.451;prop[8]=sum((49-18)/(100-18)*0.393,(64-50)/(100-18)*0.284,(100-65)/(100-18)*0.213);prop[5]=test_prop$stud_test_prop[which(week==34)[1]];prop[6]=test_prop$emply_test_prop[which(week==34)[1]]
    pop0=pop_int(prop,v[v$MMWRweek==33,],S0[k,1],comm_int_0)
    out_2=SIR_model(parms,time[which(week==51)[1]-1]-time[which(week==34)[1]],comm_inc,pop0,param2)
    out_2=out_2[-1,]
    
    param2=param;param2$tint=time[which(week==53)[1]]
    comm_int_0=as.data.frame(comm_age$data[match(S0$school_fip[k],comm_age$county_fips_code)][[1]])
    comm_int_0=comm_int_0[which(week==53)[1],6:7]
    prop[4]=mask$average_mask_adhr[which(week==53)[1]];prop[3]=mask$simulated_adhr_student[which(week==53)[1]];prop[7]=0.677;prop[8]=sum((49-18)/(100-18)*0.606,(64-50)/(100-18)*0.454,(100-65)/(100-18)*0.279);prop[5]=test_prop$stud_test_prop[which(week==53)[1]];prop[6]=test_prop$emply_test_prop[which(week==53)[1]]
    pop0=pop_int(prop,v[v$MMWRweek==52,],S0[k,1],comm_int_0)
    out_3=SIR_model(parms,endT-time[which(week==53)[1]],comm_inc,pop0,param2)
    out_3=out_3[-1,]
    out=rbind(out_1,out_2,out_3)
    N=sum(pop0)
    #match the cumulative weekly positive number
    test_match=wk_test[match(obs_inc$MMWRweek,wk_test$MMWRweek),];test_match=test_match[,-1]
    out=group_out(out,endT)
    I=out[[1]][match(obs_inc$MMWRweek,out[[1]]$week),]
    VI=out[[2]][match(obs_inc$MMWRweek,out[[2]]$week),]
    s1=I$student1;s2=I$student2;s3=VI$student1;s4=VI$student2;e1=I$employee;e2=VI$employee
    I_before[[i]]=I;VI_before[[i]]=VI
    
    s1=sapply(1:length(s1),function (x) rbinom(1,s1[x],test_match[x,1])/N)
    s3=sapply(1:length(s3),function (x) rbinom(1,s3[x],test_match[x,1])/N)
    s2=sapply(1:length(s2),function (x) rbinom(1,s2[x],test_match[x,1])/N)
    s4=sapply(1:length(s4),function (x) rbinom(1,s4[x],test_match[x,1])/N)
    e1=sapply(1:length(e1),function (x) rbinom(1,e1[x],test_match[x,2])/N)
    e2=sapply(1:length(e2),function (x) rbinom(1,e2[x],test_match[x,2])/N)
    #match the MMRweek average vaccination rate
    v_match=v[match(obs_inc$MMWRweek,v$MMWRweek),];v_match=v_match[,-1]
    
    #calculate distance of all schools
    student_prop=sum(prop[1:2])
    out_inc_I[[i]]=data.frame("student1"=s1,"student2"=s2,"employee"=e1)
    out_inc_VI[[i]]=data.frame("student1"=s3,"student2"=s4,"employee"=e2)
    
    d[i]=sqrt(sum((test_match[,1]*(obs_inc$student_rate-(s1+s2+s3+s4))*rep(c(4/111,111/111,109/111),c(17,17,12)))^2,(test_match[,2]*(obs_inc$adult_rate-(e1+e2))*rep(c(4/111,111/111,109/111),c(17,17,12)))^2))
    if(i%%1000==0){print(i)}
    N_iter=N_iter+1
  }
  w[1,]=rep(1,m)/m
  d_75=quantile(d,probs=0.75)
  output=list(results=results,w=w,d_75=d_75,d=d,out_I=out_inc_I,out_VI=out_inc_VI,I_before=I_before,VI_before=VI_before)
  saveRDS(output,file=paste(folder,'/ABC_1','_',name,'.Rdata',sep = ""))
  return(output)
}
#Second rounds and beyond----
Fit.groupn=function(output,t,obs_inc,comm_inc,N,low,high,m,name,prop,v,endT,param){
  time=seq(0,endT,by=1)
  day=as.Date("2021-01-01")+ddays(time)
  week=MMWRweek(day)$MMWRweek;week[week==53]=0;week[MMWRweek(day)$MMWRyear==2022]=week[MMWRweek(day)$MMWRyear==2022]+52
  results=output[[1]];w=output[[2]]
  w_numer=(prod(high-low))^-1
  p=length(results[,1,1])
  vr=param$vr
  vacci_e=param$vacci_e
  gamma=param$gamma
  school_days=param$school_days
  S0=param$S0
  wk_test=param$wk_test;test_prop=param$test_prop
  comm_int=param$comm_int
  comm_age=param$comm_age
  mask=param$mask;schl_test=param$schl_test
  out_inc_I=list(NULL);out_inc_VI=list(NULL);I_before=list(NULL);VI_before=list(NULL)
  #setup
  epsilon=output[[3]]
  N_iter=0;d=rep(NA,m)
  draws=results[,t-1,]
  sigmas=rep(NA,p)
  for(j in 1:p){sigmas[j]=var(draws[j,])}
  i=1
  while(i<=m){		
    part=min(which(cumsum(w[t-1,])>runif(1,0,sum(w[t-1,]))))		#choose value to sample
    setup=perturb(draws[,part],sigmas,high,low)
    parms=setup
    v0=v[1,];k=sample(1:nrow(S0),1)
    prop=c(S0[k,2],S0[k,3],mask$simulated_adhr_student[1],mask$average_mask_adhr[1],test_prop$stud_test_prop[1],test_prop$emply_test_prop[1],0.316,sum((49-18)/(100-18)*0.326,(64-50)/(100-18)*0.282,(100-65)/(100-18)*0.115))
    comm_int_0=as.data.frame(comm_age$data[match(S0$school_fip[k],comm_age$county_fips_code)][[1]])
    comm_int_0=comm_int_0[1,6:7]
    pop0=pop_int(prop,v0,S0[k,1],comm_int_0)
    out_1=SIR_model(parms,time[which(week==24)[1]-1],comm_inc,pop0,param)
    out_1=out_1[-1,]
    param2=param;param2$tint=time[which(week==34)[1]]
    comm_int_0=as.data.frame(comm_age$data[match(S0$school_fip[k],comm_age$county_fips_code)][[1]])
    comm_int_0=comm_int_0[which(week==34)[1],6:7]
    prop[4]=mask$average_mask_adhr[which(week==34)[1]];prop[3]=mask$simulated_adhr_student[which(week==34)[1]];prop[7]=0.451;prop[8]=sum((49-18)/(100-18)*0.393,(64-50)/(100-18)*0.284,(100-65)/(100-18)*0.213);prop[5]=test_prop$stud_test_prop[which(week==34)[1]];prop[6]=test_prop$emply_test_prop[which(week==34)[1]]
    pop0=pop_int(prop,v[v$MMWRweek==33,],S0[k,1],comm_int_0)
    out_2=SIR_model(parms,time[which(week==51)[1]-1]-time[which(week==34)[1]],comm_inc,pop0,param2)
    out_2=out_2[-1,]
    
    param2=param;param2$tint=time[which(week==53)[1]]
    comm_int_0=as.data.frame(comm_age$data[match(S0$school_fip[k],comm_age$county_fips_code)][[1]])
    comm_int_0=comm_int_0[time[which(week==53)[1]],6:7]
    prop[4]=mask$average_mask_adhr[which(week==53)[1]];prop[3]=mask$simulated_adhr_student[which(week==53)[1]];prop[7]=0.677;prop[8]=sum((49-18)/(100-18)*0.606,(64-50)/(100-18)*0.454,(100-65)/(100-18)*0.279);prop[5]=test_prop$stud_test_prop[which(week==53)[1]];prop[6]=test_prop$emply_test_prop[which(week==53)[1]]
    pop0=pop_int(prop,v[v$MMWRweek==52,],S0[k,1],comm_int_0)
    out_3=SIR_model(parms,endT-time[which(week==53)[1]],comm_inc,pop0,param2)
    out_3=out_3[-1,]
    out=rbind(out_1,out_2,out_3)
    N=sum(pop0)
    #match the cumulative weekly positive number
    test_match=wk_test[match(obs_inc$MMWRweek,wk_test$MMWRweek),];test_match=test_match[,-1]
    out=group_out(out,endT)
    I=out[[1]][match(obs_inc$MMWRweek,out[[1]]$week),]
    VI=out[[2]][match(obs_inc$MMWRweek,out[[2]]$week),]
    s1=I$student1;s2=I$student2;s3=VI$student1;s4=VI$student2;e1=I$employee;e2=VI$employee
    
    s1=sapply(1:length(s1),function (x) rbinom(1,s1[x],test_match[x,1])/N)
    s3=sapply(1:length(s3),function (x) rbinom(1,s3[x],test_match[x,1])/N)
    s2=sapply(1:length(s2),function (x) rbinom(1,s2[x],test_match[x,1])/N)
    s4=sapply(1:length(s4),function (x) rbinom(1,s4[x],test_match[x,1])/N)
    e1=sapply(1:length(e1),function (x) rbinom(1,e1[x],test_match[x,2])/N)
    e2=sapply(1:length(e2),function (x) rbinom(1,e2[x],test_match[x,2])/N)
    #match the MMRweek average vaccination rate
    v_match=v[match(obs_inc$MMWRweek,v$MMWRweek),];v_match=v_match[,-1]
    #calculate distance of all schools
    student_prop=sum(prop[1:2])
    d[i]=sqrt(sum((test_match[,1]*(obs_inc$student_rate-(s1+s2+s3+s4))*rep(c(4/111,111/111,109/111),c(17,17,12)))^2,(test_match[,2]*(obs_inc$adult_rate-(e1+e2))*rep(c(4/111,111/111,109/111),c(17,17,12)))^2))
    if(d[i]<=epsilon){
      results[,t,i]=setup
      pdftmp=rep(NA,p)
      for(j in 1:p){pdftmp[j]=ifelse(setup[j]>-sigmas[j]||setup[j]<sigmas[j],1/(2*sigmas[j]),0)}
      pdfK=prod(pdftmp)	#backwards filter
      w_denom=sum(w[t-1,part])*pdfK
      w[t,i]=w_numer/w_denom	#set weight
      out_inc_I[[i]]=data.frame("student1"=s1,"student2"=s2,"employee"=e1)
      out_inc_VI[[i]]=data.frame("student1"=s3,"student2"=s4,"employee"=e2)
      I_before[[i]]=I;VI_before[[i]]=VI
      if(i%%1000==0){print(c(t,i,N_iter))}
      i=i+1
    }
    N_iter=N_iter+1
  }
  d_75=quantile(d,probs=0.75)
  w[t,]=w[t,]/sum(w[t,])	#normalizing weights
  output=list(results=results,w=w,d_75=d_75,d=d,out_I=out_inc_I,out_VI=out_inc_VI,I_before=I_before,VI_before=VI_before)
  saveRDS(output,file=paste(folder,'ABC','_',t,'_',name,'.Rdata',sep = ""))
  return(output)
}
#Random perturbation function----
perturb=function (draw,sigmas,high,low) {
  for(j in 1:length(draw)){
    if (j ==3) {high[3]=draw[1]}
    if (j ==4) {high[4]=draw[2]}
    if (j ==11) {high[11]=draw[13]}
    if (j ==12) {high[12]=draw[13]}
    draw[j]=draw[j]+(sigmas[j]*runif(1,min=-1,max=1))	#move up and down
    draw[j]=ifelse(draw[j]<low[j],low[j],draw[j])
    draw[j]=ifelse(draw[j]>high[j],high[j],draw[j])
  }
  return(draw)
}
#Combined function for output----
group_out=function (out,endT) {
  time=seq(0,endT,by=1)
  day=as.Date("2021-01-01")+ddays(time)
  week=MMWRweek(day)$MMWRweek;week[week==53]=0;week[MMWRweek(day)$MMWRyear==2022]=week[MMWRweek(day)$MMWRyear==2022]+52
  out_inc_stud1=out[c(seq(1,nrow(out),by=12),seq(2,nrow(out),by=12),seq(3,nrow(out),by=12),seq(4,nrow(out),by=12)),]
  out_inc_stud2=out[c(seq(5,nrow(out),by=12),seq(6,nrow(out),by=12),seq(7,nrow(out),by=12),seq(8,nrow(out),by=12)),]
  
  out_inc_I_stud1=tapply(out_inc_stud1[,4],out_inc_stud1[,2],sum)
  out_inc_VI_stud1=tapply(out_inc_stud1[,7],out_inc_stud1[,2],sum)
  out_inc_I_stud2=tapply(out_inc_stud2[,4],out_inc_stud2[,2],sum)
  out_inc_VI_stud2=tapply(out_inc_stud2[,7],out_inc_stud2[,2],sum)
  
  out_inc_emply=out[c(seq(9,nrow(out),by=12),seq(10,nrow(out),by=12),seq(11,nrow(out),by=12),seq(12,nrow(out),by=12)),]
  out_inc_I_emply=tapply(out_inc_emply[,4],out_inc_emply[,2],sum)
  out_inc_VI_emply=tapply(out_inc_emply[,7],out_inc_emply[,2],sum)
  
  #group weekly cumulative positive case
  wk_out=week[c(seq(1,which(week==24)[1]-1),seq(which(week==34)[1],which(week==51)[1]-1),seq(which(week==53)[1],endT+1))]
  out_inc_I_stud1=tapply(out_inc_I_stud1,wk_out,sum)
  out_inc_VI_stud1=tapply(out_inc_VI_stud1,wk_out,sum)
  out_inc_I_stud2=tapply(out_inc_I_stud2,wk_out,sum)
  out_inc_VI_stud2=tapply(out_inc_VI_stud2,wk_out,sum)
  out_inc_I_emply=tapply(out_inc_I_emply,wk_out,sum)
  out_inc_VI_emply=tapply(out_inc_VI_emply,wk_out,sum)
  out_inc_I=data.frame("week"=unique(wk_out),"student1"=out_inc_I_stud1,"student2"=out_inc_I_stud2,"employee"=out_inc_I_emply)
  out_inc_VI=data.frame("week"=unique(wk_out),"student1"=out_inc_VI_stud1,"student2"=out_inc_VI_stud2,"employee"=out_inc_VI_emply)
  return(list(out_inc_I,out_inc_VI))
}