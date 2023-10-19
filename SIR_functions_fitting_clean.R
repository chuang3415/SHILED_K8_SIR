#Functions used in the SIR fitting procedures
Rates_tau=function(params,x,tstep,t,school_days){
  S=x[,1];I=x[,2];R=x[,3];V=x[,4];VI=x[,5];Q=x[,6]
  rates=matrix(NA,ncol=13,nrow = nrow(x))
  N=sum(x)
  #initialize rates
  #infection unvaccinated
  rates[,1]=params$inf[t]*school_days[t]*S*params$beta%*%(I+VI)/(N-sum(Q))
  rates[,9]=params$inf[t]*params$ext*params$c[t]*S
  rates[,11]=params$v[,t]*S
  #recovery unvaccinated
  rates[,2]=params$gamma[t]*I 
  #infection vaccinated
  rates[,3]=params$inf[t]*school_days[t]*(1-params$e[,t])*V*params$beta%*%(I+VI)/(N-sum(Q))
  rates[,10]=params$inf[t]*(1-params$e[,t])*params$ext*params$c[t]*V
  #recovery vaccinated
  rates[,4]=params$gamma_v[t]*VI
  #loss of immunity
  rates[,5]=params$omega*R
  #testing and quarantine
  rates[,6]=params$t*(VI) 
  rates[,7]=params$t*(I)
  #recovery during quarantine
  rates[,8]=params$gamma_q*Q #recovery during quarantine
  rates=rates*tstep #multiply all rates by time step
  return(rates)
}

Events_tau=function(rates,x,a){
  S=x[1];I=x[2];R=x[3];V=x[4];VI=x[5];Q=x[6]
  #drawing rates from distributions
  inf_ex=min(S,rpois(1,rates[9]))
  inf_ex_v=min(V,rpois(1,rates[10]))
  recovery=min(I,rpois(1,rates[2]))
  loss_immu=min(R,rpois(1,rates[5]))
  rec_v=min(VI,rpois(1,rates[4]))
  rec_q=min(Q,rpois(1,rates[8]))
  
  a_u=a[1];a_v=a[2]
  qua_I=a_u*min(I,rpois(1,rates[7]))
  qua_VI=a_v*min(VI,rpois(1,rates[6]))
  
  #constrains
  infection=min(S,rpois(1,rates[1])+inf_ex) #infection<=S
  inf_v=min(V,rpois(1,rates[3])+inf_ex_v) #infection vaccinated<=V
  exit_I=min(I,recovery+qua_I) #unvaccinated recovery+quarantine <=I
  exit_VI=min(VI,rec_v+qua_VI) #vaccinated recovery+quarantine <=VI
  if (rates[11]<0) {
    vaccinate_I=0 #vaccinate_I=-rpois(1,abs(rates[11])) county vaccination fluctuation
  } else {
    vaccinate_I=min(S,rpois(1,rates[11]))    
  } #vaccination rate from S to V<=S
  #if recovery+quarantine > I or VI, cap quarantine rate
  if(recovery+qua_I>=I)
  {qua_I=I-recovery}
  if(rec_v+qua_VI>=VI)
  {qua_VI=VI-rec_v}
  quarantine=qua_I+qua_VI
  #if vaccination rate and infection rate > S, cap at S
  if (infection+vaccinate_I>=S) 
  {infection=S-vaccinate_I}
  
  #simulation at each time step
  S=S-infection-vaccinate_I; #unvaccinated susceptible
  I=I+infection-exit_I; #unvaccinated infection
  R=R+recovery+rec_v-loss_immu+rec_q #recovery
  V=V+loss_immu+vaccinate_I-inf_v #vaccinated susceptible
  VI=VI+inf_v-exit_VI #vaccinated infection
  Q=Q+quarantine-rec_q #quarantine
  pop_now=c(S,I,R,V,VI,Q)
  pop_now[pop_now<0]=0
  return(pop_now)
}

Iteration_tau=function(params,tstep,endT,pop0,school_days,t0,mask_adhr,test){
  pop_tau=matrix(NA,nrow=(endT/tstep+1)*nrow(pop0),ncol=7,dimnames=list(NULL,c("time","S","I","R","V","VI","Q")))
  pop_now=pop0
  t=t0;step=12
  pop_tau[1:12,]=c(rep(t,12),pop_now)
  while(t<(t0+endT)){
    if (t!=t0) {
      pop_now_temp=pop_now
      pop_temp_1=sapply(colSums(pop_now[c(1:4),]),function (x) min(rpois(1,x*mask_adhr$simulated_adhr_student[t+1]),x))
      pop_temp_2=sapply(colSums(pop_now[c(5:8),]),function (x) min(rpois(1,x*mask_adhr$simulated_adhr_student[t+1]),x))
      pop_temp_3=sapply(colSums(pop_now[c(9:12),]),function (x) min(rpois(1,x*mask_adhr$average_mask_adhr[t+1]),x))
      pop_temp=rbind(rbind(pop_temp_1,colSums(pop_now[c(1:4),])-pop_temp_1),rbind(pop_temp_2,colSums(pop_now[c(5:8),])-pop_temp_2),rbind(pop_temp_3,colSums(pop_now[c(9:12),])-pop_temp_3))
      pop_now_temp[c(1,3),]=t(apply(pop_temp[c(1,2),],1,function (x) sapply(x,function (y) min(rpois(1,y*test[1,t+1]),y))))
      pop_now_temp[c(5,7),]=t(apply(pop_temp[c(3,4),],1,function (x) sapply(x,function (y) min(rpois(1,y*test[5,t+1]),y))))
      pop_now_temp[c(9,11),]=t(apply(pop_temp[c(5,6),],1,function (x) sapply(x,function (y) min(rpois(1,y*test[9,t+1]),y))))
      pop_now_temp[c(2,4),]=pop_temp[c(1,2),]-pop_now_temp[c(1,3),]
      pop_now_temp[c(6,8),]=pop_temp[c(3,4),]-pop_now_temp[c(5,7),]
      pop_now_temp[c(10,12),]=pop_temp[c(5,6),]-pop_now_temp[c(9,11),]
      pop_now=pop_now_temp
    }
    ratedraw=Rates_tau(params,pop_now,tstep,t+1,school_days) #calculate rates
    for(i in 1:12) {
      a_u=0;a_v=0; #Initialize with no testing
      if(i%%2==1){
        a_u=1;a_v=1
      } #indicator if the testing happen
      a_temp=c(a_u,a_v)
      eventdraw=Events_tau(ratedraw[i,],pop_now[i,],a_temp) #make events happen
      step=step+1 #update time and matrix position
      pop_tau[step,]=c(t+tstep,eventdraw)
      pop_now[i,]=eventdraw #update population
    }
    t=t+tstep
  }
  return(pop_tau)
}
#Population generation and categorization with proportion--------
pop_int=function (prop,v0,S0,comm_int) {
  g1=3 #group 1 number: students in elementary schools, middle schools, or employee
  g2=2 #group 2 number: mask vs. unmask
  g3=2 #group 3 number: test vs. no test
  gr_total=prod(g1,g2,g3)
  
  p1=c(prop[1],prop[2],1-sum(prop[1],prop[2])) #proportion of group1 student in elementary and middle schools and employee
  #proportion of group2 mask vs. unmask
  p2=c(prop[3],prop[3],prop[4]) #mask proportion student (elementary and middle schools) and employee
  #proportion of group3 testing vs. no testing in students and employees
  test=rep(c(prop[5],prop[6]),c(4,2)) #testing proportion in studnet vs. employee
  
  #population in compartments
  ST_pop=c(rmultinom(1,S0,p1))
  vint=unlist(v0[c(2,2,3)]) #vaccination in studnet vs. employee
  ST_pop_v=sapply(1:length(ST_pop),function (x) min(rpois(1,ST_pop[x]*vint[x]),ST_pop[x]))
  ST_pop=data.frame("unvaccinated"=ST_pop-ST_pop_v,"vaccinated"=ST_pop_v)
  SIR_s=matrix(NA,nrow = 6,ncol = 2);SIR_e=matrix(NA,nrow = 6,ncol = 1)
  SIR_s[1:6,]=apply(ST_pop[1:2,],1,function (y) sapply(y,function (x) rmultinom(1,x,c(max(1-prop[7]-comm_int$youth_case,0),comm_int$youth_case,prop[7])))) #0.316
  SIR_s[3,]=SIR_s[3,]+SIR_s[6,];SIR_s[6,]=0
  SIR_e[1:6,]=apply(ST_pop[3,],1,function (y) sapply(y,function (x) rmultinom(1,x,c(max(1-prop[8]-comm_int$adult_case,0),comm_int$adult_case,prop[8])))) 
  SIR_e[3,]=SIR_e[3,]+SIR_e[6,];SIR_e[6,]=0
  
  pop0_i=rbind(t(SIR_s),t(SIR_e))
  mask_pop=t(sapply(1:nrow(pop0_i),function (y) sapply(pop0_i[y,],function (x) min(rpois(1,x*p2[y]),x))))
  pop0=do.call(rbind,lapply(1:nrow(mask_pop),function (x) rbind(mask_pop[x,],pop0_i[x,]-mask_pop[x,])))
  test_pop=t(sapply(1:nrow(pop0),function (y) sapply(pop0[y,],function (x) min(rpois(1,x*test[y]),x))))
  pop0=do.call(rbind,lapply(1:nrow(test_pop),function (x) rbind(test_pop[x,],pop0[x,]-test_pop[x,])))
  
  N=sum(pop0)
  return(pop0)
}