## RATON design simulation
## No early stops for first 5 cohorts

## RATON is the old name of SSS design

install.packages(R2jags)
library(rjags)



raton_m1=function(s,prevelence,t,p10,p11,p00,p01,pil,pis,pif,sm,epi1,epi2,epi3,lambda2,mis){
  
  
  decision=NULL ## final decision
  D_n_p = NULL ## neg marker is promising
  D_p_p = NULL ## pos marker is promising
  ss_neg = NULL ## final sample size for neg marker in treatment
  ss_pos = NULL ## final sample size for pos marker in treatment
  ss_neg_1 = NULL ## final sample size for neg marker in control
  ss_pos_1 = NULL ## final sample size for pos marker in control
  
  ss_negall = NULL ## final sample size for neg marker
  ss_posall = NULL ## final sample size for pos marker
  ss_trtall = NULL ## final sample size for trt
  ss_ctrall = NULL ## final sample size for control
  ss_all = NULL ## Total final sample size 
  
  p21_postr = NULL ## response rate for missing marker in treatment
  p20_postr = NULL ## response rate for missing marker in control
  p11_postr = NULL ## response rate for pos marker in treatment
  p10_postr = NULL ## response rate for pos marker in control
  p01_postr = NULL ## response rate for neg marker in treatment
  p00_postr = NULL ## response rate for neg marker in control
  
  p21_post = NULL ## posterior probability for missing marker in treatment
  p20_post = NULL ## posterior probability for missing marker in control
  p11_post = NULL ## posterior probability for pos marker in treatment
  p10_post = NULL ## posterior probability for pos marker in control
  p01_post = NULL ## posterior probability for neg marker in treatment
  p00_post = NULL ## posterior probability for neg marker in control
  
  rate21l=NULL ## randomization rate for missing sample in trt arm
  rate11l=NULL ## randomization rate for pos sample in trt arm
  rate01l=NULL ## randomization rate for neg sample in trt arm
  
  p20=p00*(1-prevelence)+p10*prevelence ## response rate for missing in ctrl
  p21=p01*(1-prevelence)+p11*prevelence ## response rate for missing in trt
  
  RatonMCMC<-"model{


r11~dbinom(p11,n11)
r10~dbinom(p10,n10)
r00~dbinom(p00,n00)
r01~dbinom(p01,n01)



p10=phi1*alpha+(1-phi1)*p00
p11=phi2*beta+(1-phi2)*p01



p00~dbeta(0.5,0.5)
p01~dbeta(0.5,0.5)
alpha~dbeta(0.5,0.5)
beta~dunif(p01,1)
phi1~dbern(pbeta0)
phi2~dbern(pbeta1)

       r21~dbinom(p21,n21)
       p21 <-p11*phi+p01*(1-phi)
       r20~dbinom(p20,n20)
       p20 <-p10*phi+p00*(1-phi)
}"
  
  for (l in 1: sm){
    rate11=0.5 ## 1:1 randomization at the beginning
    rate01=0.5 ## 1:1 randomization at the beginning
    rate21=0.5 ## 1:1 randomization at the beginning
    s2=rbinom(1,s,mis)  ## sample size for marker missing samples
    s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
    s0=s-s2-s1 ## sample size for neg samples
    
    n21=rbinom(1,s2,rate21) ## marker missing samples in treatment
    n11=rbinom(1,s1,rate11) ## pos samples in treatment
    n01=rbinom(1,s0,rate01) ## neg samples in treatment
    n20=s2-n21## marker missing samples in control
    n10=s1-n11## pos samples in control
    n00=s0-n01 ## neg samples in control
    
    r20=rbinom(1,n20,p20) ## response missing samples  in ctrl
    r21=rbinom(1,n21,p21) ## response missing samples  in trt
    
    r10=rbinom(1,n10,p10) ## response samples (+) in ctrl
    r11=rbinom(1,n11,p11) ## response samples (+) in trt
    
    r00=rbinom(1,n00,p00) ## response samples (-) in ctrl
    r01=rbinom(1,n01,p01) ## response samples (-) in trt
    
    earlystop=0
    
    
    for(i in 1: 4){
      
      a21=r21+0.5
      b21=0.5+n21-r21
      a20=r20+0.5
      b20=0.5+n20-r20
      
      a11=r11+0.5
      b11=0.5+n11-r11
      a01=r01+0.5
      b01=0.5+n01-r01
      a10=r10+0.5
      b10=0.5+n10-r10
      a00=r00+0.5
      b00=0.5+n00-r00
      
      s2=rbinom(1,s,mis)  ## sample size for marker missing samples
      s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
      s0=s-s2-s1 ## sample size for neg samples
      
      n211=rbinom(1,s2,rate21) ## marker missing samples in treatment
      n21=n21+n211
      n20=n20+s2-n211 ## missing samples in control
      r21=r21+rbinom(1,n211,p21) ## response missing samples in trt
      r20=r20+rbinom(1,s2-n211,p20) ## response missing samples in ctrl
      
      n111 = rbinom(1,s1,rate11) ## response adaptive randomization
      n11=n11+n111 ## marker (+) samples in trt
      n10=n10+s1-n111 ## marker (+) samples in ctrl
      r11=r11+rbinom(1,n111,p11) ## response samples (+) in trt
      r10=r10+rbinom(1,s1-n111,p10) ## response samples (+) in ctrl
      
      n011 = rbinom(1,s0,rate01) ## response adaptive randomization
      n01=n01+n011 ## marker (-) samples in trt
      n00=n00+s0-n011 ## marker (-) samples in ctrl
      r01=r01+rbinom(1,n011,p01) ## response samples (-) in trt
      r00=r00+rbinom(1,s0-n011,p00) ## response samples (-) in ctrl
      
    }
    
    for(i in 5: (t-1)){
      
      a21=r21+0.5
      b21=0.5+n21-r21
      a20=r20+0.5
      b20=0.5+n20-r20
      
      a11=r11+0.5
      b11=0.5+n11-r11
      a01=r01+0.5
      b01=0.5+n01-r01
      a10=r10+0.5
      b10=0.5+n10-r10
      a00=r00+0.5
      b00=0.5+n00-r00
      cbeta11 = rbeta(10000,a11,b11)
      cbeta10 = rbeta(10000,a10,b10)
      cbeta01 = rbeta(10000,a01,b01)
      cbeta00 = rbeta(10000,a00,b00)
      pbeta1 = mean((cbeta11-cbeta01)>epi1) #pre
      pbeta0 = mean(abs(cbeta10-cbeta00)>epi2) #pro
      
      model5 <- jags.model(textConnection(RatonMCMC), 
                           data = list( n21=n21,n20=n20,n11=n11,n10=n10,n00=n00,n01=n01,
                                        r21=r21,r20=r20,r11=r11,r01=r01,r00=r00,r10=r10,
                                        pbeta1=pbeta1,pbeta0=pbeta0,
                                        phi=prevelence),quiet=T)
      update(model5, 1000, progress.bar="none")
      re <- jags.samples(model5, n.iter=10000,variable.names=c("p21", "p20","p11", "p10","p01", "p00"), progress.bar="none")
      
      p21_array <- re$p21[1:10000]
      p20_array <- re$p20[1:10000]
      p11_array <- re$p11[1:10000]
      p10_array <- re$p10[1:10000]
      p01_array <- re$p01[1:10000]
      p00_array <- re$p00[1:10000]
      
      
      if( mean((p11_array - p10_array)>epi3)<pil & mean((p01_array - p00_array)>epi3)<pil){
        earlystop=1
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01+n21
        ss_ctrall[l] = n10+n00+n20
        ss_all[l] = n11+n10+n01+n00+n21+n20
        decision[l] = 0 ## early stop for futility for (+) & (-)
        D_p_p[l] = 0
        D_n_p[l] = 0
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
        break
      }
      else if (mean((p01_array - p00_array)>epi3)>pis & mean((p11_array - p10_array)>epi3)>pis){
        earlystop=1
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01
        ss_ctrall[l] = n10+n00
        ss_all[l] = n11+n10+n01+n00
        decision[l] = 1 ## early stop for superiority for (-) & (+)
        D_p_p[l] = 1
        D_n_p[l] = 1
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
        break
        
      } else if (mean((p11_array - p10_array)>epi3)>pis & mean((p01_array - p00_array)>epi3)<pil  ){
        earlystop=1
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01
        ss_ctrall[l] = n10+n00
        ss_all[l] = n11+n10+n01+n00
        decision[l] = 11 ## early stop for superiority for (+) & futility for (-)
        D_p_p[l] = 1
        D_n_p[l] = 0
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
        break
        
      } else if (mean((p11_array - p10_array)>epi3)<pil & mean((p01_array - p00_array)>epi3)>pis  ){
        earlystop=1
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01
        ss_ctrall[l] = n10+n00
        ss_all[l] = n11+n10+n01+n00
        decision[l] = 21 ## early stop for superiority for (-) & futility for (+)
        D_p_p[l] = 0
        D_n_p[l] = 1
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
        break
        
      }
      
      ## early stop for futility for (+), recruit (-) only,not recruit marker missing 
      
      else if (mean((p11_array - p10_array)>epi3)<pil){
        earlystop=0
        s2=rbinom(1,s,mis)  ## sample size for marker missing samples
        s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
        s0=s-s2-s1 ## sample size for neg samples
        
        lambda <- lambda2
        
        rate01 = mean(p01_array> p00_array)^lambda/(mean(p01_array> p00_array)^lambda+mean(p01_array < p00_array)^lambda)
        n011 = rbinom(1,s0,rate01) ## response adaptive randomization
        n01=n01+n011 ## marker (-) samples in trt
        n00=n00+s0-n011 ## marker (-) samples in ctrl
        r01=r01+rbinom(1,n011,p01) ## response samples (-) in trt
        r00=r00+rbinom(1,s0-n011,p00) ## response samples (-) in ctrl
        
      }
      
      ## early stop for futility for (-), recruit (+) only,not recruit marker missing 
      
      else if (mean((p01_array - p00_array)>epi3)<pil){
        earlystop=0
        s2=rbinom(1,s,mis)  ## sample size for marker missing samples
        s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
        s0=s-s2-s1 ## sample size for neg samples
        
        lambda <- lambda2
        
        rate11 = mean(p11_array> p10_array)^lambda/(mean(p11_array> p10_array)^lambda+mean(p11_array < p10_array)^lambda)
        n111 = rbinom(1,s1,rate11) ## response adaptive randomization
        n11=n11+n111 ## marker (+) samples in trt
        n10=n10+s1-n111 ## marker (+) samples in ctrl
        r11=r11+rbinom(1,n111,p11) ## response samples (+) in trt
        r10=r10+rbinom(1,s1-n111,p10) ## response samples (+) in ctrl
        
      }
      
      ## early stop for superiority for (+), recruit (-) only
      
      else if (mean((p11_array - p10_array)>epi3)>pis){
        earlystop=0
        s2=rbinom(1,s,mis)  ## sample size for marker missing samples
        s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
        s0=s-s2-s1 ## sample size for neg samples
        
        lambda <- lambda2
        
        rate01 = mean(p01_array> p00_array)^lambda/(mean(p01_array> p00_array)^lambda+mean(p01_array < p00_array)^lambda)
        n011 = rbinom(1,s0,rate01) ## response adaptive randomization
        n01=n01+n011 ## marker (-) samples in trt
        n00=n00+s0-n011 ## marker (-) samples in ctrl
        r01=r01+rbinom(1,n011,p01) ## response samples (-) in trt
        r00=r00+rbinom(1,s0-n011,p00) ## response samples (-) in ctrl
        
      }
      
      ## early stop for superiority for (-), recruit (+) only
      
      else if (mean((p01_array - p00_array)>epi3)>pis ){
        earlystop=0
        s2=rbinom(1,s,mis)  ## sample size for marker missing samples
        s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
        s0=s-s2-s1 ## sample size for neg samples
        
        lambda <- lambda2
        
        rate11 = mean(p11_array> p10_array)^lambda/(mean(p11_array> p10_array)^lambda+mean(p11_array < p10_array)^lambda)
        n111 = rbinom(1,s1,rate11) ## response adaptive randomization
        n11=n11+n111 ## marker (+) samples in trt
        n10=n10+s1-n111 ## marker (+) samples in ctrl
        r11=r11+rbinom(1,n111,p11) ## response samples (+) in trt
        r10=r10+rbinom(1,s1-n111,p10) ## response samples (+) in ctrl
        
      }
      
      else { ##  neither (+) nor (-) has evidence of futility or superiority at this interim
        earlystop=0
        
        s2=rbinom(1,s,mis)  ## sample size for marker missing samples
        s1=rbinom(1,s-s2,prevelence) ## sample size for pos samples
        s0=s-s2-s1 ## sample size for neg samples
        
        lambda <- lambda2
        
        rate21 = mean(p21_array> p20_array)^lambda/(mean(p21_array> p20_array)^lambda+mean(p21_array < p20_array)^lambda)
        n211 = rbinom(1,s2,rate21) ## response adaptive randomization
        n21=n21+n211 ## marker missing samples in trt
        n20=n20+s2-n211 ## marker missing samples in ctrl
        r21=r21+rbinom(1,n211,p21) ## response samples missing in trt
        r20=r20+rbinom(1,s2-n211,p20) ## response samples missing in ctrl
        
        rate11 = mean(p11_array> p10_array)^lambda/(mean(p11_array> p10_array)^lambda+mean(p11_array < p10_array)^lambda)
        n111 = rbinom(1,s1,rate11) ## response adaptive randomization
        n11=n11+n111 ## marker (+) samples in trt
        n10=n10+s1-n111 ## marker (+) samples in ctrl
        r11=r11+rbinom(1,n111,p11) ## response samples (+) in trt
        r10=r10+rbinom(1,s1-n111,p10) ## response samples (+) in ctrl
        
        
        rate01 = mean(p01_array> p00_array)^lambda/(mean(p01_array> p00_array)^lambda+mean(p01_array < p00_array)^lambda)
        n011 = rbinom(1,s0,rate01) ## response adaptive randomization
        n01=n01+n011 ## marker (-) samples in trt
        n00=n00+s0-n011 ## marker (-) samples in ctrl
        r01=r01+rbinom(1,n011,p01) ## response samples (-) in trt
        r00=r00+rbinom(1,s0-n011,p00) ## response samples (-) in ctrl
        
      } 
    } # end of interim
    
    if (earlystop==0)   { 
      # reach maximum sample size at end of interims
      if(mean((p11_array - p10_array)>epi3)<pif & mean((p01_array - p00_array)>epi3)<pif){ # final analysis
        decision[l]=5 # neither pos nor neg is promising
        D_p_p[l] = 0
        D_n_p[l] = 0
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01+n21
        ss_ctrall[l] = n10+n00+n20
        ss_all[l] = n11+n10+n01+n00+n21+n20
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
      }
      else if (mean((p01_array - p00_array)>epi3)<pif) {
        decision[l]=6 # pos promising but neg are not promising
        D_p_p[l] = 1
        D_n_p[l] = 0
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01+n21
        ss_ctrall[l] = n10+n00+n20
        ss_all[l] = n11+n10+n01+n00+n21+n20
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
      }
      else if (mean((p11_array - p10_array)>epi3)<pif) {
        decision[l]=66 # neg promising but pos are not promising
        D_p_p[l] = 0
        D_n_p[l] = 1
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01+n21
        ss_ctrall[l] = n10+n00+n20
        ss_all[l] = n11+n10+n01+n00+n21+n20
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
      }
      else {
        decision[l]=7 # both pos and neg are promising
        D_p_p[l] = 1
        D_n_p[l] = 1
        ss_pos[l] = n11
        ss_neg[l] = n01
        ss_pos_1[l] = n10
        ss_neg_1[l] = n00
        ss_posall[l] = n11+n10
        ss_negall[l] = n01+n00
        ss_trtall[l] = n11+n01+n21
        ss_ctrall[l] = n10+n00+n20
        ss_all[l] = n11+n10+n01+n00+n21+n20
        p21_postr[l] = r21/n21
        p20_postr[l] = r20/n20
        p11_postr[l] = r11/n11
        p10_postr[l] = r10/n10
        p01_postr[l] = r01/n01
        p00_postr[l] = r00/n00
        
        p21_post[l] = mean(p21_array)
        p20_post[l] = mean(p20_array)
        p11_post[l] = mean(p11_array)
        p10_post[l] = mean(p10_array)
        p01_post[l] = mean(p01_array)
        p00_post[l] = mean(p00_array)
        
        rate21l[l]=rate21
        rate11l[l]=rate11
        rate01l[l]=rate01
      }
    } 
  } # end of simulation
  decisions <- list(table(decision))
  return(c("Raton",round(mean(D_p_p)*100,digits=2) ,round(mean(D_n_p)*100,digits=2) ,
           round(mean(ss_pos),digits=2),
           round(mean(ss_neg),digits=2),
           round(mean(ss_pos_1),digits=2),
           round(mean(ss_neg_1),digits=2),
           round(mean(ss_posall),digits=2),
           round(mean(ss_negall),digits=2),
           round(mean(ss_trtall),digits=2),
           round(mean(ss_ctrall),digits=2),
           round(mean(ss_all),digits=2),
           round(mean(p20_postr),digits=2),
           round(mean(p21_postr),digits=2),
           round(mean(p10_postr),digits=2),
           round(mean(p11_postr),digits=2),
           round(mean(p00_postr),digits=2),
           round(mean(p01_postr),digits=2),
           round(mean(p20_post),digits=2),
           round(mean(p21_post),digits=2),
           round(mean(p10_post),digits=2),
           round(mean(p11_post),digits=2),
           round(mean(p00_post),digits=2),
           round(mean(p01_post),digits=2),
           round(mean(rate11l),digits=2),
           round(mean(rate01l),digits=2),
           round(mean(rate21l),digits=2),
           decisions
  ))
}











