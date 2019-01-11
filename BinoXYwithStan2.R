#setwd("~/Documents/courbariaux/WORKLUC/BinoXY")
rm(list=ls())
library(rstan)
rstan_options(auto_write = TRUE)
knitr::opts_chunk$set(echo = TRUE)
options(digits=3,width=80)
library(tidyverse)
library(ggplot2)
seuil=0.1 #mn
printemps = (31+28+15):(31+28+31+30+31+15)
saison = printemps
horizon = 4 #jours
K=50 #membres
datecharniere="2011-01-01"
base<-function(t){
  cbind(0*t+1,sin(2*pi*t/365),cos(2*pi*t/365),
        sin(2*pi*2*t/365),cos(2*pi*2*t/365),
        sin(2*pi*3*t/365),cos(2*pi*3*t/365))
}


## ----ChargementHistorique------------------------------------------------
load("AINhist.Rdata")
bisextil=grep(pattern = '-02-29',x=t.hist)
data.histo<-tibble(Date=t.hist[-bisextil],y=y.hist[-bisextil],
                   calendaire=rep(1:365,length(Date)/365))

data.histo %>%  group_by(calendaire) %>% 
  summarize(probapluie=mean(y>seuil)) %>% 
  ggplot(mapping = aes(x = calendaire, y=probapluie))+geom_point()+geom_smooth()


data.histo %>%  filter(y>1*seuil) %>% 
  group_by(calendaire) %>%  
  summarize(pluiemean=mean(y),pluiestd=var(y)^0.5) ->datahistomoyen

data.histo %>%  filter(y>1*seuil) %>% 
  ggplot(mapping = aes(x = calendaire, y=y)) + geom_point()+
  geom_point(data=datahistomoyen,mapping = aes(x = calendaire, y=pluiemean),col='red')+
  geom_smooth(data=datahistomoyen,aes(y=pluiemean),se=0,color='blue')+
  geom_smooth(data=datahistomoyen,aes(y=pluiemean+2*pluiestd),se=0,color='blue')

## ----ChargeVouglans20052008----------------------------------------------
load("AIN20052008.Rdata")
nbech20052008=max(data$Echeance)
training20052008o=as.tibble(data)
bisextil=grep(pattern = '-02-29',x=training20052008o$Date)
l=1:dim(training20052008o)[1]
l=(-1+l)%/%nbech20052008
numjour=(l+(training20052008o$Echeance))
numjour[training20052008o$Date>"2008-12-31"]=numjour[training20052008o$Date>"2008-12-31"]-1
training20052008o$calendaire=((-1+(numjour))%%365)+1
tail(training20052008o)

## ----Charge200112015-----------------------------------------------------
load("AIN20112015.Rdata") 
##TESTER Completude base
# data %>% group_by(Echeance) %>% summarize(n())
data$Date[1]->debut
nbech=max(unique(data$Echeance))
nbjours=dim(data)[1]/nbech
rythme=rep(1:nbech,nbjours)+rep(0:(nbjours-1),each=nbech)
DateCheck=debut+rythme-1
data %>% mutate(Check=(Date!=DateCheck)) %>% 
  group_by(Echeance) %>% summarize(n(),sum(Check)) 
#pour verifier que tout est OK
which(data$Date!=DateCheck)->pb
# A comparer data$Date[pb] et DateCheck[pb]
data$Date[pb]<-DateCheck[pb] 

## ----Charge200112015suite------------------------------------------------
validation20112015o=as.tibble(data)
bisextil=grep(pattern = '-02-29',x=validation20112015o$Date)
l=1:dim(validation20112015o)[1]
l=(-1+l)%/%nbech
numjour=(l+(validation20112015o$Echeance))
numjour[validation20112015o$Date>"2012-12-31"]=numjour[validation20112015o$Date>"2012-12-31"]-1
validation20112015o$calendaire=((-1+(numjour))%%365)+1
tail(validation20112015o)
tail(training20052008o)
head(validation20112015o)
## ----fairebasetotale-----------------------------------------------------
########################################################
colnames(validation20112015o)=colnames(training20052008o)
validation20112015o %>% filter(Echeance<=nbech20052008) %>% 
  rbind(training20052008o,.) ->datao
#Nettoyage
memoire=c("seuil","horizon", "K", "printemps","base",
          "data.histo","datao")
LS=ls(); listeajeter=NULL
for (elementajeter in LS){if (!(elementajeter %in% memoire))
  listeajeter=c(listeajeter,elementajeter) }
rm(list=c(listeajeter,"listeajeter","elementajeter","LS"))

datecharniere="2011-01-01"

# datao %>% filter(Date>=datecharniere)->validation20112015o
# datao %>% filter(Date<datecharniere)->training20052008o

## ----essai---------------------------------------------------------------
is.nonzero<-function(x){ifelse(x>seuil,1,0)}
take.log<-function(x){ifelse(x>0,log(x),0)}

datao %>% 
  mutate(NX=rowSums(is.nonzero(select(.,starts_with("CEP")))),
         SX=rowSums(select(.,starts_with("CEP")))) %>% 
  mutate_at(vars(matches("CEP")), take.log) %>% 
  mutate(SLX=rowSums(select(.,starts_with("CEP"))))   %>% 
  select(-starts_with("CEP")) %>% 
  filter(Echeance==horizon) %>% 
  arrange(NX) %>% 
  filter(!is.na(NX),calendaire %in% printemps)->dtout

dtout %>% filter(Date<datecharniere, Obs>seuil)->d
dtout %>% filter(Date<datecharniere, Obs<=seuil)->d
 
 plot(log(1+d$SX),d$SLX, pch='+')
#########################################################
model_string<-"
functions {
  // Function that returns the log pdf of the Exchangeable BinoGamma model 
 real ExchBinoGamma_lp( int Nx, real Sx, real SLx, 
 real alpha, real gam) {
 return((gam-1)*SLx -alpha*Sx+gam*Nx*log(alpha)-Nx*lgamma(gam));
 }
 }
data {
int<lower=0> K;         // size of the ensemble 50
int<lower=0> T;         // number of observations 
int NX[T];         // number of non zero members 
real SX[T];    // nonzero sum of members 
real SLX[T];    // nonzero sum of members 
}
parameters {
real<lower=0> a_NX;      // first coeff of beta 
real<lower=0> b_NX;      // second coeff of beta
real<lower=0,upper=1> pi[T];   // beta distributed effect by time
real<lower=0> gamma_SX;      // second coeff of gamma
real<lower=0> alpha_SX[T];      // first coeff of gamma
real<lower=0> g_alpha_SX;      // first coeff of gamma gamma
real<lower=0> c_alpha_SX;      // second coeff of gamma gamma
}

model {
target += binomial_lpmf(NX | K, pi); // log-likelihood
target += beta_lpdf(pi | a_NX, b_NX);       // random effect
target += gamma_lpdf(alpha_SX | g_alpha_SX, c_alpha_SX);       // random effect
target += uniform_lpdf(a_NX | 0, K); // prior
target += uniform_lpdf(b_NX | 0, K); // prior
target += exponential_lpdf(gamma_SX | 1); // prior
target += uniform_lpdf(g_alpha_SX | 0, 1000*K); // prior
target += uniform_lpdf(c_alpha_SX | 0, 1000*K); // prior
  for (t in 1:T){
target += ExchBinoGamma_lp(NX[t], SX[t], SLX[t], 
alpha_SX[t], gamma_SX);
}
}
"
NX_dat <- list(K = 50, 
                    NX = d$NX,
                    T = length(d$NX),
                    SX = d$SX,
                    SLX = d$SLX
               )
fit <- stan(model_code=model_string, data = NX_dat)
print(fit)
#plot(fit)
pairs(fit, pars = c("a_NX", "b_NX","gamma_SX", "g_alpha_SX", "c_alpha_SX"))
fit->fit0
fit->fit1
save(fit0, fit1, file = "fitsPrintemps.RData")
load(file="fitsPrintemps.RData")
d=dtout
#rm(dtout)
LogExchBinoGamma<-function(Nx, Sx, SLx, p, alpha, gam) {
  res=(Nx*log(p)+(K-Nx)*log(1-p)+
         (gam-1)*SLx -alpha*Sx+
         gam*Nx*log(alpha)-Nx*lgamma(gam))
}
simu0<-rstan::extract(fit0, pars = c("a_NX", "b_NX","g_alpha_SX","c_alpha_SX","gamma_SX" ))

simu0$p<-rbeta(4000,simu0$a_NX,simu0$b_NX)
simu0$alpha<-rgamma(4000,simu0$g_alpha_SX,simu0$c_alpha_SX)
AjoutLogProb0<-function(i){
  test=LogExchBinoGamma(as.numeric(d[i,"NX"]), 
                        as.numeric(d[i,"SX"]),
                        as.numeric(d[i,"SLX"]),
                        simu0$p, simu0$alpha, simu0$gamma_SX) 
  LogProb=max(test,na.rm = T)+log(mean(exp(test-max(test,na.rm = T)),na.rm=T))
}
sapply(1:length(d$NX),AjoutLogProb0)->d$LogProb0
############
simu1<-rstan::extract(fit1, pars = c("a_NX", "b_NX","g_alpha_SX","c_alpha_SX","gamma_SX" ))
simu1$p<-rbeta(4000,simu1$a_NX,simu1$b_NX)
simu1$alpha<-rgamma(4000,simu1$g_alpha_SX,simu1$c_alpha_SX)
AjoutLogProb1<-function(i){
  test=LogExchBinoGamma(as.numeric(d[i,"NX"]), 
                        as.numeric(d[i,"SX"]),
                        as.numeric(d[i,"SLX"]),
                        simu1$p, simu1$alpha, simu1$gamma_SX) 
  LogProb=max(test,na.rm = T)+log(mean(exp(test-max(test,na.rm = T)),na.rm=T))
}
sapply(1:length(d$NX),AjoutLogProb1)->d$LogProb1
#############

## ----calculerho----------------------------------------------------------
dtout=d
d %>% filter(Date<datecharniere)->d
ro<-0.5
tol<-0.00001
repeat {
  pp<-ro*exp(d$LogProb1)/(ro*exp(d$LogProb1)+(1-ro)*exp(d$LogProb0))
  ronew<-mean(pp)
  print(c(ro,ronew))
  if (abs(ro-ronew) <= tol){
    break
  }
  ro <- ronew
}
sum(d$Obs>seuil)/length(d$Obs)
ro
d$probazeq1 <- ro*exp(d$LogProb1)/(ro*exp(d$LogProb1)+(1-ro)*exp(d$LogProb0))

zy <- matrix(
  c(100*mean((d$Obs>seuil)*d$probazeq1),100*(mean((1-(d$Obs>seuil))*d$probazeq1)),
    100*mean((d$Obs>seuil)*(1-d$probazeq1)),100*(mean((1-(d$Obs>seuil))*(1-d$probazeq1)))),    byrow=T,nr=2,nc=2 )
row.names(zy)<-c("z=1","z=0")
colnames(zy)<-c("y=1","y=0")
print(zy)
# z en ligne y en colonne
probazeg1sachantyeg1=zy[1,1]/(zy[1,1]+zy[2,1])
probazeg1sachantyeg0=zy[1,2]/(zy[1,2]+zy[2,2])
probayeg1sachantzeg1=zy[1,1]/(zy[1,1]+zy[1,2])
probayeg1sachantzeg0=zy[1,2]/(zy[2,1]+zy[2,2])

library(xtable)
tabzy= xtable(zy)
print(tabzy)
print(zy)

###
dtout$probazeq1 <- ro*exp(dtout$LogProb1)/(ro*exp(dtout$LogProb1)+(1-ro)*exp(dtout$LogProb0))
dtout$probayeq1<-probayeg1sachantzeg1*dtout$probazeq1+
  probayeg1sachantzeg0*(1-dtout$probazeq1)
dtout %>% filter(Date>=datecharniere)->d
plot(d$probayeq1,d$NX/K, col=1+(d$Obs>=seuil))
library(ROCR)
pred<-prediction(d$probayeq1,d$Obs>seuil)
perf <- performance(pred,'tpr','fpr')
par(mfrow=c(1,1))
plot(perf)
pred<-prediction(d$NX/K,d$Obs>seuil)
perf <- performance(pred,'tpr','fpr')
plot(perf,add=T,col='blue')
mean((d$NX/K-(d$Obs>seuil))^2)
mean((d$probayeq1-(d$Obs>seuil))^2)
library(verification)
verif<-verify(obs = 1*(d$Obs>seuil), pred=d$probayeq1, obs.type = "binary",pred.type='prob')
plot(verif)
verif<-verify(obs = 1*(d$Obs>seuil), pred=d$NX/K, obs.type = "binary",pred.type='prob')
plot(verif)
