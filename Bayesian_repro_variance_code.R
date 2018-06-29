##Bayesian model of variance for sea otters, southern elephant seals, and Pacific harbor seals 
##Prepare by M Tarjan
##April 28, 2016


##LOAD PACKAGES
library(gridExtra) ##required for grid.arrange
library(ggplot2)

##PREP DATA
##folder with data files
folder<-"C:/Users/max/Desktop/Tarjan/otter data/"
library(stringr)

##PREP SEA OTTER DATA
##create list of potential sampled males and matrix for years they were active
##have ped.dat.y, info, and n.prop.samp
bd.dates<-read.csv(str_c(folder, "birth_death_dates.csv"))
#bd.dates<-read.csv("birth_death_dates.csv") ##birth and death dates of all animals; different from info because I assume here that males are dead if there is no confirmation of life (only alive if confirmed alive by capture/resight)
#prop.samp<-read.csv("/Users/MaxTarjan/Dropbox/Thesis/Genetics_chapter2/Tim's ch2 repro variance analysis/prop_samp_28Apr16.csv")
prop.samp<-read.csv(str_c(folder, "prop_samp_28Apr16.csv"))
#ped.dat.y<-read.csv('pedigraph_data_years_5Oct15.csv')
ped.dat.y<-read.csv(str_c(folder,'pedigraph_data_years_5Oct15.csv'))
##terms needed for model
##(nmales, Nyrs, scalefact, Sires, proppup, propdad, Npuptot, obs, Nm, M, maleyrs, pi)
##nmales = value. number of potential fathers that were sampled
##Nyrs == value. number of years
##scalefact = proportion of "accounted for" males. Tim set this to 0.9 ##saying estimate of pop might be 10% above estimate
##Sires: vector of length year- number of males that were assigned a pup in each year??
##proppup & propdad: prop pup and prop dad by year (each a vector of length = number of years)
##Npuptot: total number of pups by year (vector)
##obs: matrix- sampled males and years they were active (Males matrix)
##Nm ==vector of length years ##number of males active each year
##M = matrix of line numbers of males matrix that are active each year ##indexing the rows of obs for males that are active ##0s when not full. each column is a year
##maleyrs = vector with length number of males (number of years that male was active)
##pi ##not in model
##maleyrsmean = mean number of years that a male of this species is reproductively active (ie alive to reproduce), given that he has already reached reproductive maturity (age 5, maybe 4...still deciding)

##create obs
ped.dat.pup<-ped.dat.y[which(ped.dat.y$sire!=0 | ped.dat.y$dam!=0),] ##restrict pedigraph data to those with parents
ped.dat.pup<-ped.dat.pup[which(ped.dat.pup$sex=='F'),] ##restrict pedigraph data to female pups
ymax<-max(ped.dat.pup$year); ymin<-min(ped.dat.pup$year)
years<-ymin:ymax
Nyrs<-length(years)

n<-0
obs<-matrix(ncol = Nyrs, nrow = length(which(bd.dates$sex=='M')))
male.bd<-bd.dates[which(bd.dates$sex=='M'),]
for (y in 1:Nyrs) {
  for (m in 1:nrow(obs)) {
    if (male.bd$birth[m]+5 <= years[y] & male.bd$death[m] >= years[y]) { ##males must be 5 years of age to be considered "active"
      obs[m,y]<-0
    }
    else {
      obs[m,y]<--1
    }
    
    ##add number of pups sired by each male
    sires.temp<-length(ped.dat.pup$sire[which(ped.dat.pup$sire==male.bd$ID[m] & ped.dat.pup$year==years[y])])
    if (sires.temp>0) {
      if (obs[m,y]==-1) {print('error')}
      n<-n+sires.temp
      obs[m,y]<-obs[m,y]+sires.temp
    }
  }
}

sum<-apply(X = obs, MARGIN = 1, FUN = sum)
sum.min<--1*ncol(obs)
row.no.yrs<-which(sum==sum.min) #rows containing males with no years of siring
obs<-obs[-c(row.no.yrs),]

males<-bd.dates$ID[which(bd.dates$sex=='M')]; males<-males[-c(row.no.yrs)]
nmales<-length(males)

obs2<-obs; obs2[obs2==-1]<-0

Sires<-apply(X = obs2, MARGIN = 2, FUN = sum)

proppup<-prop.samp$prop_pup
propdad<-prop.samp$prop_dad
Npuptot<-prop.samp$Pup_pop
Npupsamp<-prop.samp$n_pups

##number of pups sampled each year; pulled from provided table above
#for (j in 1:length(years)) {
#  pupsamp<-length(which(bd.dates$sex=='F' & bd.dates$birth-1==years[j]))
#  puppop<-n.prop.samp$Pup_pop[which(n.prop.samp$year==years[j])]
#  proppup<-c(proppup, pupsamp/puppop)
#  dadsamp<-length(which(obs[,j]!=-1))
#  print(dadsamp)
#  dadpop<-n.prop.samp$dad_pop[which(n.prop.samp$year==years[j])]
#  propdad<-c(propdad, dadsamp/dadpop)
#}

##maleyrs = vector with length number of males (number of years that male was active)
maleyrs<-dim(0)
for (j in 1:nrow(obs)) {
  maleyrs<-c(maleyrs, length(which(obs[j,]!=-1)))
}

##maleyrsmean
maleyrsmean<-3

##OVERWRITE MALEYRS WITH MEAN REPRO LIFETIME
maleyrs<-rep(maleyrsmean, length(maleyrs))

##M = matrix of line numbers of males matrix that are active each year; columns are years
M<-matrix(nrow=nrow(obs), ncol=ncol(obs))
Nm<-dim(0) ##number of males active each year
for (j in 1:length(years)) {
  M[1:length(which(obs[,j]!=-1)),j]<-which(obs[,j]!=-1)
  Nm<-c(Nm, length(which(obs[,j]!=-1)))
}

Bayes.in<-list()
Bayes.in$nmales<-nmales
Bayes.in$Nyrs<-Nyrs
#Bayes.in$scalefact<-0.9
Bayes.in$Sires<-Sires
Bayes.in$proppup<-proppup
Bayes.in$propdad<-propdad
#Bayes.in$Npuptot<-Npuptot ##number of pups in population
Bayes.in$obs<-obs
Bayes.in$Nm<-Nm
Bayes.in$M<-M
Bayes.in$maleyrs<-maleyrs
Bayes.in$Npupsamp<-Npupsamp #number of pups sampled in each year
#Bayes.in$maleyrsmean<-maleyrsmean

names(Bayes.in)

for (j in 1:length(Bayes.in)) {
  write.csv(Bayes.in[[j]], str_c(folder, "Bayes_in/", names(Bayes.in)[j], ".csv"), row.names=F)
}

##BAYESIAN MODEL
##SEA OTTER BAYESIAN MODEL
#library(rjags)
library(jagsUI) ##required for jags function
#var.names<-c('S1', 'S2', 'mu', 'sig', 'bta1', 'bta2', 'TSires')
var.names<-c('S1', 'S2', 'mu', 'sig', 'bta', 'TSires')
burnin<-10000
nsamps<-20000
##run at 20,000 and 10,000
nchains<-3
model <- jags(data = Bayes.in, parameters.to.save = var.names, model.file = "Reproskew_JAGS.txt", n.chains = nchains, n.iter = nsamps, n.burnin = burnin)

model

##PLOTS
resolution<-96

#tiff("Figures_genetics/index1_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
plot(density(model$sims.list$S1), xlab='Index 1', col="black", lwd=2, main=""); 
abline(v = model$q2.5$S1, lty="dashed", lwd=2); abline(v = model$q97.5$S1, lty="dashed", lwd=2); abline(v = model$q50$S1, lwd=2)
#dev.off()

#tiff("Figures_genetics/index2_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
plot(density(model$sims.list$S2), xlab='Index 2', col="black", lwd=2, main=""); 
abline(v = model$q2.5$S2, lty="dashed", lwd=2); abline(v = model$q97.5$S2, lty="dashed", lwd=2); abline(v = model$q50$S2, lwd=2)
#dev.off()

##plot frequency distrib of male repro success
TSires<-model$sims.list$TSires
Nreps<-nrow(TSires)
xi <- 100 #seq(0,max(TSires), 1)
Tsiresdist<-matrix(nrow=Nreps, ncol=xi)
for (j in 1:Nreps) {
  Tsiresdist[j,] <-density(TSires[j,], n=xi, from=0, to=20)$y
  x.temp<-density(TSires[j,], n=xi, from=0, to=20)$x
}
Tsiresdist_m<-apply(Tsiresdist, MARGIN = 2, FUN = mean)
Tsiresdist_CI<-apply(Tsiresdist, MARGIN = 2, FUN = quantile, probs = c(0.05, 0.95))
#xi by Tsiresdist_m

#tiff("Figures_genetics/siring_freq_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
plot(x = x.temp, y = Tsiresdist_m, type="l", ylim=c(0.009, max(Tsiresdist_CI)), xlab='Surviving female pups sired', ylab='Proportion of Males', xlim=c(0.7,19.27))
polygon(x = c(x.temp, rev(x.temp)), y = c(Tsiresdist_CI[1,], rev(Tsiresdist_CI[2,])), col="light grey", border = NA)
lines(x = x.temp, y = Tsiresdist_m, lwd=2)
#dev.off()

model.so<-model
x.temp.so<-x.temp
Tsiresdist_m.so<-Tsiresdist_m
Tsiresdist_CI.so<-Tsiresdist_CI

##PREP E SEAL DATA
es.dat<-read.csv(str_c(folder, "southern_eseal_data.csv"))
##sampled pups c(115, 77)
##sampled males c(78, 62)
##prop pups (54+90)/2

#nmales<-46
nmales<-length(unique(es.dat$male))
Nyrs<-2
Sires<-c(sum(es.dat$paternities[which(es.dat$year==1996)]), sum(es.dat$paternities[which(es.dat$year==1997)])) ##number of pups assigned per year
proppup<-c(0.72, 0.72)
propdad<-c(0.95, 0.95)
Npuptot<-c(115*100/72, 77*100/72)
Npupsamp<-c(115, 77)
Nm<-c(length(unique(es.dat$male[which(es.dat$year==1996)])), length(unique(es.dat$male[which(es.dat$year==1997)])))

obs<-matrix(nrow=length(unique(es.dat$male)), ncol=length(unique(es.dat$year)))
years<-c(1996, 1997)
for (m in 1:nrow(obs)) {
  for (y in 1:length(years)) {
    if (length(which(es.dat$male==es.dat$male[m] & es.dat$year==years[y]))==0) {
      obs[m,y]<--1
    } else {
      obs[m,y]<-sum(es.dat$paternities[which(es.dat$male==es.dat$male[m]& es.dat$year==years[y])])
    }
  }
}

##maleyrs = vector with length number of males (number of years that male was active)
maleyrs<-dim(0)
for (j in 1:nrow(obs)) {
  maleyrs<-c(maleyrs, length(which(obs[j,]!=-1)))
}

##M = matrix of line numbers of males matrix that are active each year; columns are years
M<-matrix(nrow=nrow(obs), ncol=ncol(obs))
Nm<-dim(0) ##number of males active each year
for (j in 1:length(years)) {
  M[1:length(which(obs[,j]!=-1)),j]<-which(obs[,j]!=-1)
  Nm<-c(Nm, length(which(obs[,j]!=-1)))
}

##ESTIMATE REPRO LIFETIME
age.temp<-c(5,6,8, 21, 21, 21) ##age data. must be entered manually. start at first age of repro
n.temp<-c(31.7,26.5,14,0, 0, 0) ##number of males that survive to that age
n.temp<-n.temp*1000
plot(age.temp, n.temp)
func.temp<-nls(formula = n.temp~a*age.temp^b, start = list(a=10000,b=-1))
curve(coef(func.temp)[1]*x^coef(func.temp)[2], add=T)

out.temp<-dim(0)
for (j in 5:19) { ##from min age to max age
  n.j<-coef(func.temp)[1]*j^coef(func.temp)[2] ##animals that survive to j
  n.j1<-coef(func.temp)[1]*(j+1)^coef(func.temp)[2] ##animals that survive to j+1
  n.new<-round(n.j-n.j1, digits = 0) ##new animals that died between j and j+1
  span.temp<-rep(j-(5-1),n.new) ##rep repro lifetime for that many males
  #print(span.temp); print(j)
  out.temp<-c(out.temp, span.temp)
}
mean.lifetime<-round(mean(out.temp),0) ##mean repro lifetime

##OVERWRTIE MALEYRS WITH 4 FOR EVERY MALE
maleyrs<-rep(mean.lifetime, length(maleyrs))

Bayes.in<-list()
Bayes.in$nmales<-nmales
Bayes.in$Nyrs<-Nyrs
#Bayes.in$scalefact<-0.9
Bayes.in$Sires<-Sires
Bayes.in$proppup<-proppup
Bayes.in$propdad<-propdad
#Bayes.in$Npuptot<-Npuptot
Bayes.in$Npupsamp<-Npupsamp
Bayes.in$obs<-obs
Bayes.in$Nm<-Nm
Bayes.in$M<-M
Bayes.in$maleyrs<-maleyrs

##E SEAL BAYESIAN MODEL
var.names<-c('S1', 'S2', 'mu', 'sig', 'bta', 'TSires')
burnin<-10000
nsamps<-20000
##run at 20,000 and 10,000
nchains<-3
model <- jags(data = Bayes.in, parameters.to.save = var.names, model.file = "Reproskew_JAGS.txt", n.chains = nchains, n.iter = nsamps, n.burnin = burnin)


##PLOTS
#resolution<-96

#tiff("Figures_genetics/index1_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
plot(density(model$sims.list$S1), xlab='Index 1', col="black", lwd=2, main=""); 
abline(v = model$q2.5$S1, lty="dashed", lwd=2); abline(v = model$q97.5$S1, lty="dashed", lwd=2); abline(v = model$q50$S1, lwd=2)
#dev.off()

#tiff("Figures_genetics/index2_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
plot(density(model$sims.list$S2), xlab='Index 2', col="black", lwd=2, main=""); 
abline(v = model$q2.5$S2, lty="dashed", lwd=2); abline(v = model$q97.5$S2, lty="dashed", lwd=2); abline(v = model$q50$S2, lwd=2)
#abline(v = model.so$q50$S2, col="red")
#dev.off()

##plot frequency distrib of male repro success
TSires<-model$sims.list$TSires
Nreps<-nrow(TSires)
xi <- 100 #seq(0,max(TSires), 1)
max.pat<-20
Tsiresdist<-matrix(nrow=Nreps, ncol=xi)
for (j in 1:Nreps) {
  Tsiresdist[j,] <-density(TSires[j,], n=xi, from=0, to=max.pat)$y
  x.temp<-density(TSires[j,], n=xi, from=0, to=max.pat)$x
}
Tsiresdist_m<-apply(Tsiresdist, MARGIN = 2, FUN = mean)
Tsiresdist_CI<-apply(Tsiresdist, MARGIN = 2, FUN = quantile, probs = c(0.05, 0.95))
#xi by Tsiresdist_m

#tiff("Figures_genetics/siring_freq_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
plot(x = x.temp, y = Tsiresdist_m, type="l", ylim=c(0.009, max(Tsiresdist_CI)), xlab='Surviving female pups sired', ylab='Proportion of Males', xlim=c(0.7,19.27))
polygon(x = c(x.temp, rev(x.temp)), y = c(Tsiresdist_CI[1,], rev(Tsiresdist_CI[2,])), col="light grey", border = NA)
lines(x = x.temp, y = Tsiresdist_m, lwd=2)
#lines(x = x.temp.so, y = Tsiresdist_m.so, lwd=2, col="red")
#dev.off()

model.es<-model
x.temp.es<-x.temp
Tsiresdist_m.es<-Tsiresdist_m
Tsiresdist_CI.es<-Tsiresdist_CI

##PREP HARBOR SEAL DATA
hs.dat<-read.csv(str_c(folder, 'harbor_seal_data.csv'))

nmales<-length(unique(hs.dat$male))
Nyrs<-length(unique(hs.dat$year))
Sires<-c(sum(hs.dat$paternities[which(hs.dat$year==1998)]), sum(hs.dat$paternities[which(hs.dat$year==1999)]), sum(hs.dat$paternities[which(hs.dat$year==2000)])) ##number of pups assigned per year
#sampled 70 males and 136 offspring
#sampled c(38, 52, 46) offspring across years
#assigned c(3,18,6) pups across years (male and female)
proppup<-c(0.3, 0.3, 0.3) #proportion of pups sampled (average of all years)
propdad<-c(0.5, 0.5, 0.5) #proportion of adult males sampled (average of all years)
Npupsamp<-c(38, 52, 46) #number of pups sampled
Npuptot<-round(c(38/.3, 52/.3, 46/.3),0) ##number of pups in population
Nm<-c(length(unique(hs.dat$male[which(hs.dat$year==1998)])), length(unique(hs.dat$male[which(hs.dat$year==1999)])), length(unique(hs.dat$male[which(hs.dat$year==2000)])))

obs<-matrix(nrow=length(unique(hs.dat$male)), ncol=length(unique(hs.dat$year)))
years<-unique(hs.dat$year)
for (m in 1:nrow(obs)) {
  for (y in 1:length(years)) {
    if (length(which(hs.dat$male==hs.dat$male[m] & hs.dat$year==years[y]))==0) {
      obs[m,y]<--1
    } else {
      obs[m,y]<-sum(hs.dat$paternities[which(hs.dat$male==hs.dat$male[m]& hs.dat$year==years[y])])
    }
  }
}

##maleyrs = vector with length number of males (number of years that male was active)
maleyrs<-dim(0)
for (j in 1:nrow(obs)) {
  maleyrs<-c(maleyrs, length(which(obs[j,]!=-1)))
}

##M = matrix of line numbers of males matrix that are active each year; columns are years
M<-matrix(nrow=nrow(obs), ncol=ncol(obs))
Nm<-dim(0) ##number of males active each year
for (j in 1:length(years)) {
  M[1:length(which(obs[,j]!=-1)),j]<-which(obs[,j]!=-1)
  Nm<-c(Nm, length(which(obs[,j]!=-1)))
}

##ESTIMATE REPRO LIFETIME
##ages 3 to 4, 4 to 5, 5 to 6, 6 to 7
##survival of 0.879
age.temp<-c(3, 4, 5, 6, 7, 25, 25, 25) ##age data. must be entered manually. start at first age of repro
start<-1000
n.temp<-c(start*0.879^0, start*0.879, start*0.879^2, start*0.879^3, start*0.879^4,0, 0, 0) ##number of males that survive to that age
#n.temp<-n.temp*1000
plot(age.temp, n.temp)
func.temp<-nls(formula = n.temp~a*age.temp^b, start = list(a=10000,b=-1))
curve(coef(func.temp)[1]*x^coef(func.temp)[2], add=T)
lines(x= age.temp, y=age.temp*coef(lm(n.temp~age.temp))[2]+coef(lm(n.temp~age.temp))[1])

out.temp<-dim(0)
for (j in 5:22) { ##from min age to max age
  #curve version
  n.j<-coef(func.temp)[1]*j^coef(func.temp)[2] ##animals that survive to j
  n.j1<-coef(func.temp)[1]*(j+1)^coef(func.temp)[2] ##animals that survive to j+1
  
  #linear model version
  #n.j<-j*coef(lm(n.temp~age.temp))[2]+coef(lm(n.temp~age.temp))[1]
  #n.j1<-(j+1)*coef(lm(n.temp~age.temp))[2]+coef(lm(n.temp~age.temp))[1]
  
  n.new<-round(n.j-n.j1, digits = 0) ##new animals that died between j and j+1
  span.temp<-rep(j-(5-1),n.new) ##rep repro lifetime for that many males
  #print(span.temp); print(j)
  out.temp<-c(out.temp, span.temp)
}
mean.lifetime<-round(mean(out.temp),0); mean.lifetime ##mean repro lifetime

##OVERWRTIE MALEYRS WITH 3 FOR EVERY MALE; actually the mean repro lifetime
maleyrs<-rep(mean.lifetime, length(maleyrs))

Bayes.in<-list()
Bayes.in$nmales<-nmales
Bayes.in$Nyrs<-Nyrs
#Bayes.in$scalefact<-0.9
Bayes.in$Sires<-Sires
Bayes.in$proppup<-proppup
Bayes.in$propdad<-propdad
#Bayes.in$Npuptot<-Npuptot
Bayes.in$Npupsamp<-Npupsamp
Bayes.in$obs<-obs
Bayes.in$Nm<-Nm
Bayes.in$M<-M
Bayes.in$maleyrs<-maleyrs

##HARBOR SEAL BAYESIAN MODEL
var.names<-c('S1', 'S2', 'mu', 'sig', 'bta', 'TSires')
burnin<-10000
nsamps<-20000
##run at 20,000 and 10,000
nchains<-3
model <- jags(data = Bayes.in, parameters.to.save = var.names, model.file = "Reproskew_JAGS.txt", n.chains = nchains, n.iter = nsamps, n.burnin = burnin)

#plot(model)

##PLOTS
#resolution<-96

#tiff("Figures_genetics/index1_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
#plot(density(model$sims.list$S1), xlab='Index 1', col="black", lwd=2, main=""); 
#abline(v = model$q2.5$S1, lty="dashed", lwd=2); abline(v = model$q97.5$S1, lty="dashed", lwd=2); abline(v = model$q50$S1, lwd=2)
#dev.off()

#tiff("Figures_genetics/index2_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
#plot(density(model$sims.list$S2), xlab='Index 2', col="black", lwd=2, main=""); 
#abline(v = model$q2.5$S2, lty="dashed", lwd=2); abline(v = model$q97.5$S2, lty="dashed", lwd=2); abline(v = model$q50$S2, lwd=2)
#abline(v = model.so$q50$S2, col="red")
#dev.off()

##plot frequency distrib of male repro success
TSires<-model$sims.list$TSires
Nreps<-nrow(TSires)
xi <- 100 #seq(0,max(TSires), 1)
max.pat<-20
Tsiresdist<-matrix(nrow=Nreps, ncol=xi)
for (j in 1:Nreps) {
  Tsiresdist[j,] <-density(TSires[j,], n=xi, from=0, to=max.pat)$y
  x.temp<-density(TSires[j,], n=xi, from=0, to=max.pat)$x
}
Tsiresdist_m<-apply(Tsiresdist, MARGIN = 2, FUN = mean)
Tsiresdist_CI<-apply(Tsiresdist, MARGIN = 2, FUN = quantile, probs = c(0.05, 0.95))

#tiff("Figures_genetics/siring_freq_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 14)
#plot(x = x.temp, y = Tsiresdist_m, type="l", ylim=c(0.009, max(Tsiresdist_CI)), xlab='Surviving female pups sired', ylab='Proportion of Males', xlim=c(0.7,19.27))
#polygon(x = c(x.temp, rev(x.temp)), y = c(Tsiresdist_CI[1,], rev(Tsiresdist_CI[2,])), col="light grey", border = NA)
#lines(x = x.temp, y = Tsiresdist_m, lwd=2)
#lines(x = x.temp.so, y = Tsiresdist_m.so, lwd=2, col="red")
#dev.off()

model.hs<-model
x.temp.hs<-x.temp
Tsiresdist_m.hs<-Tsiresdist_m
Tsiresdist_CI.hs<-Tsiresdist_CI

##FINAL PLOTS WITH ALL SPECIES
resolution<-96

models<-list()
models$Enhydra<-model.so
models$Mirounga<-model.es
models$Phoca<-model.hs

##PLOT INDEX 1
tiff("Figures_genetics/index1_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 18)
##get limits
par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,3,1,1))
xmax<-0
ymax<-0
for (j in 1:length(models)) {
  if (max(density(models[[j]]$sims.list$S1)$y)>ymax) {ymax<-max(density(models[[j]]$sims.list$S1)$y)}
  if (max(density(models[[j]]$sims.list$S1)$x)>xmax) {xmax<-max(density(models[[j]]$sims.list$S1)$x)}
}

color<-c("black","red","blue")
for (j in 1:length(models)) {
  if (j ==1) {
    plot(density(models[[j]]$sims.list$S1), xlab='Standardized variance in LRS', ylab="Density", col=color[j], lwd=2, main="", xlim=c(0, xmax), ylim=c(0.1,ymax)) ##index 1
  } else {
    lines(density(models[[j]]$sims.list$S1), col=color[j], lwd=2)
  }
  #abline(v = models[[j]]$q2.5$S1, lty="dashed", lwd=2, col=color[j])
  #abline(v = models[[j]]$q97.5$S1, lty="dashed", lwd=2, col=color[j])
  abline(v = models[[j]]$q50$S1, lwd=2, col=color[j], lty="dashed")
}
legend(x = xmax-xmax/2, y = ymax-ymax/10, legend = c("Sea otter", "Elephant seal", "Harbor seal"), lty = "solid", col = color, bty = "n")
dev.off()

##PLOT INDEX 2
tiff("Figures_genetics/index2_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 18)
##get limits
par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,3,1,1))
xmax<-0
ymax<-0
for (j in 1:length(models)) {
  if (max(density(models[[j]]$sims.list$S2)$y)>ymax) {ymax<-max(density(models[[j]]$sims.list$S2)$y)}
  if (max(density(models[[j]]$sims.list$S2)$x)>xmax) {xmax<-max(density(models[[j]]$sims.list$S2)$x)}
}

color<-c("black","red","blue")
for (j in 1:length(models)) {
  if (j ==1) {
    plot(density(models[[j]]$sims.list$S2), xlab=expression(italic("S")[3]), col=color[j], lwd=2, main="", xlim=c(0, xmax), ylim=c(0.7,ymax), ylab="Density")
  } else {
    lines(density(models[[j]]$sims.list$S2), col=color[j], lwd=2)
  }
  #abline(v = models[[j]]$q2.5$S1, lty="dashed", lwd=2, col=color[j])
  #abline(v = models[[j]]$q97.5$S1, lty="dashed", lwd=2, col=color[j])
  abline(v = models[[j]]$q50$S2, lwd=2, col=color[j], lty="dashed")
}
legend(x = xmax/50, y = ymax-ymax/10, legend = c("Sea otter", "Elephant seal", "Harbor seal"), lty = "solid", col = color, bty = "n")
dev.off()

x.temp<-list()
x.temp$Enhydra<-x.temp.so
x.temp$Mirounga<-x.temp.es
x.temp$Phoca<-x.temp.hs

Tsiresdist_m<-list()
Tsiresdist_m$so<-Tsiresdist_m.so
Tsiresdist_m$es<-Tsiresdist_m.es
Tsiresdist_m$hs<-Tsiresdist_m.hs

Tsiresdist_CI<-list()
Tsiresdist_CI$so<-Tsiresdist_CI.so
Tsiresdist_CI$es<-Tsiresdist_CI.es
Tsiresdist_CI$hs<-Tsiresdist_CI.hs

##PLOT SIRING FREQUENCY
tiff("Figures_genetics/siring_freq_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 18)
par(mfrow=c(1,1), mgp=c(1.5,0.5,0))
ymax<-0
for (j in 1:length(Tsiresdist_CI)) {
  if (max(Tsiresdist_CI[[j]][2,])>ymax) {ymax<-max(Tsiresdist_CI[[j]][2,])}
}
color<-c("black","red","blue")
trans<-100
poly.color<-c(rgb(0, 0, 0, trans, maxColorValue=255), rgb(255, 0, 0, trans, maxColorValue=255), rgb(0, 0, 255, trans, maxColorValue=255))
for (j in 1:length(models)) {
  if (j==1) {
    plot(x = x.temp[[j]], y = Tsiresdist_m[[j]], type="l", xlab='Pups sired', ylab='Proportion of Males', col=color[j], ylim=c(0, ymax), xlim=c(0.7, 19.3))
  }
  polygon(x = c(x.temp[[j]], rev(x.temp[[j]])), y = c(Tsiresdist_CI[[j]][1,], rev(Tsiresdist_CI[[j]][2,])), col=poly.color[j], border = NA)
  lines(x = x.temp[[j]], y = Tsiresdist_m[[j]], col=color[j], lwd=2)
}
legend(x = 12, y = ymax-ymax/10, legend = names(x.temp), col = color, lty = "solid", bty = "n")

dev.off()

##SIRING FREQUENCY PLOT 2
tiff("Figures_genetics/siring_freq_plot2.tiff", family="Arial", height=9, width=6.5, units="in", compression="lzw", res=resolution, pointsize = 18)
par(mfrow=c(3,1), mgp=c(1.5,0.5,0), mar = c(3,3,1,1))
ymax<-0
for (j in 1:length(Tsiresdist_CI)) {
  if (max(Tsiresdist_CI[[j]][2,])>ymax) {ymax<-max(Tsiresdist_CI[[j]][2,])}
}
color<-c("black","red","blue")
trans<-100
poly.color<-c(rgb(0, 0, 0, trans, maxColorValue=255), rgb(255, 0, 0, trans, maxColorValue=255), rgb(0, 0, 255, trans, maxColorValue=255))
for (j in 1:length(models)) {
  plot(x = x.temp[[j]], y = Tsiresdist_m[[j]], type="l", xlab='', ylab='Proportion of Males', col=color[j], ylim=c(.007, ymax), xlim=c(0.7, 19.3), axes = F)
  axis(side = 1, at = seq(0,20, 1))
  axis(side = 2, at = seq(0,.2,.05))
  if (names(x.temp)[j]=="Enhydra") {
    title(xlab="Surviving female pups sired")
  }
  else {
    title(xlab="Pups sired")
  }
  polygon(x = c(x.temp[[j]], rev(x.temp[[j]])), y = c(Tsiresdist_CI[[j]][1,], rev(Tsiresdist_CI[[j]][2,])), col=poly.color[j], border = NA)
  lines(x = x.temp[[j]], y = Tsiresdist_m[[j]], col=color[j], lwd=2)
  legend(x = 15, y = ymax-ymax/10, legend = names(x.temp)[j], col = color[j], lty = "solid", bty = "n")
}

dev.off()

##PLOTS OF BETA 1 AND BETA 2
tiff("Figures_genetics/beta_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 18)
#par(mfrow=c(2,1), mar = c(3,3,0.25,0.25), mgp=c(1.5,0.5,0))
#plot(density(model.so$sims.list$bta1), main="", xlab=expression(beta[1]), lwd=2, ylab="Density")
#abline(v = model.so$q50$bta1, lwd=2)
#abline(v = model.so$q2.5$bta1, lty="dashed", lwd=2)
#abline(v = model.so$q97.5$bta1, lty="dashed", lwd=2)

#plot(density(model.so$sims.list$bta2), main="", xlab=expression(beta[2]), lwd=2, ylab="Density")
#abline(v = model.so$q50$bta2, lwd=2)
#abline(v = model.so$q2.5$bta2, lty="dashed", lwd=2)
#abline(v = model.so$q97.5$bta2, lty="dashed", lwd=2)

plot(density(model.so$sims.list$bta), main="", xlab=expression(beta), lwd=2, ylab="Density")
abline(v = model.so$q50$bta, lwd=2)
abline(v = model.so$q2.5$bta, lty="dashed", lwd=2)
abline(v = model.so$q97.5$bta, lty="dashed", lwd=2)

#mtext(text = "Probability density", outer=T, line = 0, side = 2)
dev.off()
par(mfrow=c(1,1))

##PLOT SIRING FREQU FOR SEA OTTERS
tiff("Figures_genetics/enhydra_siring_freq_plot.tiff", family="Arial", height=6.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 18)
par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,3,1,1))
ymax<-0
for (j in 1) {
  if (max(Tsiresdist_CI[[j]][2,])>ymax) {ymax<-max(Tsiresdist_CI[[j]][2,])}
}

trans<-100
poly.color<-rgb(0, 0, 0, trans, maxColorValue=255)
color<-c("black","red","blue")
j<-1
plot(x = x.temp[[j]], y = Tsiresdist_m[[j]], type="l", xlab='Surviving female pups', ylab='Proportion of Males', col=color[j], ylim=c(0.007, ymax), xlim=c(0.7, 19.4))
polygon(x = c(x.temp[[j]], rev(x.temp[[j]])), y = c(Tsiresdist_CI[[j]][1,], rev(Tsiresdist_CI[[j]][2,])), col=poly.color[j], border = NA)
lines(x = x.temp[[j]], y = Tsiresdist_m[[j]], col=color[j], lwd=2)

dev.off()

##PLOT SKEW INDICES AS MEDIAN WITH CI
skew.summary<-data.frame(species=c("E. lutris", "M. leonina", "P. vitulina"), index= c(rep("S1", 3), rep("S2", 3)), q50=c(model.so$q50$S1, model.es$q50$S1, model.hs$q50$S1, model.so$q50$S2, model.es$q50$S2, model.hs$q50$S2), q2.5=c(model.so$q2.5$S1, model.es$q2.5$S1, model.hs$q2.5$S1, model.so$q2.5$S2, model.es$q2.5$S2, model.hs$q2.5$S2), q97.5=c(model.so$q97.5$S1, model.es$q97.5$S1, model.hs$q97.5$S1, model.so$q97.5$S2, model.es$q97.5$S2, model.hs$q97.5$S2))

plot1<-ggplot(data=subset(skew.summary, index=="S1"), aes(species, q50)) + geom_point(size=4) + geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=0.1) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.x=element_text(vjust=0.4, size=14), axis.title.y=element_text(vjust=0.4, size=12)) +xlab("") +ylab("Standardized variance in LRS")
plot1 <- plot1 + theme(axis.title.y = element_text(vjust=1)) ##move y axis title away from axis
plot1 <- plot1 + scale_y_continuous(breaks=seq(0,15,2)) ##add more breaks in y axis
plot1 <- plot1 + theme(axis.text.x = element_text(angle = 45, vjust = 1, face = "italic", hjust=1)) ##rotate x axis labels
plot1

plot2<-ggplot(data=subset(skew.summary, index=="S2"), aes(species, q50)) + geom_point(size=4) + geom_errorbar(aes(ymin=q2.5, ymax=q97.5), width=0.1) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.x=element_text(vjust=0.4, size=14), axis.title.y=element_text(vjust=0.4, size=12)) +xlab("") 
plot2 <- plot2 + ylab(expression(italic(S[3])~ "of repro. variance"))
plot2 <- plot2 + theme(axis.title.y = element_text(vjust=1.1)) ##move y axis title away from axis
plot2 <- plot2 + scale_y_continuous(breaks=seq(0,1, 0.1), limits=c(0,1)) ##add more breaks in y axis
plot2 <- plot2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, face = "italic", hjust = 1)) ##rotate x axis labels
plot2

# Define grid layout to locate plots and print each graph
grid.arrange(plot1, plot2, ncol=2)

tiff("Figures_genetics/variance_indices_plot.tiff", family="Arial", height=3.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 12)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

##SUMMARY STATS
print('sea otter: lifetime repro success'); quantile(model.so$sims.list$TSires, c(0.25, 0.5, 0.75, 0.90))

print('elephant seal: lifetime repro success'); quantile(model.es$sims.list$TSires, c(0.25, 0.5, 0.75, 0.9))

print('harbor seal: lifetime repro success'); quantile(model.hs$sims.list$TSires, c(0.25, 0.5, 0.75, 0.9))

print('sea otter index 1'); model.so$q50$S1; model.so$q2.5$S1; model.so$q97.5$S1
print('sea otter index 2'); model.so$q50$S2; model.so$q2.5$S2; model.so$q97.5$S2
print('sea otter index 2 mean'); model.so$mean$S2; model.so$sd$S2

print('e seal index 1'); model.es$q50$S1; model.es$q2.5$S1; model.es$q97.5$S1
print('e seal index 2'); model.es$q50$S2; model.es$q2.5$S2; model.es$q97.5$S2

print('harbor seal index 1'); model.hs$q50$S1; model.hs$q2.5$S1; model.hs$q97.5$S1
print('harbor seal index 2'); model.hs$q50$S2; model.hs$q2.5$S2; model.hs$q97.5$S2

#print('beta 1 mean sd'); mean(model.so$sims.list$bta1); sd(model.so$sims.list$bta1)
#print('beta 2 mean sd'); mean(model.so$sims.list$bta2); sd(model.so$sims.list$bta2)

print('beta mean sd'); mean(model.so$sims.list$bta); sd(model.so$sims.list$bta)

##TEST EFFECT OF DIFFERENT POPULATION SIZE
##AUTOMATED version with multiple percent changes
per.change<-c(-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50)

Bayes.in<-list()
files<-str_c(folder, 'Bayes_in/', list.files(str_c(folder,'Bayes_in/')))
for (j in 1:length(files)) {
  file.temp<-as.matrix(read.csv(files[j]))
  if (ncol(file.temp)==1) {
    file.temp<-as.vector(file.temp)
  }
  #file.name.temp<-str_sub(files[j], 10, -5)
  file.name.temp<-str_split(files[j], "/")[[1]]
  file.name.temp<-file.name.temp[length(file.name.temp)] %>% str_sub(start = 0, end=-5)
  Bayes.in[[file.name.temp]]<-file.temp
}
names(Bayes.in)
propdad.original<-Bayes.in$propdad

#out<-list()
out<-dim(0)
for (j in 1:length(per.change)) {
  Bayes.in$propdad<-propdad.original*per.change[j]/100+propdad.original ##assign percent change to proportion of dads sampled
  
  var.names<-c('S1', 'S2')
  burnin<-10000
  nsamps<-20000
  ##run at 20,000 and 10,000
  nchains<-3
  model <- jags(data = Bayes.in, parameters.to.save = var.names, model.file = "Reproskew_JAGS.txt", n.chains = nchains, n.iter = nsamps, n.burnin = burnin)
  
  ##NEED TO SAVE EACH VERSION OF THE MODEL
  out<-rbind(out, c(per.change[j], mean(Bayes.in$propdad), model$q50$S1, model$q2.5$S1, model$q97.5$S1, model$q50$S2, model$q2.5$S2, model$q97.5$S2))
  #out[[j]]<-model
}
out<-data.frame(out); colnames(out)<- c("per.change", "mean.prop", "S1q50", "S1q2.5", "S1q97.5", "S2q50", "S2q2.5", "S2q97.5")

##plot estimate in variance as a function of percent change
plot1 <- ggplot(data=out, aes(mean.prop, S1q50)) 
plot1 <- plot1 + geom_point(size=4) + geom_errorbar(aes(ymin=S1q2.5, ymax=S1q97.5), width=0.01) 
plot1 <- plot1 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.x=element_text(vjust=0.4, size=12), axis.title.y=element_text(vjust=0.4, size=12)) 
plot1 <- plot1 + xlab("Estimated prop. of males sampled") +ylab("Standardized variance in LRS")
plot1 <- plot1 + theme(axis.title.y = element_text(vjust=1)) ##move y axis title away from axis
plot1 <- plot1 + scale_y_continuous(breaks=seq(0,20, 1)) ##add more breaks in y axis
plot1 <- plot1 + theme(axis.text.x = element_text(vjust = 1)) ##rotate x axis labels
plot1 <- plot1 + scale_x_continuous(breaks = seq(0,1,0.1))
plot1

plot2 <- ggplot(data=out, aes(mean.prop, S2q50)) 
plot2 <- plot2 + geom_point(size=4) + geom_errorbar(aes(ymin=S2q2.5, ymax=S2q97.5), width=0.01) 
plot2 <- plot2 + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), panel.background=element_blank(), axis.line=element_line(colour="black"), axis.ticks=element_line(colour="black"), axis.text=element_text(colour="black"), axis.text.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.x=element_text(vjust=0.4, size=12), axis.title.y=element_text(vjust=0.4, size=12)) 
plot2 <- plot2 + xlab("Estimated prop. of males sampled") + ylab(expression(italic(S[3])~ "of repro. variance"))
plot2 <- plot2 + theme(axis.title.y = element_text(vjust=1)) ##move y axis title away from axis
plot2 <- plot2 + scale_y_continuous(breaks=seq(0,1, 0.1)) ##add more breaks in y axis
plot2 <- plot2 + theme(axis.text.x = element_text(vjust = 1)) ##rotate x axis labels
plot2 <- plot2 + scale_x_continuous(breaks = seq(0,1,0.1))
plot2

# Define grid layout to locate plots and print each graph
grid.arrange(plot1, plot2, ncol=2)


##SAVE PLOT
tiff("Figures_genetics/sample_size_indices_plot.tiff", family="Arial", height=3.5, width=7, units="in", compression="lzw", res=resolution, pointsize = 12)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

