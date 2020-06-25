###############################################################
# JAGS model for fitting functional responses - written by Bart DiFiore
###############################################################

source(here::here("code", "0_libraries.R"))
source(here::here("code","8a_JAGS_functions.R"))

##############################
# import foraging assay data
##############################

foraging_data <- read.csv(here::here("data", "functional_response", "raw", "foraging_assay_data.csv"))

##############################
# some wrangling
##############################

# give each lobster a unique ID number
foraging_data$id <- as.numeric(foraging_data$lobster_id)

# add "t" in front of each temperature and coerce to factor
foraging_data$temp2 <- as.factor(paste("t", foraging_data$temp, sep = ""))

# make new df with only necessary variables
s <- foraging_data[,c( "temp2", "lobster_id","Initial", "Killed")]
# rename them
names(s) <- c("temp", "id", "initial", "killed")

# arrange observations by temp, id, then initial
s <- arrange(s, temp, id, initial)

##############################
# plot raw foraging assay data for each lobster
##############################

# plot raw data where killed is a function of initial; facet by id
lattice:: xyplot(killed~initial|as.factor(id), data = s)

##############################
# define heirarchical model
##############################

jagsscript = cat("
                 
                 model{
                 
                 
                 # PRIORS
                 
                 # Population level estimates (hyperprior)
                 
                 mu.logit.a ~ dnorm(0, 0.01)
                 mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))
                 
                 mu.log.h ~ dnorm(0, 0.01)
                 mu.h <- exp(max(min(mu.log.h, 20), -20))
               
            
                 # Treatment level effects
                 for(i in 1:Ntreats){
                 #a
                 t.logit.a[i] ~ dnorm(mu.logit.a, t.tau.a) # tau is a measure of variance (going to vary for each interaction of the model)
                 t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))
                 
                 #h
                 t.log.h[i] ~ dnorm(mu.log.h, t.tau.h)
                 t.h[i] <- exp(max(min(t.log.h[i],20),-20))
                 
                 }
                 
                 
                 # Individual level effects
                 for(i in 1:num.ind){
                 
                 # a (make logit.a)
                 loga[i] ~ dnorm(t.logit.a[tind[i]], tau_int.a) # priors come from normal dist
                 a[i] <- exp(loga[i])/(1+exp(loga[i]))
                 
                 # h (make logit.h)
                 logh[i] ~ dnorm(t.log.h[tind[i]], tau_int.h)
                 h[i] <- exp(max(min(logh[i],10),-20))
                 }
                 
                 # functional response likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
                 
                 killed[i] ~ dbin(prob[i],initial[i]) # probability of being killed given any combo of a & h
                 
                 }
                 
                 
                 
                 # Variances for all levels 
                 sigma_int.a ~dunif(0,10) # variance drawn from 0-10
                 tau_int.a <- 1/(sigma_int.a*sigma_int.a) # transform var into precision and feed back into precision estimate for treatment level 
                 sigma_int.h ~dunif(0,10)
                 tau_int.h <- 1/(sigma_int.h*sigma_int.h)
                 
                 sigma_t.a ~ dunif(0,10)
                 t.tau.a <- 1/(sigma_t.a*sigma_t.a)
                 sigma_t.h ~ dunif(0,10)
                 t.tau.h <- 1/(sigma_t.h*sigma_t.h)
                 
                 }
                 
                 ", file = here::here("code", "heirarchical_jagsSAM.txt"))

##############################
# ?
##############################

model.loc=here::here("code", "heirarchical_jagsSAM.txt")
jags.params=c("a", "h", "t.a", "t.h", "mu.a", "mu.h") # can also include sigmas here to track

##############################
# set factor levels to ensure curves match up with the right data
##############################

# find distinct lobster
tind <- distinct(s, temp, id) %>%
  arrange(id)

temps <- sort(unique(foraging_data$temp))
unique(as.numeric(tind$temp))

# set factor levels of individuals by treatment
tind <- as.factor(as.vector(tind$temp))

# Test against simulated dataset.
jags.data = list("initial"= s$initial,
                 "killed" = s$killed,
                 "id" = s$id,
                 "num.ind" = length(unique(s$id)),
                 "n" = length(s$initial), 
                 "t" = s$temp, 
                 "Ntreats" = length(unique(s$temp)), 
                 "tind" = tind, 
                 "T" = 24) # named list

##############################
# set number of model iterations
##############################

n.chains = 3 # apparently pretty standard
n.burnin = 250000 # how many estimates to throw out; # 250,000; 10,000
n.thin = 10 # only keep every other estimate; # 10; 2
n.iter = 500000 # number of interations ; # 500,000; 25,000

##############################
# run model 
##############################

model = R2jags::jags(jags.data, parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

##############################
# run model 
##############################

# NOTE if Rhat < 1.2 then chains have
a = MCMCsummary(model,params=c('a'), round = 4) 
h = MCMCsummary(model,params=c('h'), round = 4)

# treatment-level estimates
t.a = MCMCsummary(model,params=c('t.a'), round = 4)
t.h = MCMCsummary(model,params=c('t.h'), round = 4)

# order associated with treatments: 21, 11, 26, 16 
# unique(as.numeric(foraging_data$temp))

# population-level estimates
mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')

# lobster_ids
ids <- sort(unique(foraging_data$lobster_id))
# numeric factor levels associated with each lobster--I matched lobster ID with numeric idenitifier and then confirmed
unique(as.numeric(s$id))

# temperature ids
t.ids <- sort(unique(foraging_data$temp))
# numeric factor levels associated with each lobster--I matched lobster ID with numeric idenitifier and then confirmed
unique(as.numeric(s$temp))

##############################
# plot curves and their associated data at the individual level
##############################

# png("figures/6_individual_fits_JAGS_T=24_NEW.png", width = 2000, height = 2000, res = 150)
d <- par(mfrow = c(5,5), mar = c(4,4,1,1))
for(i in 1:22){
  plot(I(killed/24) ~ jitter(initial),data = s[as.numeric(s$id) == i, ], xlab="Number of prey offered",ylab="Consumption rate", main = paste(ids[i]), ylim = c(0,1.6))
  curve(holling2(x,a[i,1],h[i,1],P=1,T=1),add=TRUE,lty=1) #true curve
  text(x = 10, y = 60, label = paste("a =", round(a[i,1], 3), "\n", "h = ", round(h[i,1], 3)))
  abline(a = 0,b = 1, lty = 4)
}
par(d)
# dev.off()

##############################
# plot curves and their associated data at the treatment level
##############################

# png("figures/6_treatment_fits_JAGS_T=24_NEW.png", width = 2000, height = 2000, res = 150)
d <- par(mfrow = c(2,2), mar = c(4,4,1,1))
for(i in 1:4){
  plot(I(killed/24) ~ jitter(initial),data = s[as.numeric(s$temp) == i, ], xlab="Number of prey offered",ylab="Consumption rate", ylim = c(0,1.6), main = paste(levels(s$temp)[i]))
  curve(holling2(x,t.a[i,1],t.h[i,1],P=1,T=1),add=TRUE,lty=1) #true curve
  # text(x = 10, y = 60, label = paste("a =", round(a[i,1], 3), "\n", "h = ", round(h[i,1], 3)))
  # abline(a = 0,b = 1, lty = 4)
}
par(d)
# dev.off()

##############################
# visualize mixing of chains
##############################

MCMCplot(model, params = c("a"), rank= T)
MCMCplot(model, params = c("h"), rank= T)

MCMCtrace(model)