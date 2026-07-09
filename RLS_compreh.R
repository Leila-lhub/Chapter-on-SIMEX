
#RLS - 

#setwd("C:/Users/nombo/OneDrive/Documents/~/Code R pour article 1")
source("Subroutines_comprendre.R")

#Simulation part
Epsilons = c(1,2,5,10,15,20,30,40,60)
for(val in Epsilons){
  set.seed(1836)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2), perturbY = TRUE, epsilon =val , n = 500, k=20, nsim = 500))
  #assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2), perturbY = TRUE, epsilon =val , n = 500, k=20, nsim = 5))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLS/RLS_sim",val,".Rdata", sep = ""))
  #save(list = paste("out",val,sep =""), file = paste("FichiersRLS/RLS_tests",val,".Rdata", sep = ""))
  print(val)
}

#I tried it again with RLS2 as output instead of RLS
#The difference is that B=500 instead of 200
#Didn't seem to change anything, so we'll keep B=200; it's much faster...

#Extrapolation part
for(val in Epsilons){
  load(paste("./FichiersRLS/FichiersRLS/RLS2_sim",val,".Rdata", sep = ""))
  #load(paste("FichiersRLS/RLS_tests",val,".Rdata", sep = ""))
  assign(paste("extrap",val, sep = ""),sim_extrap(get(paste("out",val,sep = "")),ps=seq(0,2,length.out = 20)))
  print(val)
}


#Pour faire les graphiques facilement, mettre le tout dans une grande matrice
nsim=500
Epsilons = c(1,2,5,10,15,20,30,40,60)
results = data.frame(matrix(NA, nrow = length(Epsilons)*4*nsim, ncol = 4))
names(results) = c("Epsilon", "Method", "Iter", "Estimate")
results$Epsilon = rep(Epsilons, each = 4*nsim)
results$Method = rep(c("true", "noisy", "linear", "quadratic"), each = nsim)
results$Iter = rep(1:nsim, length.out = length(Epsilons)*4*nsim)

#Ajouter ici les estimés
k=0
for(val in Epsilons){
  results$Estimate[(k+1):(k+(4*nsim))] = c(t(get(paste("extrap",val,sep = ""))))
  print(val)
  k = k + (4*nsim)
}

results$Method = factor(results$Method, levels = c("true", "noisy", "linear", "quadratic", "non-linear"))


#Figures - 

library(ggplot2)
head(results)
p <- ggplot(results, aes(x=as.factor(Epsilon), y=Estimate, fill = Method)) + 
  geom_boxplot() + labs(x="Epsilon", y = "Estimate")
p

#Now, let's compute the averages 

library(dplyr)
mean_data <- group_by(results, Epsilon, Method) %>%
  summarise(MeanEst = mean(Estimate))

sum(results$Estimate == 999) #0; no extrapolation issues without non-linear

ggplot(mean_data, aes(x = Epsilon, y = MeanEst, colour = Method)) +
  geom_point() + geom_line()

#Now, look at the relative bias - by comparing to 0.293

mean_data$RelBias = (mean_data$MeanEst - 0.293)/0.293*100
ggplot(mean_data, aes(x = Epsilon, y = RelBias, colour = Method)) +
  geom_point() + geom_line()

#Do it more properly - by comparing to true value at each iteration

trueVals = c()
k=0
nsim=500
for(i in 1:length(Epsilons)){
  trueVals = c(trueVals, rep(results$Estimate[(k+1):(k+nsim)],4))
  k = k + nsim*4
}

results$RelBiais2 = (results$Estimate - trueVals)/trueVals*100
mean_dataRelBiais2 <- group_by(results, Epsilon, Method) %>%
  summarise(MeanRelBiais = mean(RelBiais2))

p = ggplot(
  subset(mean_dataRelBiais2, Method != "true"), aes(x = Epsilon, y = MeanRelBiais, colour = Method)) +
  geom_point(size = 2) + 
  geom_line(aes(lty = Method), lwd = 0.8) +
  geom_hline(yintercept = 0) +
  ylab("Mean relative bias of slope estimate (in %)") +
  xlab(expression(paste("Privacy parameter ", epsilon))) + 
  theme(text = element_text(size = 15), legend.position = c(0.8, 0.4))

#Compliqué de modifier la légende...
#On utilise un petit tour de passe-passe...

results2 = results
levels(results2$Method) = c("true", "Naive estimate from noisy data", "SIMEX with linear extrapolation", "SIMEX with quadratic extrapolation", "non-linear")

mean_dataRelBiais3 <- group_by(results2, Epsilon, Method) %>%
  summarise(MeanRelBiais = mean(RelBiais2))

ggplot(subset(mean_dataRelBiais3, Method != "true"), aes(x = Epsilon, y = MeanRelBiais, colour = Method)) +   geom_point(size = 2) + 
  geom_line(aes(lty = Method), lwd = 0.8) + geom_hline(yintercept = 0) + 
  ylab("Mean relative bias of slope estimate (in %)") + xlab(expression(paste("Privacy parameter ", epsilon))) + 
  theme(text = element_text(size = 15), legend.position = c(0.7, 0.4), legend.title = element_text(size = 14), legend.text = element_text(size = 14), axis.text = element_text(size = 14))

#Save as figure 1

#Publish a table in appendix, with the variance of each estimates, and the non-linear method incorporated

#Also, compute theoretical values for naive estimates and comment on that in paper

##Boxplot graphe
results$Methode=results$Method
levels(results$Methode) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")

p <- ggplot(results, aes(x=as.factor(Epsilon), y=Estimate, fill = Methode)) + 
  geom_boxplot() + labs(x=expression(paste("Paramètre de confidentialité différentielle ", epsilon)), y = "Estimés de la pente",fill="Méthode:")+
  theme(text = element_text(size = 13), legend.position = c(0.8,0.2), legend.title = element_text(size = 12), legend.text = element_text(size = 12), axis.text = element_text(size = 10))

p


