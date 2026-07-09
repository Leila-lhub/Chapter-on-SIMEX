
#RLS - 

#setwd("C:/Users/nombo/OneDrive/Documents/~/Code R pour article 1")
source("Subroutines-c-gaussien.R")

#Simulation part
Epsilons = c(0.99,0.9,0.8,0.7,0.6,0.5)
Deltas=c(0.001,0.01,0.05,0.1)
for(val in Epsilons){
  for(valdelta in Deltas){
  set.seed(1836)
  assign(paste("out",val,"_",valdelta,sep=""), my.simulation.sim(betas = c(2), perturbY = TRUE, epsilon =val ,delta=valdelta, n = 500, k=20, nsim = 500))
  save(list = paste("out",val,"_",valdelta,sep =""), file = paste("FichiersRLS/RLSGauss_sim",val,"_",valdelta,".Rdata", sep = ""))
  print(valdelta)
  }
}

#I tried it again with RLS2 as output instead of RLS
#The difference is that B=500 instead of 200
#Didn't seem to change anything, so we'll keep B=200; it's much faster...

#Extrapolation part
for(val in Epsilons){
  for(valdelta in Deltas){
  load(paste("FichiersRLS/RLSGauss_sim",val,"_",valdelta,".Rdata", sep = ""))
  #load(paste("FichiersRLS/RLS_tests",val,".Rdata", sep = ""))
  assign(paste("extrap",val,"_",valdelta, sep = ""),sim_extrap(get(paste("out",val,"_",valdelta,sep = "")),ps=seq(0,2,length.out = 20)))
  print(c(val,valdelta))
  }
}


#Pour faire les graphiques facilement, mettre le tout dans une grande matrice
nsim=500
#Epsilons = c(0.99, 0.95, 0.9, 0.5)
#Deltas=c(0.0015,0.001,0.01,0.1)
results = data.frame(matrix(NA, nrow = length(Epsilons)*length(Deltas)*4*nsim, ncol = 5))
names(results) = c("Epsilon","Delta", "Method", "Iter", "Estimate")
results$Epsilon = rep(Epsilons, each = 4*nsim*length(Deltas))
results$Delta = rep(Deltas, each = 4*nsim)
results$Method = rep(c("true", "noisy", "linear", "quadratic"), each = nsim)
results$Iter = rep(1:nsim, length.out = length(Epsilons)*length(Deltas)*4*nsim)



head(results)
#Ajouter ici les estimés
k=0
for(val in Epsilons){
  for(valdelta in Deltas){
  results$Estimate[(k+1):(k+(4*nsim))] = c(t(get(paste("extrap",val,"_",valdelta,sep = ""))))
  print(val)
  k = k + (4*nsim)
  }
}

results$Method = factor(results$Method, levels = c("true", "noisy", "linear", "quadratic", "non-linear"))


results$Methode=results$Method
levels(results$Methode) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")

#Figures - 

library(ggplot2)
head(results)

p <- ggplot(results, aes(x=as.factor(Epsilon), y=Estimate, fill = Methode)) + 
  geom_boxplot() + labs(x="Epsilon", y = "Estimés")+
  theme(text = element_text(size = 13), legend.position = c(0.7, 0.4), legend.title = element_text(size = 14), legend.text = element_text(size = 14), axis.text = element_text(size = 10))

p

ggplot(subset(results, Epsilon= "0.99"), aes(x = as.factor(Delta), y = Estimate, colour = Methode)) +  
  geom_boxplot() + geom_hline(yintercept = 0) + 
  labs(x="Delta", y = "Estimés", legend.position = c(0.8, 0.4))


lab.delta=paste("delta=",Deltas,sep="")
results$Deltas <- factor(results$Delta, labels =lab.delta)

lab.epsilon=paste("Epsilon=",Epsilons,sep="")
results$Epsilons<-factor(results$Epsilon, labels =lab.epsilon)

p2 <- ggplot(subset(results, Method!= "true"), aes(x=as.factor(Epsilon), y=Estimate, fill=Methode)) + 
  geom_boxplot() + geom_hline(yintercept = 0.8) +
  facet_wrap(~Deltas) +  labs(fill="Méthode:")+
  ylab("Biais relatif moyen de l'estimé de la pente (en %)") + xlab(expression(paste("Paramètre de confidentialité différentielle ", epsilon))) + 
  theme(text = element_text(size = 13), legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), axis.text = element_text(size = 10))

p2

#########
p3 <- ggplot(results, aes(x=as.factor(Epsilon), y=Estimate, fill=Methode)) + 
  geom_boxplot() + geom_hline(yintercept = 0.8) +
  facet_wrap(~Delta) +  labs(fill="Méthode:")+
  ylab("Biais relatif moyen de l'estimé de la pente (en %)") + xlab(expression(paste("Paramètre de confidentialité différentielle ", epsilon))) + 
  theme(text = element_text(size = 13), legend.position = "bottom", legend.title = element_text(size = 12), legend.text = element_text(size = 12), axis.text = element_text(size = 10))

p3
############
#Now, let's compute the averages 

library(dplyr)
mean_data <- group_by(results, Epsilon,Delta, Method) %>%
  summarise(MeanEst = mean(Estimate))

sum(results$Estimate == 999) #0; no extrapolation issues without non-linear

ggplot(mean_data, aes(x =Epsilon, y = MeanEst, colour = Method)) +
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
mean_dataRelBiais2 <- group_by(results, Epsilon, Delta, Method) %>%summarise(MeanRelBiais = mean(RelBiais2))

p = ggplot(
  subset(mean_dataRelBiais2, Method != "true"), aes(x = Epsilon, y = MeanRelBiais, colour = Method)) +
  geom_point(size = 2) + 
  geom_line(aes(lty = Method), lwd = 0.8) +
  geom_hline(yintercept = 0) +
  ylab("Mean relative bias of slope estimate (in %)") +
  xlab(expression(paste("Privacy parameter ", epsilon))) + 
  theme(text = element_text(size = 15), legend.position = c(0.8, 0.4))
p
#Compliqué de modifier la légende...
#On utilise un petit tour de passe-passe...
results2 = results
colnames(results2)=c("Epsilon","Delta","Method","Iter","Estimate","Methode","Deltas","Epsilons","RelBiais2")
levels(results2$Methode) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")
head(results2)
levels(results2$Epsilons)
mean_dataRelBiais3 <- group_by(results2, Epsilons, Delta, Methode) %>%
  summarise(MeanRelBiais = mean(RelBiais2))
#head(results); colnames(results)
head(mean_dataRelBiais3)

ggplot(subset(mean_dataRelBiais3, Methode != "Estimateur avec données sans bruit"), aes(x = Delta, y = MeanRelBiais, colour = Methode)) + 
  geom_point(size = 2) + geom_line(aes(lty = Methode), lwd = 0.8) +facet_wrap(~Epsilons)+
  geom_hline(yintercept = 0) + ylim(-100,-95)+
  ylab("Biais relatif moyen de l'estimé de la pente (en %)") + 
  xlab(expression(paste("Paramètre de confidentialité différentielle ", delta))) +
  theme(text = element_text(size = 13), legend.title = element_text(size = 11), axis.text = element_text(size = 10),legend.position = "bottom")

ggplot(subset(mean_dataRelBiais3, Methode != "Estimateur avec données sans bruit"), aes(x = Delta, y = MeanRelBiais, colour = Methode)) + 
  geom_point(size = 2) + geom_line(aes(lty = Methode), lwd = 0.8) +facet_wrap(~Epsilons)+
  geom_hline(yintercept = 0)+
  ylab("Biais relatif moyen de l'estimé de la pente (en %)") + 
  xlab(expression(paste("Paramètre de confidentialité différentielle ", delta))) +
  theme(text = element_text(size = 13), legend.title = element_text(size = 11), axis.text = element_text(size = 10),legend.position = "bottom")

###Graphe finale 
lab.epsilons=paste(expression(epsilon),"=",Epsilons,sep="")

mean_dataRelBiais3 <- group_by(results2, Epsilon, Delta, Methode) %>%
  summarise(MeanRelBiais = mean(RelBiais2))
#head(results); colnames(results)
head(mean_dataRelBiais3)

library(devtools)
library(ggplot2)
library(latex2exp)
levels(factor(mean_dataRelBiais3$Epsilon))=lab.epsilons

appender <- function(string)TeX(paste("$\\epsilon = $", string)) 

gf=ggplot(subset(mean_dataRelBiais3, Methode != "Estimateur avec données sans bruit"), aes(x = Delta, y = MeanRelBiais, colour = Methode)) + 
  geom_line(aes(lty = Methode), lwd = 0.8) + geom_point(size = 2) +facet_wrap(~Epsilon, scales = "fixed", labeller = as_labeller(appender, default = label_parsed))+
  geom_hline(yintercept = 0) + ylim(-100,-95)+
  ylab("Biais relatif moyen de l'estimé de la pente (en %)") + 
  xlab(expression(paste("Paramètre de confidentialité différentielle ", delta))) +
  theme(text = element_text(size = 13), legend.title = element_text(size = 11), axis.text = element_text(size = 10),legend.position = "bottom")

library(scales)
gf+guides(fill=guide_legend(nrow=3,byrow=TRUE))

+scale_x_discrete(labels=c("0.001" = "0.001", "0.01" = "0.01","0.05"="0.05", "0.1"="0.1"))

levels(factor(mean_dataRelBiais3$Delta))
