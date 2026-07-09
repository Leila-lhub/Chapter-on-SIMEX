
#RLM - 2 variables

#Figure 1 - Want to show the impact of the correlation between the variables
#Consider a fixed value of epsilon (20) (tried 10 too, but poor results...)
#Test various correlations
#Keep the same coefficient for both variables

#setwd("C:/Users/nombo/OneDrive/Documents/~/Code R pour article 1")
source("Subroutines_comprendre.R")

Rhos = c(0,0.1,0.25,0.4,0.5,0.6,0.75,0.9)

##La simulation du modele y=beta0+2x1+2x2+rnorm(0,1); on etudie l'influence du niveau de correlation entre x1 et x2
##epsilon=10;  (y, x1 et x2 respectent la conf. diff.)
for(val in  Rhos){
  set.seed(26597)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2,2), perturbY = TRUE, epsilon = 10 , SigmaX = matrix(c(1,val,val,1),nrow = 2), n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_2var_10_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}

#Try again with variables with different effects
##La simulation du modele y=beta0+2x1-2x2+rnorm(0,1); on etudie l'influence du niveau de correlation entre x1 et x2
##epsilon=10; (y, x1 et x2 respectent la conf. diff.)
##beta0=1
for(val in  Rhos){
  set.seed(26597)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2,-2), perturbY = TRUE, epsilon = 10 , SigmaX = matrix(c(1,val,val,1),nrow = 2), n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_2var_10b_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}


#Another attempt
##La simulation du modele y=beta0+2x1+4x2+rnorm(0,1); on etudie l'influence du niveau de correlation entre x1 et x2
##epsilon=10; (y, x1 et x2 respectent la conf. diff.)
for(val in  Rhos){
  set.seed(26597)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2,4), perturbY = TRUE, epsilon = 10 , SigmaX = matrix(c(1,val,val,1),nrow = 2), n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_2var_10c_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}

#Yet another, this time increase epsilon
##La simulation du modele y=beta0+2x1+2x2+rnorm(0,1); on etudie l'influence du niveau de correlation entre x1 et x2
##epsilon=20; (y, x1 et x2 respectent la conf. diff.)
for(val in  Rhos){
  set.seed(26597)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2,2), perturbY = TRUE, epsilon = 20 , SigmaX = matrix(c(1,val,val,1),nrow = 2), n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_2var_20_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}

#Try again with variables with different effects
##La simulation du modele y=beta0+2x1-2x2+rnorm(0,1); on etudie l'influence du niveau de correlation entre x1 et x2
##epsilon=20; (y, x1 et x2 respectent la conf. diff.)
for(val in  Rhos){
  set.seed(26597)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2,-2), perturbY = TRUE, epsilon = 20 , SigmaX = matrix(c(1,val,val,1),nrow = 2), n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_2var_20b_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}


#Another attempt
##La simulation du modele y=beta0+2x1+4x2+rnorm(0,1); on etudie l'influence du niveau de correlation entre x1 et x2
##epsilon=20; (y, x1 et x2 respectent la conf. diff.)
for(val in  Rhos){
  set.seed(26597)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = c(2,4), perturbY = TRUE, epsilon = 20 , SigmaX = matrix(c(1,val,val,1),nrow = 2), n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_2var_20c_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}



#OK, making the real figure

#First setup - 20 
for(val in Rhos){
  load(paste("FichiersRLM/RLM_sim_2var_20_",gsub("[.]","",val),".Rdata", sep = ""))
  assign(paste("extrap",val, sep = ""),sim_extrap(get(paste("out",val,sep = "")),ps=seq(0,2,length.out = 20)))
  print(val)
}

nsim=500
results = data.frame(matrix(NA, nrow =length(Rhos)*4*nsim*2, ncol = 5))
names(results) = c("Rho", "Method", "Var", "Iter", "Estimate")
results$Rho = rep(Rhos, each = 4*500*2)
results$Method = rep(c("true", "noisy", "linear", "quadratic"), each = 500*2)
results$Var = rep( c("X1","X2"), each = 500)
results$Iter = rep(1:500, length.out = length(Rhos)*4*500*2)

k=0
for(val in Rhos){
  results$Estimate[(k+1):(k+4000)] = c(t(get(paste("extrap",val,sep = ""))))
  print(val)
  k = k + 4000
}
results$Method = factor(results$Method, levels = c("true", "noisy", "linear", "quadratic"))
#levels(results$Method) = c("true", "Naive estimate from noisy data","SIMEX - linear extrapolation", "SIMEX - quadratic extrapolation")
levels(results$Method) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")
head(results)
colnames(results)=c("Rho", "Methode", "Var", "Iter", "Estimate")

results1 = results

#Second setup - 20b
for(val in Rhos){
  load(paste("FichiersRLM/RLM_sim_2var_20b_",gsub("[.]","",val),".Rdata", sep = ""))
  assign(paste("extrap",val, sep = ""),sim_extrap(get(paste("out",val,sep = "")),ps=seq(0,2,length.out = 20)))
  print(val)
}

results = data.frame(matrix(NA, nrow = 8*4*500*2, ncol = 5))
names(results) = c("Rho", "Method", "Var", "Iter", "Estimate")
results$Rho = rep(Rhos, each = 4*500*2)
results$Method = rep(c("true", "noisy", "linear", "quadratic"), each = 500*2)
results$Var = rep( c("X1","X2"), each = 500)
results$Iter = rep(1:500, length.out = 8*4*500*2)

k=0
for(val in Rhos){
  results$Estimate[(k+1):(k+4000)] = c(t(get(paste("extrap",val,sep = ""))))
  print(val)
  k = k + 4000
}
results$Method = factor(results$Method, levels = c("true", "noisy", "linear", "quadratic"))
#levels(results$Method) = c("true", "Naive estimate from noisy data","SIMEX - linear extrapolation", "SIMEX - quadratic extrapolation")
levels(results$Method) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")
head(results)
colnames(results)=c("Rho", "Methode", "Var", "Iter", "Estimate")

results2 = results

#Third setup - 20c
for(val in Rhos){
  load(paste("FichiersRLM/RLM_sim_2var_20c_",gsub("[.]","",val),".Rdata", sep = ""))
  assign(paste("extrap",val, sep = ""),sim_extrap(get(paste("out",val,sep = "")),ps=seq(0,2,length.out = 20)))
  print(val)
}

results = data.frame(matrix(NA, nrow = 8*4*500*2, ncol = 5))
names(results) = c("Rho", "Method", "Var", "Iter", "Estimate")
results$Rho = rep(Rhos, each = 4*500*2)
results$Method = rep(c("true", "noisy", "linear", "quadratic"), each = 500*2)
results$Var = rep( c("X1","X2"), each = 500)
results$Iter = rep(1:500, length.out = 8*4*500*2)

k=0
for(val in Rhos){
  results$Estimate[(k+1):(k+4000)] = c(t(get(paste("extrap",val,sep = ""))))
  print(val)
  k = k + 4000
}
results$Method = factor(results$Method, levels = c("true", "noisy", "linear", "quadratic"))
#levels(results$Method) = c("true", "Naive estimate from noisy data","SIMEX - linear extrapolation", "SIMEX - quadratic extrapolation")
levels(results$Method) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")
colnames(results)=c("Rho", "Methode", "Var", "Iter", "Estimate")

results3 = results



#Figure

library(ggplot2)


library(dplyr)
mean_data1 <- group_by(results1, Rho, Methode, Var) %>%
  summarise(MeanEst = mean(Estimate))
mean_data2 <- group_by(results2, Rho, Methode, Var) %>%
  summarise(MeanEst = mean(Estimate))
mean_data3 <- group_by(results3, Rho, Methode, Var) %>%
  summarise(MeanEst = mean(Estimate))


ggplot(mean_data1, aes(x = as.numeric(Rho), y = MeanEst, colour = Methode)) +
  geom_point() + geom_line() +facet_wrap(~Var)

#C'est un peu surprenant, mais les vrais coefficients varient avec la valeur du coefficient de corrélation entre les prédicteurs. J'ai vérifié séparément dans un autre code, et c'est bien le cas. C'est sûrement à cause de la normalization de Y. 

#Now, look at the relative bias

trueVals = c()
k=0
for(i in 1:8){
  trueVals = c(trueVals, rep(results1$Estimate[(k+1):(k+1000)],4))
  k = k + 4000
}
results1$RelBiais = (results1$Estimate - trueVals)/trueVals*100
results1$Var = as.factor(results1$Var)
#levels(results1$Var) = c("Slope for X1", "Slope for X2")
levels(results1$Var) = c("Coefficient estimé pour X1", "Coefficient estimé pour X2")
mean_data1 <- group_by(results1, Rho, Methode, Var) %>%
  summarise(MeanRelBiais = mean(RelBiais))
head(mean_data1)

trueVals = c()
k=0
for(i in 1:8){
  trueVals = c(trueVals, rep(results2$Estimate[(k+1):(k+1000)],4))
  k = k + 4000
}
results2$RelBiais = (results2$Estimate - trueVals)/trueVals*100
results2$Var = as.factor(results2$Var)
#levels(results2$Var) = c("Slope for X1", "Slope for X2")
levels(results2$Var) =c("Coefficient estimé pour X1", "Coefficient estimé pour X2")
mean_data2 <- group_by(results2, Rho, Methode, Var) %>%
  summarise(MeanRelBiais = mean(RelBiais))


trueVals = c()
k=0
for(i in 1:8){
  trueVals = c(trueVals, rep(results3$Estimate[(k+1):(k+1000)],4))
  k = k + 4000
}
results3$RelBiais = (results3$Estimate - trueVals)/trueVals*100
results3$Var = as.factor(results3$Var)
#levels(results3$Var) = c("Slope for X1", "Slope for X2")
levels(results3$Var) =c("Coefficient estimé pour X1", "Coefficient estimé pour X2")
mean_data3 <- group_by(results3, Rho, Methode, Var) %>%
  summarise(MeanRelBiais = mean(RelBiais))



p1 = ggplot(subset(mean_data1, Methode != "Estimateur avec données sans bruit"), 
            aes(x = Rho, y = MeanRelBiais, colour = Methode)) +
  geom_point() + geom_line(aes(lty = Methode), lwd = 0.8) +facet_wrap(~Var) + 
  theme( text = element_text(size = 15), legend.title = element_text(size=14), 
         axis.text = element_text(size =14)) + geom_hline(yintercept = 0, lwd = 0.8)+ 
  geom_point(size = 2) +ylab("Biais relatif moyen des pentes (en %)") + 
  xlab("Coefficient de corrélation entre X1 et X2")


p2 = ggplot(subset(mean_data2, Methode != "Estimateur avec données sans bruit"), 
            aes(x = Rho, y = MeanRelBiais, colour = Methode)) +
  geom_point() + geom_line(aes(lty = Methode), lwd = 0.8) +facet_wrap(~Var) + 
  theme( text = element_text(size = 15), legend.title = element_text(size=14), 
         axis.text = element_text(size =14)) + geom_hline(yintercept = 0, lwd = 0.8)+ 
  geom_point(size = 2) +ylab("Biais relatif moyen des pentes (en %)") + 
  xlab("Coefficient de corrélation entre X1 et X2")


p3 = ggplot(subset(mean_data3, Methode != "Estimateur avec données sans bruit"), 
            aes(x = Rho, y = MeanRelBiais, colour = Methode)) +
  geom_point() + geom_line(aes(lty = Methode), lwd = 0.8) +facet_wrap(~Var) +
  theme( text = element_text(size = 15), legend.title = element_text(size=14), 
         axis.text = element_text(size =14)) + geom_hline(yintercept = 0, lwd = 0.8)+
  geom_point(size = 2) +ylab("Biais relatif moyen des pentes (en %)") + 
  xlab("Coefficient de corrélation entre X1 et X2")

library(ggpubr)
p1
p2
p3
ggarrange(p1, p2 , p3,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 2, common.legend = T)


#Essayer une autre stratégie - 

newout = rbind(mean_data1, mean_data2, mean_data3)
newout$context = factor(rep(c("Modèle 1", "Modèle 2", "Modèle 3"),each = 64))

library(latex2exp)
gf=ggplot(subset(newout, Methode != "Estimateur avec données sans bruit"), aes(x = Rho, y = MeanRelBiais, colour = Methode)) +
  geom_point() + geom_line(aes(lty = Methode), lwd = 0.8) +facet_wrap(~Var*context) +
  theme( text = element_text(size = 13), legend.title = element_text(size=12), axis.text = element_text(size =10)) + geom_hline(yintercept = 0, lwd = 0.8)+
  geom_point(size = 2) +ylab("Biais relatif moyen des pentes (en %)") +
  xlab(TeX("Coefficient de corrélation entre X1 et X2 , $\\epsilon$ = 20")) + theme(legend.position = "bottom")


gf


#SAVED AS FIGURE RLM2correlatedVar
