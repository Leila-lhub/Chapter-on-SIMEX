
#RLM - p variables

#Figure 1 - Want to show the impact of the number of variables
#Consider a fixed value of epsilon (10)

#setwd("C:/Users/nombo/OneDrive/Documents/~/Code R pour article 1")
source("Subroutines_comprendre.R")

#Start with uncorrelated...
  
#Did 1 to 5 first; added 6 to 10 after
for(val in 1:10){
  set.seed(54974)
  assign(paste("out",val,sep=""), my.simulation.sim(betas = rep(2,val), perturbY = TRUE, epsilon = 20 , SigmaX = NULL, n = 500, k=20, nsim = 500))
  save(list = paste("out",val,sep =""), file = paste("FichiersRLM/RLM_sim_pvar_20_",gsub("[.]","",val),".Rdata", sep = ""))
  print(val)
}


#Extrapolation part
for(val in 1:10){
  load(paste("FichiersRLM/RLM_sim_pvar_20_",gsub("[.]","",val),".Rdata", sep = ""))
  assign(paste("extrap",val, sep = ""),sim_extrap(get(paste("out",val,sep = "")),ps=seq(0,2,length.out = 20)))
  print(val)
}

#A bit tricky since the number of parameters changes with p...

p = 9  #Change p to see different graphics

#Pour faire les graphiques facilement, mettre le tout dans une grande matrice
results = data.frame(matrix(NA, nrow = 4*500*p, ncol = 4))
names(results) = c("Method", "Var", "Iter", "Estimate")
results$Method = rep(c("true", "noisy", "linear", "quadratic"), each = 500*p)
results$Var = rep( 1:p, each = 500)
results$Iter = rep(1:500, length.out = 4*500*p)

#Populate dataframe ; no loop needed since only one value extrap object
results$Estimate = c(t(get(paste("extrap",p,sep = ""))))
results$Method = factor(results$Method, levels = c("true", "noisy", "linear", "quadratic"))


#Figures - 

library(ggplot2)

p <- ggplot(results, aes(x=Method, y=Estimate, fill = Method)) + 
  geom_boxplot() + labs(x="Method", y = "Estimate")
p


#Ok, now get the necessary information for the plot I want to publish

out = data.frame(matrix(NA, nrow = 40, ncol = 3))
names(out) = c("NbrVar", "Method", "MeanRelBias")
out$NbrVar = rep(1:10,each = 4)
out$Method = rep(c("true", "noisy", "linear", "quadratic"), 10)

#Need to loop over the extrap objects and compute the 4 quantities needed
#For method = true, the result should be zero

temp = c()
for( i in 1:10 ){
  
  obj = get(paste("extrap",i,sep = ""))
  obj2 = sweep(obj, 2, obj[1,])
  obj3 = sweep(obj2,2, obj[1,], FUN = "/")
  vec = rowMeans(obj3)*100
  temp = c(temp,vec)
}

out$MeanRelBias = temp
out #Ok, ça c'est les biais moyens sur toutes les simulations

ggplot(out, aes(x = NbrVar, y = MeanRelBias, colour = Method)) +
  geom_point() + geom_line()

#Ok, et si je voulais plutôt un boxplot

out2 = data.frame(matrix(NA, nrow = 40*500, ncol = 4))
names(out2) = c("NbrVar", "Method", "Iter", "MeanRelBias")
out2$NbrVar = rep(1:10,each = 4*500)
out2$Method = rep(rep(c("true", "noisy", "linear", "quadratic"), each = 500), 10)
out2$Iter = rep(1:500, 40)

#Need to loop over the extrap objects and compute the 4 quantities needed
#For method = true, the result should be zero

temp = c()
for( i in 1:10 ){
  
  obj = get(paste("extrap",i,sep = ""))
  obj2 = sweep(obj, 2, obj[1,])
  obj3 = sweep(obj2,2, obj[1,], FUN = "/")
  vec1 = sapply( split( obj3[1,], rep(1:500, each = i)), mean)*100
  vec2 = sapply( split( obj3[2,], rep(1:500, each = i)), mean)*100
  vec3 = sapply( split( obj3[3,], rep(1:500, each = i)), mean)*100
  vec4 = sapply( split( obj3[4,], rep(1:500, each = i)), mean)*100
  temp = c(temp,vec1,vec2,vec3,vec4)
}

out2$MeanRelBias = temp
out2

ggplot(subset(out2, Method %in% c("noisy", "quadratic")), aes(x=as.factor(NbrVar), y=MeanRelBias, fill = Method)) +   geom_boxplot() + labs(x="Number of (uncorrelated) predictors", y = "Mean relative bias over all estimates slopes") + geom_hline(yintercept = 0)

#Graph final
#Je pense que je vais plutôt mettre le premier...

out3 = out
out3$Method = factor(out$Method, levels = c("true", "noisy", "linear", "quadratic"))
#levels(out3$Method) = c("true", "Naive estimate from noisy data","SIMEX with linear extrapolation", "SIMEX with quadratic extrapolation")

out3$Methode=out3$Method
levels(out3$Methode) = c("Estimateur avec données sans bruit", "Estimateur naif avec données bruitées", "SIMEX avec extrapolation linéaire", "SIMEX avec extrapolation quadratique", "non-linéaire")
head(out3)
ggplot(subset(out3, Methode != "Estimateur avec données sans bruit"), aes(x = NbrVar, y = MeanRelBias, colour = Methode)) +
  geom_point() + geom_line(aes(lty = Methode), lwd = 0.8) + theme( text = element_text(size = 15),legend.position = c(0.7, 0.6), legend.title = element_text(size=14), axis.text = element_text(size =14)) + geom_hline(yintercept = 0)+ geom_point(size = 2) +
  ylab("Moyenne du bias relatif sur les pentes (en %)")  +scale_x_continuous(name = "Nombre de prédicteurs (non correlés)", breaks = 1:10, labels = 1:10)

#SAVED AS FIGURE 2
