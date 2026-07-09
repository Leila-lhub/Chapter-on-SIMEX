
#New global functions which can do all simulations - 

#Y perturbed or not
#Any number of variables
#Any correlation between independent variables
#Let's keep X Uniform though

#Calculer les pentes de régression

my.slopes = function(x,y){
  #x is the matrix of predictors, without the intercept
  n = length(y)
  X = cbind(rep(1,n), x)
  parms = solve( t(X) %*% X) %*% t(X) %*% y
  return(parms[-1])
}

my.slopes2 = function(x,y){
  #x is the matrix of predictors, without the intercept
  n = length(y)
  X = cbind(rep(1,n), x)
  reg = lm(y ~ x)
  s_reg = summary(reg)$coefficients
  parms = c(coef=as.vector(s_reg[,"Estimate"]),var=as.vector(s_reg[,"Std. Error"]^2))
  return(parms)
}

#Une itération de simex pour bruit gaussien
#Need to add noise to all variables, including Y
#ceci nous donne les coefficients - sans intercept- obtenus via le jeu de donnees bruité
##ici p c'est lambda
coeff.simex.gauss.iter = function(x,nbrvar,y,n,p,sd_bruit,perturbY = TRUE){
  x_bruit = x + matrix( sqrt(p)*rnorm(n*nbrvar,mean=0, sd = sd_bruit), nrow = n, ncol = nbrvar ) 
  if(perturbY == TRUE){
    y_bruit = y + sqrt(p)*rnorm(n,mean = 0, sd = sd_bruit)  }else{ y_bruit = y }
  return( my.slopes(x_bruit,y_bruit) )
}

## Ajout de Leila le 2025-05-09 ##
coeff.simex.gauss.iter2 = function(x,nbrvar,y,n,p,sd_bruitX, sd_bruitY,perturbY = TRUE){
  x_bruit = x + matrix( sqrt(p)*rnorm(n*nbrvar,mean=0, sd = sd_bruitX), nrow = n, ncol = nbrvar ) 
  if(perturbY == TRUE){
    y_bruit = y + sqrt(p)*rnorm(n,mean = 0, sd = sd_bruitY)  }else{ y_bruit = y }
  return( my.slopes(x_bruit,y_bruit) )
}
## Fin ajout de Leila ##


# SIMEX pour bruit Gaussien - partie simulation
#Takes as input a matrix (or vector) for x and a vector for y and a lot 
#of parameters 

##la fonction my.simex.gauss est pour la phase de la simulation de la SIMEX
my.simex.gauss = function(x,y,n,ps,k,sd_bruit,B,perturbY){
  #x=xnoisy
  nbrvar = dim(x)[2] 
  if(is.null(nbrvar)) nbrvar = 1  #nbrvar = nbre de var explicatives dans le jeu de données
  slopes = matrix(NA, nrow = k, ncol = nbrvar)
  #k:la taille du vecteur p; p est le parametre 'shape' de la distribution gamma
  #my.simex.lap(xnoisy,ynoisy,n,ps,k,sd_u,B=200, perturbY)
  #y=ynoisy
  for( i in 1:k ){
    temp = replicate(B,coeff.simex.gauss.iter(x,nbrvar,y,n,p=ps[i],sd_bruit, perturbY=perturbY)) 
    if( nbrvar == 1){ 
      slopes[i] = mean(temp)
    }else{
      slopes[i,] = rowMeans(temp)
      }
  }
  return(slopes) #slopes: la matrice des estimateurs des paramètres pour les différentes valeurs de p
}

## Ajout de Leila le 2025-05-09 ##
my.simex.gauss2 = function(x,y,n,ps,k,sd_bruitX, sd_bruitY,B,perturbY){
  #x=xnoisy
  nbrvar = dim(x)[2] 
  if(is.null(nbrvar)) nbrvar = 1  #nbrvar = nbre de var explicatives dans le jeu de données
  slopes = matrix(NA, nrow = k, ncol = nbrvar)
  #k:la taille du vecteur p; p est le parametre 'shape' de la distribution gamma
  #my.simex.lap(xnoisy,ynoisy,n,ps,k,sd_u,B=200, perturbY)
  #y=ynoisy
  for( i in 1:k ){
    temp = replicate(B,coeff.simex.gauss.iter2(x,nbrvar,y,n,p=ps[i],sd_bruitX, sd_bruitY, perturbY=perturbY)) 
    if( nbrvar == 1){ 
      slopes[i] = mean(temp)
    }else{
      slopes[i,] = rowMeans(temp)
    }
  }
  return(slopes) #slopes: la matrice des estimateurs des paramètres pour les différentes valeurs de p
}
## Fin ajout de Leila ##

#Ok, this is the first part
#I do the simex part but not the extrapolation
#That way, I can test different extrapolations later, since simulation is slow

#I hard-coded several things in there, given the simulations I needed to publish
my.simulation.sim = function(betas, perturbY, epsilon,delta, n = 500, SigmaX = NULL,  k=20, nsim = 500){
  
 # betas = c(2)
 # perturbY = TRUE
 # epsilon =15
 # n = 500
 # k=20
 #  nsim = 5
  
  require(mvtnorm)
  
  ps = seq(0,2,length.out = k)
  nvar = length(betas)
  
  slopes_true = matrix(NA,nrow = nsim, ncol = nvar)
  slopes_noisy = matrix(NA,nrow = nsim, ncol = nvar)
  out_simex = array(NA, c(k,nvar,nsim))
  
  for(iter in 1:nsim){
    
    if(is.null(SigmaX)){ #uncorrelated X's
      x = matrix(runif(nvar*n,0,1),ncol = nvar)
    }else{ #use given covariance matrix
      x = matrix(pnorm(rmvnorm(n, mean = c(rep(0,nvar)), sigma = SigmaX)), ncol = nvar) ##ici cest uniquement pour 2 variables?
    }
    
    #n=10;SigmaX=matrix(c(1,0.8,0.8,1),ncol=2)
    
    
    y = 1 + x%*%betas + rnorm(n,0,1)
    
    c=sqrt(2*log(1.25/delta))+0.0001
    Df_2=sqrt(1^2+1^2)
    sd_u = Df_2*c/epsilon 
    
    ystand = (y-min(y)) / (max(y) - min(y)) #Standardise y into 0,1
    
    ##Normale erreur
    xnoisy = x + matrix(rnorm(n*nvar,mean=0,sd=sd_u), ncol = nvar)
    #ynoisy = ystand + rnorm(n,sqrt(2)/sd_u) #2022-11-24 I think this is an error?
    ynoisy = ystand + rnorm(n,mean=0, sd=sd_u)
    
    if(perturbY == TRUE){
        #True slope
        slopes_true[iter,] = my.slopes(x,ystand)
    
        #Naive estimate
        slopes_noisy[iter,] = my.slopes(xnoisy,ynoisy)
    
        #Simex
        out_simex[,,iter] = my.simex.gauss(xnoisy,ynoisy,n,ps,k,sd_u,B=200, perturbY) 
    } else {
        #True slope
        slopes_true[iter,] = my.slopes(x,ystand)
      
        #Naive estimate
        slopes_noisy[iter,] = my.slopes(xnoisy,ystand)
      
        #Simex
        out_simex[,,iter] = my.simex.gauss(xnoisy,ystand,n,ps,k,sd_u,B=200, perturbY)
      
    }

    print(iter)
    
  }
  
  return(list(slopes_true, slopes_noisy, out_simex))
  
}


#Takes as input an object as returned by the function above
#Returns a matrix; each row a method of estimation - each column an iteration
sim_extrap = function(obj, ps){
  #obj=out15

  nsim = dim(obj[[1]])[1]
  nvar = dim(obj[[1]])[2]
  res = matrix(NA, nrow = 4, ncol = nsim*nvar)
  
  #Slopes with real data
  res[1,] = c(obj[[1]]) 
  
  #Slopes with noisy data (naive estimation)
  res[2,] = c(obj[[2]])
  
  #Extrapolations pour le SIMEX -
  #Loop over the number of parameters of interest
  
  #Extrapolation linéaire
  
  betas = matrix(NA,nsim,nvar)
  for(i in 1:nsim){
   for(j in 1:nvar){
     temp = lm(obj[[3]][,j,i] ~ ps)$coeff
     betas[i,j] = temp[1] - 1*temp[2]
       
   }
  }
  res[3,] = c(betas)
  
  
  #Extrapolation quadratique
  betas = matrix(NA,nsim,nvar)
  for(i in 1:nsim){
    for(j in 1:nvar){
      temp = lm(obj[[3]][,j,i] ~ ps + I(ps^2))$coeff
      betas[i,j] = temp[1] - 1*temp[2] + 1*temp[3]
      
    }
  }
  res[4,] = c(betas)
  
  #Removed non-linear extrapolation
  
  return(res)
  
}







