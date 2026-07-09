
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

#Une itération de simex pour bruit Laplace
#Need to add noise to all variables, including Y
#ceci nous donne les coefficients - sans intercept- obtenus via le jeu de donnees bruité
coeff.simex.lap.iter = function(x,nbrvar,y,n,p,sd_bruit,perturbY = TRUE){
  x_bruit = x + matrix( rgamma(n*nbrvar,shape = p, scale = sd_bruit/sqrt(2)), nrow = n, ncol = nbrvar) - matrix( rgamma(n*nbrvar,shape = p, scale = sd_bruit/sqrt(2)), nrow = n, ncol = nbrvar)  
  if(perturbY == TRUE){
    y_bruit = y + rgamma(n,shape = p, scale = sd_bruit/sqrt(2)) - rgamma(n,shape = p, scale = sd_bruit/sqrt(2))
  }else{
    y_bruit = y
  }
  return( my.slopes(x_bruit,y_bruit) )
}

# SIMEX pour bruit Laplace - partie simulation
#Takes as input a matrix (or vector) for x and a vector for y and a lot 
#of parameters 

my.simex.lap = function(x,y,n,ps,k,sd_bruit,B,perturbY){
  
  nbrvar = dim(x)[2] 
  if(is.null(nbrvar)) nbrvar = 1  #nbrvar = nbre de var explicatives dans le jeu de données
  slopes = matrix(NA, nrow = k, ncol = nbrvar)
  #k:la taille du vecteur p; p est le parametre 'shape' de la distribution gamma
  #my.simex.lap(xnoisy,ynoisy,n,ps,k,sd_u,B=200, perturbY)
  
  for( i in 1:k ){
    temp = replicate(B,coeff.simex.lap.iter(x,nbrvar,y,n,ps[i],sd_bruit, perturbY)) 
    if( nbrvar == 1){ 
      slopes[i] = mean(temp)
    }else{
      slopes[i,] = rowMeans(temp)
      }
  }
  return(slopes) #slopes: la matrice des estimateurs des paramètres pour les différentes valeurs de p
}

#Ok, this is the first part
#I do the simex part but not the extrapolation
#That way, I can test different extrapolations later, since simulation is slow

#I hard-coded several things in there, given the simulations I needed to publish
my.simulation.sim = function(betas, perturbY, epsilon, n = 500, SigmaX = NULL,  k=20, nsim = 500){
  
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
    
    sd_u = 1*(nvar+1)*sqrt(2)/epsilon 
    ystand = (y-min(y)) / (max(y) - min(y)) #Standardise y into 0,1
    #plusieurs gamma à la fois --> exponentielle
    ##laplace erreur
    xnoisy = x + matrix(rexp(n*nvar,sqrt(2)/sd_u), ncol = nvar) - matrix(rexp(n*nvar,sqrt(2)/sd_u), ncol = nvar)
    ynoisy = ystand + rexp(n,sqrt(2)/sd_u)-rexp(n,sqrt(2)/sd_u)
    
    #True slope
    slopes_true[iter,] = my.slopes(x,ystand)
    
    #Naive estimate
    slopes_noisy[iter,] = my.slopes(xnoisy,ynoisy)
    
    #Simex
    out_simex[,,iter] = my.simex.lap(xnoisy,ynoisy,n,ps,k,sd_u,B=200, perturbY) 

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







