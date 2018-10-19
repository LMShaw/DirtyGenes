#drichlet_LRT_test Single piece of code looking to perform the Dirichlet Likelihood Ratio Test with the option of 
#goodness-of-fit testing on multitype population data. The output will determine whether types play a significant role
#in population structure and whether the Dirichlet model provides a reasonable fit to the observed data, as described in
#Shaw et al. (2017). The power of the test under H_0 can also be estimated.
library(gtools)
library(plyr)
dirichlet.LRT.test = function(df,randomise = 0, ever.present = TRUE, min.proportion = 0, col = 1, gof.sims = NULL, 
                              power.sims = NULL, power.threshold = 0.05, power.plot.obs = NULL){
  #Input should be df: a data frame in which rows denote a population observation, one of the columns 
  #(default column 1) denotes the environment and the remaining columns dentoe number/proportions
  #of each class "type" within the population.
  #randomise: Set to some value > 0 here if there are not enough replications in your dataset to use the 
  #chi-squared test directly. Values represents the number of randomisation trials.
  #ever.present: should only classes that are present in all observations be used to prevent boundary issues
  #with MLE calculation. If TRUE, any classes not meeting this criterion will be aggregated into an "other"
  #class type. Default = TRUE
  #min.proportion sets a a proportion which a class must hit in at least one observation to be considered in the output.
  #Any classes not meeting this criterion with the aggregated into the "other" group. Default = 0.
  #col: The column number in df denoting the environment in which an observed population is found. Default = 1
  #gof.sims: If given stipulates the number of simulations used to obtain goodness-of-fit p-values. Default = NULL 
  #power.sims: If given stipulates the number of simulations used to estimate the power of the test 
  #if MLE under H_1 is true.
  #power.threshold: Determines the default level of significance of the test for power estimations. Default = 0.05
  #power.plot.obs: If given and power.sims is also given, a plot will be produced estimating the power 
  #of the test if the MLE under H_1 is true for each value in the vector power.plot.obs observations of each environment.
  
  optimloglik=function(param,x){
    #Function uses dirichlet loglikelihood for optimisation purposes (negative uses min rather than max)
    names(param) = names(x)
    -sum(log(ddirichlet(x,param)))}
  
  dirLRT=function(df){
  #Function to perform the LRT when data is in required format
  x = df[,-1]
  dummy <- constrOptim(rep(1,ncol(x)), optimloglik, NULL, ui = diag(ncol(x)), ci = rep(0,ncol(x)),outer.eps = 1e-5,x=x)
  L_null = -dummy$value #Log-liklihood at MLE
  null.parameters = dummy$par #MLE paramters
  names(null.parameters) = names(x) #Class names associated with each parameter
  names(df)[1] = "colnam" #Give column a name for using later on
  L_alt = 0 #Starting point for alternative log-likelihood value
  alt.parameters = matrix(nrow = 0, ncol = ncol(x))
  for (i in levels(df$colnam)){
    x = df[df$colnam == i,-1] #Extract data for specific type
    dummy = constrOptim(null.parameters, optimloglik, NULL, ui = diag(ncol(x)), ci = rep(0,ncol(x)),outer.eps = 1e-5,x=x)
    L_alt = L_alt - dummy$value
    alt.parameters = rbind(alt.parameters,dummy$par)
    #Finds log-likelihood under MLE for each factor and sums across each factor to acquire L_alt. 
    #Uses negative to account for negative in optimloglik
  }
  colnames(alt.parameters) = names(null.parameters) #Class names for each paramter
  rownames(alt.parameters) = levels(df$colnam) #Environment of each parameterisation
  LRT.stat = 2*(L_alt-L_null)
  return(list("LRT.stat" = LRT.stat,
              "chisq.p" = pchisq(LRT.stat,(nrow(alt.parameters)-1)*(ncol(alt.parameters)),lower.tail=FALSE),
              "L_null" = L_null,
              "L_alt" = L_alt,
              "null.parameters" = null.parameters,
              "alt.parameters" = alt.parameters)) #Returns the necessary parameters
  }
  
  #Get the data into compositional form, meeting the ever present and minimium proportion criteria set out in the input
  df[is.na(df)] <- 0
  x = df[,-col] #Extracts the data
  x <- x/rowSums(x) #Makes data compositional
  x["LRT.other"] <- 0
  x$LRT.other = x$LRT.other + rowSums(x[, sapply(x, max) < min.proportion])
  x = x[sapply(x, max) >= min.proportion  | colnames(x) == "LRT.other"] #Remove columns not meeting the minimum proportion requirement
  #Deal with ever-present issue
  if (ever.present == TRUE){
    #Want aggregate non-ever-present classes into the "other group"
    x$LRT.other = rowSums(x[lapply(x,min) == 0 | colnames(x) == "LRT.other"]) #Aggregation
    x = x[sapply(x,min) != 0 | colnames(x) == "LRT.other"] #Remove non-ever presents. If LRT.other falls into this catergory it is also removed
  }
  x = x[sapply(x,max) != 0] #Removes LRT.other if nothing is in it
  x <- x/rowSums(x) #Makes data compositional
  x[x == 0] <- sqrt(.Machine$double.eps) #Removes 0s if necessary
  x <- x/rowSums(x) #Reruns composition
  df <- cbind(df[,col],x) #Rebinds with factor column putting factors in col 1 for running the test
  #Can now run the LRT. We put this into a separate function to make life easir with the power testing later
  main = dirLRT(df) #Perform the main test and store results 
  if (randomise){ #Do we need to use randomisation?
    randstat = rep(0,randomise) #Create a vector for storing test stats from randomisation procedure
    for (j in 1:randomise){
      #print(j) Uncomment to track progress
      randstat[j] = dirLRT(transform(df, `df[, col]` = sample(`df[, col]`) ))$LRT.stat
    }
    main$rand.p = sum(randstat > main$LRT.stat)/randomise
  }
  else {
    main$rand.p = NULL
  }
  counts = count(df[,1])$freq #Finds frequency of each environment in the data. Useful in both future tests.
  
  #Have now performed the LRT, move on to gof testing
  if (is.null(gof.sims)){
    #If no entry in the input, do not perform the test or give a result
    null.gof.p = NULL
    alt.gof.p = NULL
  } else {
    #Want to simulate gof.sims copies of the data under both sets of parameters and calculate their likelihood
    null.sims = rep(0,gof.sims)
    alt.sims = null.sims #Sets up vectors to store simulated log-likelihoods
    for (j in 1:gof.sims){
      null.sims[j] = sum(log(ddirichlet(rdirichlet(nrow(df),main$null.parameters),main$null.parameters))) #Null sim
      alt.sims[j] = 0 #See LRT, same method but with a simulation
      for (i in 1:length(counts)){
        alt.sims[j] = alt.sims[j] + sum(log(ddirichlet(rdirichlet(counts[i],main$alt.parameters[i,]),main$alt.parameters[i,])))
        #Finds log-likelihood under MLE for each factor and sums across each factor to acquire simulated MLE
      }
      #Repeat the process the set number of times
    }
    #Determines p-values for gof under both hypotheses
    null.gof.p = sum(main$L_null > null.sims)/gof.sims
    alt.gof.p = sum(main$L_alt > alt.sims)/gof.sims
  }
  
  #Finally consider the power of the test. If power.sims is given, simulate the LRT under that many observations.
  #If power.plot.obs is also given, simulate for up to that number of observations in each environment an produce
  #a plot showing the increasing power of the test
  if (is.null(power.sims)){
    #Dismiss this part if NULL
    power.estimate <- NULL
    p.est <- NULL
  } else {
    #Test the power of the current observation format
    pvals = rep(0,power.sims) #A vector to stor p-values from each simulation
    for (j in 1:power.sims){
      #print(j) #Uncomment to monitor progress
      simdata = as.data.frame(matrix(nrow=0,ncol = ncol(x))) #Set up simulated data frame
      for (i in 1:length(counts)){
        simdata = rbind(simdata,rdirichlet(counts[i],main$alt.parameters[i,]))
        #Simulates the data under the alternative hypothesis
      }
      pvals[j] = dirLRT(cbind(df[,1],simdata))$LRT.p
    }
    power.estimate = sum(pvals < power.threshold)/power.sims
    p.est <- NULL
    if (!is.null(power.plot.obs)){
      #If power.plot.obs is given then also produce a plot
      p.est = rep(0,length(power.plot.obs)) #Place to store estimator of p
      for (k in power.plot.obs){
        #Need two observations of each environment
        pvals = rep(0,power.sims) #A vector to stor p-values from each simulation
        for (j in 1:power.sims){
          #print(c(k,j)) #Uncomment to print current part ofsimulation
          simdata = as.data.frame(matrix(nrow=0,ncol = ncol(x))) #Set up simulated data frame
          for (i in 1:length(counts)){
            simdata = rbind(simdata,rdirichlet(k,main$alt.parameters[i,]))
            #Simulates the data under the alternative hypothesis
          }
          pvals[j] = dirLRT(cbind(rep(levels(df[,1]),each = k),simdata))$LRT.p #Final simulated data frame
        }
        p.est[k-1] = sum(pvals < power.threshold)/power.sims
      }
      plot(power.plot.obs,p.est)
    }
  }
    
  
  
return(list("LRT.stat" = main$LRT.stat, 
            "chisq.p" = main$chisq.p,
            "rand.p" = main$rand.p,
        "null.parameters" = main$null.parameters,
        "alt.parameters" = main$alt.parameters,
        "null.gof.p" = null.gof.p,
        "alt.gof.p" = alt.gof.p,
        "power.estimate" = power.estimate,
        "power.est.vector" = p.est))
  #Returns the test stat and p-value from the LRTs, MLEs for Dirichlet paramters under the null and alternative hypotheses
  #and goodness-of-fit p-values for the Dirichlet distribution under the null and alternative hypotheses. Also gives the results of the power testing
}