require(coda)
source("read.admb.R")

mat2mcmclist<-function(x,burnin=0)
{
    #takes a simple mcmc matrix (x) from admb and converts to an mcmc.list object 
    #to use the coda package functions. Burning discards the specified number of 
    # initial values. THe function currently creates two chains
    #there is no error handling for only 1 parameter
    titer<-length(x[,1])
    useiter<-titer-burnin-1
    jmp<-floor(useiter-1)/2
    break1A<-burnin+1
    break1B<-burnin+1+jmp
    break2A<-break1B+1
    break2B<-break2A+jmp
  
    obj<-list(NA)
    npar<-length(x[1,])
    obj[[1]]<-mcmc(x[break1A:break1B,])
    obj[[2]]<-mcmc(x[break2A:break2B,])
    obj <- mcmc.list(obj)
    return(obj)
}
test <- array(rnorm(10000,0,1),c(1000,10))
gelman.diag(mat2mcmclist(test))
plot(mat2mcmclist(test))


run.Simulation=function(N=10)
{
    theta <<- NULL
    for(i in sample(1:1000,N))
    {
         arg = paste("./vonB -sim", i)
        system(arg)
        print(arg)
        P=read.fit("vonB")
        theta<<-rbind(theta, P$est[1:3])
    }
}
    run.Simulation(10)
    boxplot(exp(theta))    
