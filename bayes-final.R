##Import data
data= read.delim("C:/Users/ADMIN/Desktop/Bayes_test.txt", header=FALSE)
###Model M0
t0=Sys.time()#timing
##MH algorithm
#iterate 800 times
n=800
#Set four different starting points
beta=array(dim=c(4,n,3))
beta[1,1,]=c(1,0,0)
beta[2,1,]=c(1,1,0.1)
beta[3,1,]=c(0.5,0.5,0.5)
beta[4,1,]=c(-0.1,-0.1,0)

set.seed(9)
 for (k in 1:4) {
      for (t in 2:n) {
        #beta0,beta1,beta2 are independent, let normal distribution be jumping distribution,
        #the mean is the previous beta and variance is 0.5^2
        beta[k,t,1]=rnorm(1,beta[k,t-1,1],0.5)
        beta[k,t,2]=rnorm(1,beta[k,t-1,2],0.5)
        beta[k,t,3]=rnorm(1,beta[k,t-1,3],0.5)
        beta_post=array(dim=c(1000,1))
        beta_prev=array(dim=c(1000,1))
        #compute p(beta*|y)
        for (i in 1:1000){
          beta_post[i]=(exp(beta[k,t,1]+beta[k,t,2]*data$V1[i]+beta[k,t,3]*data$V2[i]))^data$V3[i]/(1+exp(beta[k,t,1]+beta[k,t,2]*data$V1[i]+beta[k,t,3]*data$V2[i]))
        }
       like_post=prod(beta_post)
      
       #compute p(beta^(t-1)|y)
        for (i in 1:1000){
          beta_prev[i]=(exp(beta[k,t-1,1]+beta[k,t-1,2]*data$V1[i]+beta[k,t-1,3]*data$V2[i]))^data$V3[i]/(1+exp(beta[k,t-1,1]+beta[k,t-1,2]*data$V1[i]+beta[k,t-1,3]*data$V2[i]))
          }
        like_prev=prod(beta_prev)
        #ratio of density
        r=like_post/like_prev
        #draw r from uniform(1)
        u=runif(1)
        if (u>r) {
          beta[k,t,1:3]=beta[k,t-1,1:3]
        }
      }
 }
t1=Sys.time()
t_m0=t1-t0#time of iterating process
#plot the last 200 iterations to see to convergence.
par(mfrow=c(1,2))
plot(as.vector(beta[1:4,601:800,2])~as.vector(beta[1:4,601:800,1]),xlab="beta0",ylab="beta1",xlim=c(-4,4),ylim=c(-4,4))
plot(as.vector(beta[1:4,601:800,3])~as.vector(beta[1:4,601:800,1]),xlab="beta0",ylab="beta2",xlim=c(-4,4),ylim=c(-4,4))
#use the mean of last 200 iterations as approximate beta

beta0_post=mean(beta[1:4,601:800,1])
beta1_post=mean(beta[1:4,601:800,2])
beta2_post=mean(beta[1:4,601:800,3])

###DIC
t2=Sys.time()#timing
#compute log(p(y|beta_hat))
beta_hat=array(c(1000,1))
for (i in 1:1000){
  beta_hat[i]=(exp(beta0_post+beta1_post*data$V1[i]+beta2_post*data$V2[i]))^data$V3[i]/(1+exp(beta0_post+beta1_post*data$V1[i]+beta2_post*data$V2[i]))
}
log_like_beta_hat=log(prod(beta_hat))
#compute pDIC, use the last 200 iterations as simulation data.
beta_sim=array(c(1000,1))
log_like_beta_sim=array(c(200,1))
sim_mean_k=c(1,1,1,1)
for (k in 1:4) {
  for (t in 601:800) {
    for (i in 1:1000){
     beta_sim[i]=(exp(beta[k,t,1]+beta[k,t,2]*data$V1[i]+beta[k,t,3]*data$V2[i]))^data$V3[i]/(1+exp(beta[k,t,1]+beta[k,t,2]*data$V1[i]+beta[k,t,3]*data$V2[i]))
    }
    log_like_beta_sim[t-600]=log(prod(beta_sim))
  }
  sim_mean_k[k]=mean(log_like_beta_sim)
}
sim_mean=mean(sim_mean_k)

pdic=2*(log_like_beta_hat-sim_mean)
#DIC of M0
dic=-2*log_like_beta_hat+2*pdic
t3=Sys.time()
t_dic_m0=t3-t2#time of computing dic
####Model M1
##MH of M1
n=800
beta_m1=array(dim=c(4,n,4))
beta_m1[1,1,]=c(1,1,1,1)
beta_m1[2,1,]=c(1,1,0.1,0.1)
beta_m1[3,1,]=c(0.5,0.5,0.5,0.5)
beta_m1[4,1,]=c(-0.1,-0.1,0,1)

set.seed(9)

for (k in 1:4) {
  for (t in 2:n) {
    #beta0,beta1,beta2, beta12 are independent, let normal distribution be jumping distribution,
    #the mean is the previous beta and variance is 0.5^2
    beta_m1[k,t,1]=rnorm(1,beta_m1[k,t-1,1],0.5)
    beta_m1[k,t,2]=rnorm(1,beta_m1[k,t-1,2],0.5)
    beta_m1[k,t,3]=rnorm(1,beta_m1[k,t-1,3],0.5)
    beta_m1[k,t,4]=rnorm(1,beta_m1[k,t-1,4],0.5)
    beta_post_m1=array(dim=c(1000,1))
    beta_prev_m1=array(dim=c(1000,1))
    #compute p(beta*|y)
    for (i in 1:1000){
      beta_post_m1[i]=(exp(beta_m1[k,t,1]+beta_m1[k,t,2]*data$V1[i]+beta_m1[k,t,3]*data$V2[i]+beta_m1[k,t,4]*data$V1[i]*data$V2[i]))^data$V3[i]/
        (1+exp(beta_m1[k,t,1]+beta_m1[k,t,2]*data$V1[i]+beta_m1[k,t,3]*data$V2[i]+beta_m1[k,t,4]*data$V1[i]*data$V2[i]))
    }
    like_post_m1=prod(beta_post_m1)
    
    #compute p(beta^(t-1)|y)
    for (i in 1:1000){
      beta_prev_m1[i]=(exp(beta_m1[k,t-1,1]+beta_m1[k,t-1,2]*data$V1[i]+beta_m1[k,t-1,3]*data$V2[i]+beta_m1[k,t-1,4]*data$V1[i]*data$V2[i]))^data$V3[i]/
        (1+exp(beta_m1[k,t-1,1]+beta_m1[k,t-1,2]*data$V1[i]+beta_m1[k,t-1,3]*data$V2[i]+beta_m1[k,t-1,4]*data$V1[i]*data$V2[i]))
    }
    like_prev_m1=prod(beta_prev_m1)
    #ratio of density
    r_m1=like_post_m1/like_prev_m1
    #draw r from uniform(1)
    u=runif(1)
    if (u>r_m1) {
      beta_m1[k,t,1:4]=beta_m1[k,t-1,1:4]
    }
  }
}
#plot the last 200 iterations to see to convergence.
par(mfrow=c(1,3))
plot(as.vector(beta_m1[1:4,601:800,2])~as.vector(beta_m1[1:4,601:800,1]),xlab="beta0",ylab="beta1",xlim=c(-4,4),ylim=c(-4,4))
plot(as.vector(beta_m1[1:4,601:800,3])~as.vector(beta_m1[1:4,601:800,1]),xlab="beta0",ylab="beta2",xlim=c(-4,4),ylim=c(-4,4))
plot(as.vector(beta_m1[1:4,601:800,4])~as.vector(beta_m1[1:4,601:800,1]),xlab="beta0",ylab="beta12",xlim=c(-4,4),ylim=c(-4,4))
#use the mean of last 200 iterations as approximate beta

beta0_post_m1=mean(beta_m1[1:4,601:800,1])
beta1_post_m1=mean(beta_m1[1:4,601:800,2])
beta2_post_m1=mean(beta_m1[1:4,601:800,3])
beta12_post_m1=mean(beta_m1[1:4,601:800,4])
##DIC of M1
#compute log(p(y|beta_hat))
beta_hat_m1=array(c(1000,1))
for (i in 1:1000){
  beta_hat_m1[i]=(exp(beta0_post_m1+beta1_post_m1*data$V1[i]+beta2_post_m1*data$V2[i]+beta12_post_m1*data$V1[i]*data$V2[i]))^data$V3[i]/
    (1+exp(beta0_post_m1+beta1_post_m1*data$V1[i]+beta2_post_m1*data$V2[i]+beta12_post_m1*data$V1[i]*data$V2[i]))
}
log_like_beta_hat_m1=log(prod(beta_hat_m1))
#compute pDIC, use the last 200 iterations as simulation data.
beta_sim_m1=array(c(1000,1))
log_like_beta_sim_m1=array(c(200,1))
sim_mean_k_m1=c(1,1,1,1)
for (k in 1:4) {
  for (t in 601:800) {
    for (i in 1:1000){
      beta_sim_m1[i]=(exp(beta_m1[k,t,1]+beta_m1[k,t,2]*data$V1[i]+beta_m1[k,t,3]*data$V2[i]+beta_m1[k,t,4]*data$V1[i]*data$V2[i]))^data$V3[i]/
        (1+exp(beta_m1[k,t,1]+beta_m1[k,t,2]*data$V1[i]+beta_m1[k,t,3]*data$V2[i]+beta_m1[k,t,4]*data$V1[i]*data$V2[i]))
    }
    log_like_beta_sim_m1[t-600]=log(prod(beta_sim_m1))
  }
  sim_mean_k_m1[k]=mean(log_like_beta_sim_m1)
}
sim_mean_m1=mean(sim_mean_k_m1)

pdic_m1=2*(log_like_beta_hat_m1-sim_mean_m1)

#DIC of M1
dic_m1=-2*log_like_beta_hat_m1+2*pdic_m1

###glm model
#model0
t4=Sys.time()#timing
summary(glm(data$V3~data$V1+data$V2,family="binomial"))
t5=Sys.time()
t_glm_m0=t5-t4
#model1
summary(glm(data$V3~data$V1+data$V2+data$V1*data$V2,family="binomial"))

