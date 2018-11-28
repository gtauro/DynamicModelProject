rm(list = ls())

#pt. 2 of final project
library(ggplot2)
library(tidyr)
library(deSolve)

ddSim=function(t,y,p){
  H = y[1]
  P = y[2]
  b=p[1]
  a=p[2]
  w=p[3]
  d=p[4]
  e=p[5]
  s=p[6]
  
  
  dHdt=(b*H)*(1-a*H)-(w*P)*(H/(d+H))
  dPdt=(e*w*P)*(H/(d+H))-s*P
  
  return(list(c(dHdt,dPdt)))
}
# Specific given initial parameters
params = c(.8,.001,5,400,.07,.2)
NO = c(2,2)
times = seq(0,30, by = .1)

modelSim = ode(y=NO, times = times, func = ddSim, parms = params)
out= data.frame(data=modelSim, time = modelSim[,1], prey=modelSim[,2], pred=modelSim[,3])

out %>%
  gather(key,value, prey, pred) %>%
  ggplot(aes(x=time, y=value, color=key)) + geom_line()


# modeled with changing parameters


# params =c(b,a,w,d,e,s)
# b=list(c(b,b/m,b*m,b,b,b,b,b,b,b,b,b,b))
# a=list(c(a,a,a,a/m,a*m,a,a,a,a,a,a,a,a))
# w=list(c(w,w,w,w,w,w,w,w,w,w,w,w/m,w*m))
# d=list(c(d,d,d,d,d,d,d,d,d,d/m,d*m,d,d))
# e=list(c(e,e,e,e,e,e/m,e*m,e,e,e,e,e,e))
# s=list(c(s,s,s,s,s,s,s,s/m,s*m,s,s,s,s))

parameters2 = data.frame("b"=c(.8,1.6,.4,.8,.8,.8,.8,.8,.8,.8,.8,.8,.8), 
                         "a"=c(.001,.001,.001,.002,.0005,.001,.001,.001,.001,.001,.001,.001,.001),
                         "w"=c(.5,.5,.5,.5,.5,1,.25,.5,.5,.5,.5,.5,.5),
                         "d"=c(400,400,400,400,400,400,400,800,200,400,400,400,400),
                         "e"=c(.07,.07,.07,.07,.07,.07,.07,.07,.07,1.4,.035,.07,.07),
                         "s"=c(.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.2,.4,.1))

modelH_Output=data.frame("t", "H", "bHhig", "bHlow", "aHhigh", "aHlow", "eHhigh", "eHlow", "sHhigh", "sHlow", "dHhigh", "dHlow", "wHhigh", "wHlow")
modelP_Output=data.frame("t", "P", "bPhigh", "bPlow", "aPhigh", "aPlow", "ePhigh", "ePlo", "sPhigh", "sPlow", "dPhigh", "dPlow", "wPhigh", "wPlow")

modelH_Output[modelH_Output$X.t.==times]
modelP_Output[modelP_Output$X.t.==times]

modelSimList = list()

for(i in 1:nrow(parameters2)){
  params = parameters2[i,]
  modelSim2 = ode(y=NO, times = times, func = ddSim, parms = params)
  modelSimList[[i]] = data.frame(time = modelSim2[,1], prey=modelSim2[,2], pred=modelSim2[,3])
}

for(j in 1:length(modelSimList)){
  plot <- modelSimList[[j]] %>%
    gather(key,value, prey, pred) %>%
    ggplot(aes(x=time, y=value, color=key)) + geom_line()
  print(plot)
}