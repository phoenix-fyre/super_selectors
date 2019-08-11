### Packages to load ###
library(deSolve)
library(tidyverse)
library("reshape2")

### The model of resistant bacteria emergence within hosts put into a R function
model <- function(t,state,parameters){ # input is time, state of variables, and parameter values
	with(as.list(c(state,parameters)),{
		
		dS <- S*g*(1-m) -d*S*I -e*A*S # dynamics of sensitive cells
		dR <- R*(g-c) -d*I*R -(1-x)*e*A*R + g*m*S # dynamics of resistant cells
		dI <- b - f*I + r*I*(S+R) # dynamics of the immune system
		# return rates of change
		list(c(dS,dR,dI))
	})
}

### This is a non-deterministic model ###
timestep<-0.001
t1<-seq(0,20,by=timestep) # t1 is time pre-treatment
t2<-seq(max(t1)+timestep,120,by=timestep) # t2 is time during treatment
times<-c(t1,t2) # make one time vector for plotting

d_vals <- seq(0.8,0.9, by = .02)
e_vals <- rev(seq(0.6,0.8, by = .04))
c_vals <- seq(0.3,0.5, by = .04)

for(i in 1:length(d_vals))
{
# parameter values
  g<-1.2 # bacterial growth rate
  d<-d_vals[i] # immune induced death rate
  e<-e_vals[i] # efficacy of antibiotic
  c<-c_vals[i] # growth cost of resistance
  r<-0.4 # immune system recruitment rate
  m<-0.01 # mutation probability
  A_pre<-0 # antibiotic concentration before treatment
  A_post<-0.5 # antibiotic concentration during treatment
  b<-0.15 # baseline immune system recruitment
  f<-0.3 # immune system decay rate
  x<-1 # extent of resistance (1=complete, 0=none)
  I_baseline<-b/f # baseline immune system level in absence of pathogen

S_init<-0.01 # initial number of sensitive cells
R_init<-0 # initial number of resistant cells
I_init<- I_baseline # initial number of immune cells

parameters<-c(g=g,d=d,e=e,c=c,r=r,A=A_pre,m=m,x=x,b=b,f=f) # put parameters into a vector
state<-c(S= S_init,R=R_init,I=I_init) # put initial states in vector
output_pre_treatment<-ode(y=state,times=t1,func=model,parms=parameters) # run the model for pre-treatment period

parameters<-c(g=g,d=d,e=e,c=c,r=r,A=A_post,m=m,x=x,b=b,f=f) # put parameters into a vector, now adding antibiotic
state<-c(S=as.numeric(output_pre_treatment[length(t1),2]),R=as.numeric(output_pre_treatment[length(t1),3]),I=as.numeric(output_pre_treatment[length(t1),4])) # initial state is final state of pre-treatment
output_during_treatment<-ode(y=state,times=t2,func=model,parms=parameters) # run the model again


results<-rbind(output_pre_treatment, output_during_treatment) # link all the results together


####### Plot our results utilizing ggplot2 ########
##Make our data into a tibble for graphing and manipulation
results_tibble <- tibble(times = results[,1],
                        sensitive = results[,2],
                        resistant = results[,3],
                        immune = results[,4]) 
##Convert this tibble into 'long' format for easier graphing using fnct melt from pckg reshape2
rslt_tbl_long <- melt(results_tibble, id= "times") ##Convert data to long format
tail(results_tibble, n= 100)
### the ggplot figure itself
ggplot(data = rslt_tbl_long, aes(x = times, y = value, colour =variable)) +
  geom_line() +
  labs(x = "Time Steps",  y = "Cell Density") +
  scale_color_manual(values = c("dodgerblue1","firebrick1","goldenrod1"))
  
name <- paste('comp_plot',i,'.jpeg',sep = '_')

#file.path("C:","Users","Zachary","Documents","Edinburgh","bioinf","Plots","gif_pics",name)

ggsave(filename = file.path("C:","Users","Zachary","Documents","Edinburgh","bioinf","Plots","gif_pics",name),
       device = "jpeg")

}