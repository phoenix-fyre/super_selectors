## Empty workspace ##
rm(list = ls())

## Necessary R packages to load ##
library(deSolve)
library(tidyverse)

# here is your model put into a function
model <- function(t,state,parameters){ # input is time, state of variables, and parameter values
	with(as.list(c(state,parameters)),{
		dS <- S*g*(1-m) -d*S*I -e*A*S # dynamics of sensitive cells
		dR <- R*(g-c) -d*I*R -(1-x)*e*A*R + g*m*S # dynamics of resistant cells
		dI <- b - f*I + r*I*(S+R) # dynamics of the immune system
		# return rates of change
		list(c(dS,dR,dI))
	})
}


timestep<-0.001
t1<-seq(0,20,by=timestep) # t1 is time pre-treatment
t2<-seq(max(t1)+timestep,120,by=timestep) # t2 is time during treatment
times<-c(t1,t2) # make one time vector for plotting

### Values for our loop to measure total resistance
r_vals <- seq(0.01,1,length.out=50) # values to loop over
total_resistance<-NULL # to store results
max_resistance <- NULL # to store the maximum value results

for (i in 1:length(r_vals)){
  
# parameter values
 g <- 1 # bacterial growth rate
  d <- 0.5 # immune induced death rate
  e <- 1 # efficacy of antibiotic
  c <- 0.30 # growth cost of resistance
  r <- r_vals[i] # immune system recruitment rate
  m <- 0.01 # mutation probability
  A_pre <- 0 # antibiotic concentration before treatment
  A_post <- 0.5 # antibiotic concentration during treatment
  b <- 0.1 # baseline immune system recruitment
  f <- 0.3 # immune system decay rate
  x <- 1 # extent of resistance (1=complete, 0=none)
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

total_resistance[i]<-sum(results[,3]) # total resistance across the infection
max_resistance[i] <- max(results[,3])

}
#temp_tibble <- total_resistance_tibble
resistance_tibble <- tibble(total_resistance = total_resistance,
                                  max_resistance = max_resistance,
                                  r_vals = r_vals)
#all_results <- tibble()
#all_results <- rbind(all_results, temp_tibble)

tail(resistance_tibble)

ggplot(data = resistance_tibble, aes(x= r_vals, y = log10(total_resistance))) +
  geom_point(aes(color = log(total_resistance)), size = 1.5) +
  scale_color_continuous(name = "Resistance", low = "#2100d9", high = "#4ecfe6", guide = "colourbar") +
  labs(x = "Immune Cell Recruitment Rate", y = "Log 10 of Total Resistance")

ggplot(data = resistance_tibble, aes(x= A_vals, y = max_resistance)) +
  geom_point(aes(color = max_resistance), size = 1.5) +
  scale_color_continuous(name = "Resistance", low = "#46dd3e", high = "#113d0e") +
  labs(x = "Antibiotic Concentration", y = "Maximum Resistance")

##########Using a different model integrating differential immune death and antibiotic degradation ##########
model <- function(t,state,parameters){ # input is time, state of variables, and parameter values
  with(as.list(c(state,parameters)),{
    
    dS <- S*g*(1-m) -d_S*S*I -e*(A*sin(2*pi*freq*t+phase)+A)*S # dynamics of sensitive cells
    dR <- R*(g-c) -d_R*I*R -(1-x)*e*(A*sin(2*pi*freq*t+phase)+A)*R + g*m*S # dynamics of resistant cells
    dI <- b - f*I + r*I*(S+R) # dynamics of the immune system
    # return rates of change
    list(c(dS,dR,dI))
  })
}


timestep<-0.001
t1<-seq(0,20,by=timestep) # t1 is time pre-treatment
t2<-seq(max(t1)+timestep,72,by=timestep) # t2 is time during treatment
times<-c(t1,t2) # make one time vector for plotting

#immune_strength <- 
r_vals <- seq(0.1,0.9, length.out = 17)
for(j in 1:length(r_vals))
{
### Values for our loop to measure resistance against baseline immune system recruitment
b_vals <- seq(0.1,0.3,length.out=21) # values of r to loop over
total_resistance<-NULL # to store results
max_resistance <- NULL # to store the maximum value results
avg_resistance <- NULL #To store average frequency per run
relative_resistance <- NULL
max_relative_resistance <-NULL
relative_sensitive <- NULL

for (i in 1:length(b_vals)){
  
  # parameter values
  g <- 1 # bacterial growth rate
  d_S <- 0.6 # immune induced death rate for sensitive cells
  d_R <- 0.4 # immune induced death rate for resistant cells
  e <- 1 # efficacy of antibiotic
  c <- 0.50 # growth cost of resistance
  r <- 0.75 # immune system recruitment rate
  m <- 0.01 # mutation probability
  A_pre <- 0 # antibiotic concentration before treatment
  A_post <- 0.5 # antibiotic concentration during treatment
  b <- b_vals[i] # baseline immune system recruitment
  f <- 0.3 # immune system decay rate
  x <- 1 # extent of resistance (1=complete, 0=none)
  I_baseline<-b/f # baseline immune system level in absence of pathogen
  freq <- 0.2
  phase <- 0
  
  S_init<-0.01 # initial number of sensitive cells
  R_init<-0 # initial number of resistant cells
  I_init<- I_baseline # initial number of immune cells
  
  parameters<-c(g=g,d_S=d_S,d_R=d_R,e=e,c=c,r=r,A=A_pre,m=m,x=x,b=b,f=f) # put parameters into a vector
  state<-c(S= S_init,R=R_init,I=I_init) # put initial states in vector
  output_pre_treatment<-ode(y=state,times=t1,func=model,parms=parameters) # run the model for pre-treatment period
  
  parameters<-c(g=g,d_S=d_S,d_R=d_R,e=e,c=c,r=r,A=A_post,m=m,x=x,b=b,f=f) # put parameters into a vector, now adding antibiotic
  state<-c(S=as.numeric(output_pre_treatment[length(t1),2]),R=as.numeric(output_pre_treatment[length(t1),3]),I=as.numeric(output_pre_treatment[length(t1),4])) # initial state is final state of pre-treatment
  output_during_treatment<-ode(y=state,times=t2,func=model,parms=parameters) # run the model again
  
  
  results<-rbind(output_pre_treatment, output_during_treatment) # link all the results together
  
  total_resistance[i]<-sum(results[,3]) # total resistance across the infection
  max_resistance[i] <- max(results[,3])
  avg_resistance[i] <- mean(results[,3])
  relative_resistance[i] <- mean(results[,3]/(results[,3] + results[,2]))
  max_relative_resistance[i] <- max(results[,3]/(results[,3] + results[,2]))
  relative_sensitive[i] <- mean(results[,2]/(results[,3] + results[,2]))
}


total_resistance_tibble <- tibble(total_resistance = total_resistance,
                                  max_resistance = max_resistance,
                                  avg_resistance = avg_resistance,
                                  relative_resistance = relative_resistance,
                                  max_relative_resistance=max_relative_resistance,
                                  relative_sensitive=relative_sensitive,
                                  b_vals = b_vals)

(total_plot <- ggplot(data = total_resistance_tibble, aes(x= b_vals, y = log(total_resistance))) +
  geom_point(aes(color = log(total_resistance)), size = 1.5) +
  scale_color_continuous(name = "Resistance", low = "#ff9b9b", high = "#ff3030") +
  labs(x = "Baseline Immune Recruitment", y = "Log 10 of Total Resistance"))

pdf(paste0("immune_recruitment_",r_vals[j], "_total_resistance",".pdf"))
print(total_plot)
dev.off()

(freq_plot <- ggplot(data = total_resistance_tibble, aes(x= b_vals, y = relative_resistance)) +
  geom_point(color = "firebrick1", size = 1.5) +
  geom_point(aes(x = b_vals, y = relative_sensitive), color = "dodgerblue1") +
  #scale_color_continuous(name = "Resistance", low = "#46dd3e", high = "#113d0e") +
  labs(x = "Baseline Immune Recruitment", y = "Average Relative Frequency"))

pdf(paste0("immune_recruitment_",r_vals[j],"_avg_rel_freq",".pdf"))
print(freq_plot)
dev.off()

}

ggplot(data = total_resistance_tibble, aes(x= A_vals, y = avg_resistance)) +
  geom_point(aes(color = avg_resistance), size = 1.5) +
  scale_color_continuous(name = "Resistance", low = "#ff9b9b", high = "#ff3030") +
  labs(x = "Antibiotic Concentration", y = "Average Resistance Frequency")


ggplot(data = total_resistance_tibble, aes(x= A_vals, y = max_relative_resistance)) +
  geom_point(aes(color = max_resistance), size = 1.5) +
  scale_color_continuous(name = "Resistance", low = "#46dd3e", high = "#113d0e") +
  labs(x = "Antibiotic Concentration", y = "Average Relative Resistance Frequency")


### Unequal immune effect causing increased resistance with increased immunity






