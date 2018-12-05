#/*************************************************/
#/*      STABILITY OF PREDATOR-PREY DYNAMICS      */
#/*                                               */
#/* MEGAN HARTLE, ANDREW MACKINNON, GABRIEL TAURO */
#/*                                               */
#/*        BIOCOMPUTING, DR. STUART JONES         */
#/*************************************************/

# PART 1 
# LOTKA-VOLTERRA MODEL


# PART 1B
# SIMULATE DYNAMICS USING LV MODEL

library(deSolve)  # Includes deSolve library for ode functions
library(ggplot2)  # Includes ggplot2 library for graphing functions

# Creates simulation that models Lotka-Volterra equation
LVSim <- function(t, y, p) {
  H = y[1]  # Herbivore variable
  P = y[2]  # Predator variable

  b = p[1]  # Prey birth rate
  a = p[2]  # Predator attack rate
  e = p[3]  # Conversion efficiency of prey to predators
  s = p[4]  # Predator death rate
  
  dHdt = (b * H) - (a * P * H)  # Differential equations for LV
  dPdt = (e * a * P * H) - (s * P)
  
  return(list(c(dHdt, dPdt)))  # Returns results of derivations 
}

N0 = c(25, 5)  # Initial populations (prey, predator)
params = c(0.5, 0.02, 0.1, 0.2)  # Initial parameters
times = seq(0, 300, by = .1)  # Time steps being used for LV simulation

# Generates ODE function that utilizes LV model
LVModelSim = ode(y = N0, times = times, func = LVSim, parms = params)

out = data.frame(time = LVModelSim[, 1], Prey = LVModelSim[, 2], Predator = LVModelSim[, 3])  # Assigns LVModelSim to a data frame

# Creates a graph comparing predator (red) and prey (cyan) populations
out %>%
  gather(key, value, Prey, Predator) %>%
  ggplot(aes(x = time, y = value, color = key)) + 
  geom_line() + 
  theme_classic() + 
  ggtitle("Lotka-Volterra Model Using Given Initial Values") +
  ylab("Count") + 
  xlab("Time")


# PART 1C
# RUN ADDITIONAL LV SIMULATIONS

# Creates a data frame that includes 9 extra sets of parameters for simulation
LVParameters = data.frame(
  "b" = c(.5, .25, 1, .5, .5, .5, .5, .5, .5),
  "a" = c(.02, .02, .02, .01, .04, .02, .02, .02, .02),
  "e" = c(.1, .1, .1, .1, .1, .05, .2, .1, .1),
  "s" = c(.2, .2, .2, .2, .2, .2, .2, .1, .4)
)

LVModelSimList = list()  # Creates a list for storage of LV simulations

# Generates a for loop that runs the LV simulation using LVParameters
for (i in 1:nrow(LVParameters)) {
  params = unlist(LVParameters[i,])
  LVModelSim2 = ode(y = N0, times = times, func = LVSim, parms = params)
  LVModelSimList[[i]] = data.frame(time = LVModelSim2[, 1], Prey = LVModelSim2[, 2], Predator = LVModelSim2[, 3])
}

# Creates a list of names for respective LV plots
LVNameList = list("Base - LV", "Low B - LV", "High B - LV", "Low A - LV", "High A - LV", "Low E - LV", "High E - LV", "Low S - LV", "High S - LV")

# Generates a for loop that creates and prints a plot for each simulation in LVModelSimList
for (j in 1:length(LVModelSimList)) {
  plot <- LVModelSimList[[j]] %>%
    gather(key, value, Prey, Predator) %>%
    ggplot(aes(x = time, y = value, color = key)) + 
      geom_line() + 
      theme_classic() + 
      ggtitle(LVNameList[[j]]) +  # Calls on LVNameList for plot titles
      ylab("Count") + 
      xlab("Time")
  
  print(plot)
}


# PART 2
# ROSENZWEIG-MACARTHUR MODEL


# PART 2B
# SIMULATE DYNAMICS USING RM MODEL

# Creates simulation that models Rosenzweig-MacArthur equation
RMSim = function(t, y, p) {
  H = y[1]  # Herbivore variable
  P = y[2]  # Predator variable
  
  b = p[1]  # Prey birth rate
  a = p[2]  # Alpha, tied to carrying capacity of prey
  w = p[3]  # Saturating response
  d = p[4]  # Prey density variable at which predator attack rate reaches half its maximum
  e = p[5]  # Conversion efficiency of prey to predators
  s = p[6]  # Predator death rate
  
  dHdt = (b * H) * (1 - a * H) - (w * P) * (H / (d + H))  # Differential equations for RM
  dPdt = (e * w * P) * (H / (d + H)) - s * P
  
  return(list(c(dHdt, dPdt)))  # Returns results of derivations 
}

N0 = c(500, 120)  # Initial populations (prey, predator)
params = c(.8, .001, 5, 400, .07, .2)  # Initial parameters
times = seq(0, 200, by = .1)  # Time steps being used for RM simulation

# Generates ODE function that utilizes RM model
RMModelSim = ode(y = N0, times = times, func = RMSim, parms = params)

out = data.frame(data = RMModelSim, time = RMModelSim[, 1], Prey = RMModelSim[, 2], Predator = RMModelSim[, 3])  # Assigns RMModelSim to a data frame

# Creates a graph comparing predator (red) and prey (cyan) populations
out %>%
  gather(key, value, Prey, Predator) %>%
  ggplot(aes(x = time, y = value, color = key)) + 
    geom_line() + 
    theme_classic() +
    ggtitle("Rosenzweig-MacArthur Model Using Given Initial Values") +  
    ylab("Count") + 
    xlab("Time")


# PART 2C
# RUN ADDTIONAL RM SIMULATIONS

# Creates a data frame that includes 13 extra sets of parameters for simulation
RMParameters = data.frame(
  "b" = c(.8, 1.6, .4, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8),
  "a" = c(.001, .001, .001, .002, .0005, .001, .001, .001, .001, .001, .001, .001, .001),
  "w" = c(5, 5, 5, 5, 5, 10, 2.5, 5, 5, 5, 5, 5, 5),
  "d" = c(400, 400, 400, 400, 400, 400, 400, 800, 200, 400, 400, 400, 400),
  "e" = c(.07, .07, .07, .07, .07, .07, .07, .07, .07, .14, .035, .07, .07),
  "s" = c(.2, .2, .2, .2, .2, .2, .2, .2, .2, .2, .2, .4, .1)
)

RMModelSimList = list()  # Creates a list for storage of RM simulations

# Generates a for loop that runs the LV simulation using RMParameters
for (i in 1:nrow(RMParameters)) {
  params = unlist(RMParameters[i, ])
  RMModelSim2 = ode(y = N0, times = times, func = RMSim, parms = params)
  RMModelSimList[[i]] = data.frame(time = RMModelSim2[, 1], Prey = RMModelSim2[, 2], Predator = RMModelSim2[, 3])
}

# Creates a list of names for respective RM plots
RMNameList = list("Base - RM", "High B - RM", "Low B - RM", "High A - RM", "Low A - RM", "High W - RM", "Low W - RM", "High D - RM", "Low D - RM", "High E - RM", "Low E - RM", "High S - RM", "Low S - RM")

# Generates a for loop that creates and prints a plot for each simulation in RMModelSimList
for (j in 1:length(RMModelSimList)) {
  plot <- RMModelSimList[[j]] %>%
    gather(key, value, Prey, Predator) %>%
    ggplot(aes(x = time, y = value, color = key)) + 
      geom_line() + 
      theme_classic() + 
      ggtitle(RMNameList[[j]]) +  # Calls on RMNameList for plot titles
      ylab("Count") + 
      xlab("Time")
  
  print(plot)
}


# PART 3
# PARADOX OF ENRICHMENT

# Creates a data frame that includes 5 sets of parameters for RM simulations to illustrate Paradox of Enrichment
PoEParameters = data.frame(
  b = c(.8, .8, .8, .8, .8),
  a = c(.00125, .0009, .00075, .0006, .0005),
  w = c(5, 5, 5, 5, 5),
  d = c(400, 400, 400, 400, 400),
  e = c(.1, .1, .1, .1, .1),
  s = c(.2, .2, .2, .2, .2)
)

PoEModelSimList = list()  # Creates a list for storage of PoE simulations


# Generates a for loop that runs the PoE simulation using PoEParameters
for (i in 1:nrow(PoEParameters)) {
  params = unlist(PoEParameters[i,])
  PoEModelSim = ode(y = N0, times = times, func = RMSim, parms = params)
  PoEModelSimList[[i]] = data.frame(time = PoEModelSim[, 1], Prey = PoEModelSim[, 2], Predator = PoEModelSim[, 3])
}

# Creates a list of names for respective PoE plots
PoENameList = list("High Alpha - PoE", "Mid-High Alpha - PoE", "Mid Alpha - PoE", "Mid-Low Alpha - PoE", "Low Alpha - PoE")

# Generates a for loop that creates and prints a plot for each simulation in PoEModelSimList
for (j in 1:length(PoEModelSimList)) {
  plot <- PoEModelSimList[[j]] %>%
    gather(key, value, Prey, Predator) %>%
    ggplot(aes(x = time, y = value, color = key)) + 
      geom_line() + 
      theme_classic() + 
      ggtitle(PoENameList[[j]]) +
      ylab("Count") + 
      xlab("Time")
  
  print(plot)
}