library(tidyverse)
library(truncnorm)
library(gridExtra)

#set working directory to location of folder
setwd("~/GitHub/oer-simulation-study/")

#based on this code
#http://egap.org/content/power-analysis-simulations-r


possible.ns <- c(100,500,1000,2000,3000,4000,5000,10000,15000,20000) #  The sample sizes we'll be considering
powers <- rep(NA, length(possible.ns))           # Empty object to collect simulation estimates
effects <- rep(NA, length(possible.ns))
alpha <- 0.05                                    # Standard significance level
sims <- 10000                                      # Number of simulations to conduct for each N
set.seed(1000)


simulateOER <- function(has.book = .6,access.effect = .25, mean.grade = 70, sd.grade = 20, allocation.ratio = .5){
  for (j in 1:length(possible.ns)){
    N <- possible.ns[j]                              # Pick the jth value for N
    
    significant.experiments <- rep(NA, sims)         # Empty object to count significant experiments
    
    #### Inner loop to conduct experiments "sims" times over for each N ####
    for (i in 1:sims){
      Y0 <-  rnorm(n=N,mean=mean.grade, sd=sd.grade)             # OER group potential outcome
      Y0[Y0 > 100] <- 100 # Trim down anything over 100 to 100. 
      Y0[Y0 < 0] <- 0 # Trim down anything over 100 to 100. 
      
      b.rate <- rbinom(n=N,size = 1,prob = has.book) # determine how many people in the sample have the book
      
      # Control group
      Y1 <- Y0 - (1-b.rate)*(access.effect*sd.grade)            # treatment potential outcome. For only students who didn't have access to book, improve by access effect
      Y1[Y1 > 100] <- 100 # Trim down anything over 100 to 100. 
      
      condition <- rbinom(n=N, size=1, prob=allocation.ratio)          # Do a random assignment to condition (1 = treatment, 0 = control)
      Y.sim <- Y1*condition + Y0*(1-condition)               # Reveal outcomes according to assignment
      fit.sim <- lm(Y.sim ~ condition)                   # Do analysis (Simple regression)
      positive <- summary(fit.sim)$coefficients[2,1] < 0 # indicates that the non-oer group was worse than oer group
      p.value <- summary(fit.sim)$coefficients[2,4]  # Extract p-values
      significant.experiments[i] <- ((p.value <= alpha)*positive) # Determine significance according to p <= 0.05. Only successful if the direction is positive. 
      # significant.experiments[i] <- (p.value <= alpha) # Determine significance according to p <= 0.05
    }
    
    powers[j] <- mean(significant.experiments)       # store average success rate (power) for each N
  }
  return(powers)
}


simulateSingleExperiment <- function(N = 1000, has.book = .6,access.effect = .25, mean.grade = 70, sd.grade = 20, allocation.ratio = .5, sims = 1000){
  significant.experiments <- rep(NA, sims)         # Empty object to count significant experiments
    #### Inner loop to conduct experiments "sims" times over for each N ####
    for (i in 1:sims){
      Y0 <-  rnorm(n=N,mean=mean.grade, sd=sd.grade)             # control potential outcome
      Y0[Y0 > 100] <- 100 # Trim down anything over 100 to 100. 
      Y0[Y0 < 0] <- 0 # Trim down anything over 100 to 100. 
      b.rate <- rbinom(n=N,size = 1,prob = has.book) # determine how many people in the sample have the book
      Y1 <- Y0 - (1-b.rate)*(access.effect*sd.grade)            # treatment potential outcome. For only students who didn't have access to book, improve by access effect
      Y1[Y1 < 0] <- 0 #adjust anyone who happened to drop below 0
      condition <- rbinom(n=N, size=1, prob=allocation.ratio)          # Do a random assignment to condition (1 = treatment, 0 = control)
      Y.sim <- Y1*condition + Y0*(1-condition)               # Reveal outcomes according to assignment
      fit.sim <- lm(Y.sim ~ condition)                   # Do analysis (Simple regression)
      positive <- summary(fit.sim)$coefficients[2,1] < 0
      p.value <- summary(fit.sim)$coefficients[2,4]  # Extract p-values
      significant.experiments[i] <- ((p.value <= alpha)*positive) # Determine significance according to p <= 0.05. Only successful if the direction is positive. 
    }
    
  power <- mean(significant.experiments)       # store average success rate (power) for each N
  return(power)
}


# run simulations across range of access levels
run_sims <- function(){
  sim90_2 <- simulateOER(has.book = .9,access.effect = .25)
  print("90% done")
  sim80_2 <- simulateOER(has.book = .8,access.effect = .25)
  print("80% done")
  sim70_2 <- simulateOER(has.book = .7,access.effect = .25)
  print("70% done")
  sim60_2 <- simulateOER(has.book = .6,access.effect = .25)
  print("60% done")
  sim50_2 <- simulateOER(has.book = .5,access.effect = .25)
  print("50% done")
  sim40_2 <- simulateOER(has.book = .4,access.effect = .25)
  print("40% done")
  sim0_2 <- simulateOER(has.book = 0,access.effect = .25)
  print("0% done")

  results_1 <- data.frame(possible.ns,sim90_2,sim80_2,sim70_2,sim60_2,sim50_2,sim40_2,sim0_2)
  write_csv(results_1, path = "simulation_results.csv")
}


# Run simulations of past studies
evaluate_past_studies <- function(){
  papers <- read_csv(file = "OER papers.csv")
  
  papers <- papers %>% 
    mutate(allocation = `n treatment`/(`n control`+ `n treatment`))
  
  papers$N <- ifelse(is.na(papers$`n treatment`),papers$`n control`,papers$`n control`+papers$`n treatment`)
  papers$allocation[is.na(papers$allocation)] <- .5 # assume equal allocation if not reported
  
  for (i in 1:nrow(papers)){
     papers$power_40[i] <- simulateSingleExperiment(N = papers$N[i], allocation.ratio = papers$allocation[i], has.book = .4)
     papers$power_60[i] <- simulateSingleExperiment(N = papers$N[i], allocation.ratio = papers$allocation[i], has.book = .6)
     papers$power_80[i] <- simulateSingleExperiment(N = papers$N[i], allocation.ratio = papers$allocation[i], has.book = .8)
    print(paste(i,"of",nrow(papers)))  
  }
  write_csv(papers, path = "paper_analysis.csv")
}


# create the Simulation Plot
results_1 <- read_csv(file = "simulation_results.csv")
names(results_1) <- c("possible.ns","90%", "80%", "70%", "60%", "50%", "40%")
rates <- c("0%","40%", "50%", "60%", "70%", "80%", "90%")

plot1 <- results_1 %>% 
  filter(possible.ns <= 5000) %>% 
  gather(key = sim,value = "power",-possible.ns) %>% 
  ggplot(., aes(x = possible.ns,y = power,group = sim)) +
  geom_point(aes(shape = sim), size = 3) +
  geom_line(aes(linetype = sim)) +
  ylim(0,1) +
  scale_x_continuous(breaks = seq(from=0, to=5000, by= 500)) +
  ylab("Probability of Rejecting Null (Power)")+ xlab("Sample Size (n)") +
  # theme_minimal() + 
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 10)) +
  # theme( legend.position = "bottom",legend.direction = "horizontal") +
  guides(shape=guide_legend(title="Access Rate"),
         # guides(shape=guide_legend(expression(beta)),
         linetype = FALSE)


########### Sensitivity Analysis: ############
# Find the minimum effect size to achieve 80% power
# across varying sample sizes, increment d until power of 80% is achieved. 
getSensitivity <- function(has.book = .6){
  sens <- c()
  for(y in 1:length(possible.ns)){
    for(i in 1:300){ # 300 = 3 SD. 
      power <- simulateSingleExperiment(N = possible.ns[y], has.book = has.book, access.effect = i/100, sims = 1000)
      if(power >= .80){
        print(i)
        sens[y] <- i/100
        break
      }
    }
  }
  return(sens)
}

getSensitivityBinary <- function(has.book = .6){
  sens <- c()
  U <- 3
  for(y in 1:length(possible.ns)){
    print(possible.ns[y])
    print(U)
    print('')
    V <- seq(0, U, .01)
    L <- 1
    R <- length(V)
    while (L < R) {
      m <- floor((L + R) / 2)
      power <- simulateSingleExperiment(N = possible.ns[y], has.book = has.book, access.effect = V[m], sims = 1000)
      if (power < .8) {
        L <- m+1
      }
      else if (power >= .8) {
        R <- m
      }
    }
    sens[y] <- V[L]
    U <- V[L]
  }
  return(sens)
}
    

runSens <- function(){
  sens_90 <- getSensitivity(has.book = .90)
  sens_80 <- getSensitivity(has.book = .80)
  sens_70 <- getSensitivity(has.book = .70)
  sens_60 <- getSensitivity(has.book = .60)
  sens_50 <- getSensitivity(has.book = .50)
  sens_40 <- getSensitivity(has.book = .40)
  sens_results <- data.frame(possible.ns,sens_90,sens_80,sens_70,sens_60,sens_50,sens_40)
  write_csv(sens_results, path = "sensitivity_results.csv")
}

sens_results <- read_csv(file = "sensitivity_results.csv")
names(sens_results) <- c("possible.ns","90%", "80%", "70%", "60%", "50%", "40%")
plot_sens <- sens_results %>% 
  filter(possible.ns <= 5000) %>% 
  gather(key = sim,value = "ES",-possible.ns) %>% 
  ggplot(., aes(x = possible.ns,y = ES,group = sim)) +
  geom_point(aes(shape = sim), size = 3) +
  geom_line(aes(linetype = sim)) +
  ylim(0,2) +
  scale_x_continuous(breaks = seq(from=0, to=5000, by= 500)) +
  ylab("Minimum Effect Size")+ xlab("Sample Size (n)") +
  # theme_minimal() + 
  theme_classic()+
  theme(aspect.ratio = 1,
        text = element_text(size = 10)) +
  # theme( legend.position = "bottom",legend.direction = "horizontal") +
  guides(shape=guide_legend(title="Access Rate"),
         # guides(shape=guide_legend(expression(beta)),
         linetype = FALSE)





# for(i in 1:length(rates)){
#   plot <- results_1 %>% 
#     filter(possible.ns <= 5000) %>% 
#     gather(key = sim,value = "power",-possible.ns) %>%
#     filter(sim %in% rates[1:i]) %>% 
#     ggplot(., aes(x = possible.ns,y = power,group = sim)) +
#     geom_hline(yintercept = .8, linetype = 7,color = "tomato") + 
#     geom_point(aes(shape = sim), size = 3) +
#     geom_line(aes(linetype = sim)) +
#     ylim(0,1) +
#     scale_x_continuous(breaks = seq(from=0, to=5000, by= 500)) +
#     ylab("Probability of Rejecting Null (Power)")+ xlab("Sample Size (n)") +
#     # theme_minimal() + 
#     theme_classic()+
#     theme(aspect.ratio = 1) +
#     guides(shape=guide_legend(title="Access Rate"),
#            linetype = FALSE) 
#   ggsave(filename = paste0("plot",i,".png"),plot = plot,device = "png")
# }
