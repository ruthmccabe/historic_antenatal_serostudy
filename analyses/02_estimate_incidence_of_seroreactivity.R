# fit estimate of spline and then use this to derive incidence and incidence among suspectible persons

library(tidyverse)
library(scam)
library(matrixStats)

prevalence <- read.csv("Data/output/prevalence_of_seroreactivity.csv") %>% 
  mutate(date_begin = as.Date(date_begin,"%Y-%m-%d")) %>%
  arrange(date_begin)

# fit shape constrained P-splines using scam package

prevalence$smooth <- scam(mean ~ s(as.numeric(date_begin), bs="mpi", k=24), data=prevalence)$fitted.values

# look at fit

ggplot(prevalence,aes(x=date_begin + 7,y=mean*100,ymin=lower*100,ymax=upper*100))+
  geom_line(aes(y=smooth*100,col="Shape constrained P-spline"),lwd=1.5)+
  geom_point(aes(fill="Estimated seroprevalence"))+
  geom_errorbar()+
  theme_bw()+
  scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+
  labs(x="Date",y="Estimated prevalence \nof seroreactivity (%) \nwith 95% CIs",col="",fill="")+
  theme(legend.position="bottom")+
  scale_colour_manual(values="red")

# derive the incidence using bootstrapping procedure as explained in paper

data <- read.csv("Data/output/clean_data_for_analysis.csv") %>% mutate(date_begin = as.Date(date_begin,"%Y-%m-%d"))

iter <- 1000
spline_bootstrap <- matrix(NA,nrow=nrow(prevalence),ncol=iter)
incidence_bootstrap <- matrix(NA,nrow=nrow(prevalence),ncol=iter)
incidence_among_susceptibles_bootstrap <- matrix(NA,nrow=nrow(prevalence),ncol=iter)

for(b in 1:iter){
  
  row_sample <- c()
  # sample from each fortnight according to number of observations 
  for(i in unique(data$time_fortnight)){
    rows_fortnights <- which(data$time_fortnight==i)
    rows_fortnights_sample <- sample(rows_fortnights,replace=TRUE,size=length(rows_fortnights))
    row_sample <- c(row_sample,rows_fortnights_sample)
  }
  
  data_sample <- data[row_sample,]
  
  prevalence_sample <- data_sample %>% group_by(time_fortnight) %>% mutate(total = length(seroreactive),
                                                                    reactive = sum(seroreactive)) %>%
    select(time_fortnight,date_begin,total,reactive) %>% unique() %>% 
    mutate(binom.confint(x = reactive,n = total,method="exact")) %>% data.frame() %>%
    select(time_fortnight,date_begin,total,reactive,mean,lower,upper) %>% arrange(date_begin)
  
  prevalence_sample$smooth <- spline_bootstrap[,b] <- scam(mean ~ s(as.numeric(date_begin), bs = "mpi",k=24), data = prevalence_sample)$fitted.values

  prevalence_sample$incidence <- NA
  prevalence_sample$incidence_among_susceptibles <- NA
  
  for(i in 2:nrow(prevalence_sample)){
    prevalence_sample$incidence[i] <- (prevalence_sample$smooth[(i)]-prevalence_sample$smooth[(i-1)])
    prevalence_sample$incidence_among_susceptibles[i] <- (prevalence_sample$smooth[(i)]-prevalence_sample$smooth[(i-1)])/(1-prevalence_sample$smooth[(i-1)])
  }
  
  incidence_bootstrap[,b] <- prevalence_sample$incidence
  incidence_among_susceptibles_bootstrap[,b] <- prevalence_sample$incidence_among_susceptibles
  
  print(b)
}

incidence_estimate <- cbind("date_begin"=prevalence$date_begin,"time_fortnight"=prevalence$time_fortnight,
                            data.frame(rowQuantiles(incidence_bootstrap,probs=c(0.025,0.5,0.975))) %>% 
  rename(lower = X2.5.,median = X50.,upper = X97.5.))
write.csv(incidence_estimate,"Data/output/incidence_of_seroreactivity.csv",row.names=FALSE)


incidence_among_susceptibles_estimate <- cbind("date_begin"=prevalence$date_begin,"time_fortnight"=prevalence$time_fortnight,
                                               data.frame(rowQuantiles(incidence_among_susceptibles_bootstrap,
                                                                       probs=c(0.025,0.5,0.975))) %>% 
                                                 rename(lower = X2.5.,median = X50.,upper = X97.5.))
write.csv(incidence_among_susceptibles_estimate,"Data/output/incidence_of_seroreactivity_among_susceptibles.csv",
          row.names=FALSE)



