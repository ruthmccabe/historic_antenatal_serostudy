# fit incidence to daily hospitalisations and deaths from across sample period

library(tidyverse)
library(Hmisc)
library(cowplot)

incidence <- read.csv("Data/output/incidence_of_seroreactivity.csv") %>% 
  mutate(date_begin = as.Date(date_begin,"%Y-%m-%d")) %>%
  rename(date = date_begin) 


# deaths

# format official death data by replacing NA with 0 and taking rolling weekly mean
deaths <- read.csv("Data/input/covid_official/daily_deaths.csv") %>% mutate(date=as.Date(date)) %>% 
  rename(deaths=newDeaths28DaysByDeathDate) %>% dplyr::select(date,deaths) %>% filter(date<=as.Date("2021-01-07")) %>%
  complete(date = seq.Date(min(incidence$date), max(date), by="day")) %>% 
  mutate(deaths=ifelse(is.na(deaths),0,deaths)) %>%
  mutate(deaths_smooth = zoo::rollmean(deaths, 7, na.pad = TRUE,align="center")) %>% data.frame()

# merge to get time series with deaths and incidence estimates where we have them
deaths_merge <- merge(deaths,incidence,by="date",all=TRUE) %>% arrange(date)

# define lag prior to take between 0 and 40 (days)
lag_prior <- seq(0,40,1)
deaths_merge[,(ncol(deaths_merge)+1):(ncol(deaths_merge)+length(lag_prior))] <- NA

# for each prior value calculate the appropriate lag
for(p in lag_prior){
  deaths_merge[,(7+p+1)] <- Lag(deaths_merge$deaths_smooth,-p)
}

# now reduce to data set dates in which we have incidence estimate
deaths_incidence_model_data <- deaths_merge %>% filter(!is.na(median))
deaths_incidence_model_data$median_100 <- deaths_incidence_model_data$median*100

# estimate variance of each estimate
deaths_incidence_model_data <- deaths_incidence_model_data %>% mutate(variance_est = ((upper-lower)/(2*1.96))^2)

# run regression for each lag
weighted_deaths_lag_results <- c()

for(p in lag_prior){
  model_p <- lm(formula=median_100~-1+deaths_incidence_model_data[,(7+p+1)],
                weights=1/deaths_incidence_model_data$variance_est,
                data=deaths_incidence_model_data)
  results_p <- data.frame("p"=p,"SSE"=deviance(model_p),
                          "scale"=summary(model_p)$coefficients[1,1])
  weighted_deaths_lag_results <- rbind(weighted_deaths_lag_results,results_p)
}


# optimal lag and scale based on minimising SSE
weighted_lag_deaths <- weighted_deaths_lag_results %>% filter(SSE==min(SSE)) %>% select(p) %>% as.numeric()
weighted_scale_deaths <- weighted_deaths_lag_results %>% filter(SSE==min(SSE)) %>% select(scale) %>% as.numeric()


# plot
deaths_plot <- ggplot(deaths_merge %>% filter(date>as.Date("2019-11-04"),date<=as.Date("2020-10-01")),
       aes(x=date))+
  geom_point(aes(y=median*100))+
  geom_errorbar(aes(ymin=lower*100,ymax=upper*100),width=13)+
  theme_bw()+
  geom_line(aes(x=date-weighted_lag_deaths,y=deaths_smooth*weighted_scale_deaths),col="#da68a0",lwd=1)+
  scale_y_continuous(
    name = "Estimated incidence \nof seroreactivity (%)",
    sec.axis = sec_axis(~./weighted_scale_deaths, name="Shifted daily deaths \nin London \n(7 day average)")
  )+
  #ggtitle("No intercept")+
  labs(x="Date",tag="B")+
  scale_x_date(breaks="2 month",date_labels = "%b-%y")


# hospitalisations

# format official hospitalisations data by replacing NA with 0 and taking rolling weekly mean
hosps <- read.csv("Data/input/covid_official/daily_admitted_to_hospital.csv") %>% mutate(date=as.Date(date)) %>% 
  rename(hosps=newAdmissions) %>% dplyr::select(date,hosps) %>% filter(date<=as.Date("2021-01-01")) %>%
  complete(date = seq.Date(min(incidence$date), max(date), by="day")) %>% mutate(hosps=ifelse(is.na(hosps),0,hosps)) %>%
  mutate(hosps_smooth = zoo::rollmean(hosps, 7, na.pad = TRUE,align="center")) %>% data.frame()


# merge to get time series with deaths and incidence estimates where we have them
hosps_merge <- merge(hosps,incidence,by="date",all=TRUE) %>% arrange(date)

# define lag prior to take between 0 and 40 (days)
lag_prior <- seq(0,40,1)
hosps_merge[,(ncol(hosps_merge)+1):(ncol(hosps_merge)+length(lag_prior))] <- NA

# for each prior value calculate the appropriate lag
for(p in lag_prior){
  hosps_merge[,(7+p+1)] <- Lag(hosps_merge$hosps_smooth,-p)
}

# now reduce to data set dates in which we have incidence estimate
hosps_incidence_model_data <- hosps_merge %>% filter(!is.na(median))
hosps_incidence_model_data$median_100 <- hosps_incidence_model_data$median*100

# estimate variance of each estimate
hosps_incidence_model_data <- hosps_incidence_model_data %>% mutate(variance_est = ((upper-lower)/(2*1.96))^2)

# run regression for each lag
weighted_hosps_lag_results <- c()

for(p in lag_prior){
  model_p <- lm(formula=median_100~-1+hosps_incidence_model_data[,(7+p+1)],
                weights=1/hosps_incidence_model_data$variance_est,
                data=hosps_incidence_model_data)
  results_p <- data.frame("p"=p,"SSE"=deviance(model_p),
                          "scale"=summary(model_p)$coefficients[1,1])
  weighted_hosps_lag_results <- rbind(weighted_hosps_lag_results,results_p)
}


# optimal lag and scale based on minimising SSE
weighted_lag_hosps <- weighted_hosps_lag_results %>% filter(SSE==min(SSE)) %>% select(p) %>% as.numeric()
weighted_scale_hosps <- weighted_hosps_lag_results %>% filter(SSE==min(SSE)) %>% select(scale) %>% as.numeric()


# plot
hosps_plot <- ggplot(hosps_merge %>% filter(date>as.Date("2019-11-04"),date<=as.Date("2020-10-01")),
                      aes(x=date))+
  geom_point(aes(y=median*100))+
  geom_errorbar(aes(ymin=lower*100,ymax=upper*100),width=13)+
  theme_bw()+
  geom_line(aes(x=date-weighted_lag_hosps,y=hosps_smooth*weighted_scale_hosps),col="#77c593",lwd=1)+
  scale_y_continuous(
    name = "Estimated incidence \nof seroreactivity (%)",
    sec.axis = sec_axis(~./weighted_scale_hosps, name="Shifted daily hospital \nadmissions in London \n(7 day average)")
  )+
  #ggtitle("No intercept")+
  labs(x="Date",tag="A")+
  scale_x_date(breaks="2 month",date_labels = "%b-%y")


# plot together

plot_grid(hosps_plot,deaths_plot,align="v",nrow=2)
ggsave("data/output/figures/figure3.tiff",height=5,width=6)

