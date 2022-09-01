# calculate estimates of seroreactivity per fortnight with 95% exact binomial confidence intervals 

library(tidyverse)
library(binom)

clean_data <- read.csv("Data/output/clean_data_for_analysis.csv")

prevalence <- clean_data %>% group_by(time_fortnight) %>% mutate(total = length(seroreactive),
                                                                     reactive = sum(seroreactive)) %>%
  select(time_fortnight,date_begin,total,reactive) %>% unique() %>% 
  mutate(binom.confint(x = reactive,n = total,method="exact")) %>% data.frame() %>%
  select(time_fortnight,date_begin,total,reactive,mean,lower,upper)

write.csv(prevalence,"Data/output/prevalence_of_seroreactivity.csv",row.names=FALSE)

# prevalence per 4 week period (for table only)

prevalence_4_weeks <- cbind(prevalence,"four_week"=rep(seq(1,12,1),each=2)) %>% select(date_begin,total,reactive,four_week) %>%
  group_by(four_week) %>% mutate(total_4_week=sum(total),reactive_4_week=sum(reactive)) %>% 
  select(date_begin,four_week,total_4_week,reactive_4_week) %>%
  rename(total=total_4_week,reactive=reactive_4_week) %>% filter(!duplicated(four_week)) %>% 
  mutate(binom.confint(x = reactive,n = total,method="exact")) %>% data.frame() %>%
  select(date_begin,total,reactive,mean,lower,upper)

write.csv(prevalence_4_weeks,"Data/output/prevalence_of_seroreactivity_per_four_weeks.csv",row.names=FALSE)
