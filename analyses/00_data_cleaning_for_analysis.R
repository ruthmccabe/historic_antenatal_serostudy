## prepare data for analysis

library(tidyverse)
library(readxl)

# read in data
data_raw <- read.csv("data/input/data_raw.csv") %>% mutate(result = as.numeric(result)) #warning message is just about introducing NAs

# split study_data to get the necessary demographic information into distinct columns
data <- cbind(data_raw,plyr::ldply(strsplit(data_raw$study_data,"-"), rbind))
colnames(data)[(ncol(data)-4):ncol(data)] <- c("year","fortnight","age","ethnicity","IMD")

# create binary variables for seroreactivity and ambiguity

data <- data %>% mutate(seroreactive = ifelse(result>=1,1,0),
                        ambiguous = ifelse(result>=0.8&result<=1.2,1,0),
                        seroreactive_low_threshold = ifelse(result>=0.8,1,0),
                        seroreactive_high_threshold = ifelse(result>=1.2,1,0))


# only want observations that have all demographic information and a binding ratio (result)

data <- data %>% filter(!is.na(result),
                        year=="2019"|year=="2020",
                        fortnight!="",
                        age=="1"|age=="2"|age=="3",
                        ethnicity=="1"|ethnicity=="2"|ethnicity=="3"|ethnicity=="4",
                        IMD=="1"|IMD=="2"|IMD=="3") 


# create time_fortnight variable
data <- data %>% mutate(time_fortnight = paste0(year,"-",fortnight))

# merge 2019-27 (just one day) with 2019-26
data$time_fortnight[which(data$time_fortnight=="2019-27")] <- "2019-26"

# read in fortnight lookup
date_lookup <- data.frame(read_excel("Data/input/date_lookup.xlsx",sheet="lookup"))
date_lookup$date <- as.Date(date_lookup$date,"%d/%m/%Y") #warning message not relevant

# merge fortnight lookup with data frame to get dates for the fortnights
# date is the first date of the fortnight
data <- merge(data,date_lookup,by="time_fortnight",sort=FALSE) %>% rename(date_begin=date) %>%
  mutate(time_fortnight = fct_reorder(time_fortnight,date_begin,min)) %>% arrange(date_begin)

# label variables for plotting 

clean_data <- data %>% mutate(seroreactive_label = factor(ifelse(seroreactive==1,"Seroreactive","Non-seroreactive"),
                                                          levels=c("Seroreactive","Non-seroreactive")),
                              age_group = factor(ifelse(age==1,"18 - 29",
                                                        ifelse(age==2,"30 - 34","35 - 44"))),
                              ethnicity_group = factor(ifelse(ethnicity==1,"All black",
                                                              ifelse(ethnicity==2,"All Asian",
                                                                     ifelse(ethnicity==3,"All white","Other")))),
                              IMD_group = factor(ifelse(IMD==1,"IMD decile group 1 - 2",
                                                        ifelse(IMD==2,"IMD decile group 3 - 6","IMD decile group 7 - 10"))))


# data is ready for analysis 

write.csv(clean_data,"Data/output/clean_data_for_analysis.csv",row.names=FALSE)









