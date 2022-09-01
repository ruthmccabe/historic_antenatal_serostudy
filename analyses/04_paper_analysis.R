# analysis in paper

library(tidyverse)
library(binom)

# read in data

data <- read.csv("data/output/clean_data_for_analysis.csv") %>% 
  mutate(date_begin = as.Date(date_begin,"%Y-%m-%d"),
         age_group = factor(age_group,levels=c("18 - 29","30 - 34","35 - 44")),
         ethnicity_group = factor(ethnicity_group,levels=c("All black","All Asian","All white","Other")),
         IMD_group = factor(IMD_group,levels=c("IMD decile group 1 - 2","IMD decile group 3 - 6","IMD decile group 7 - 10")),
         time_fortnight_label = factor(paste0(fortnight,"-",substring(year,3,4))),
         time_fortnight_label = fct_reorder(time_fortnight_label,date_begin,min))

# fortnight 27-19 is only 1 day so merge with 26-19
data$time_fortnight_label[which(data$time_fortnight_label=="27-19")] <- "26-19"

june_2019_data <- read.csv("data/input/june_2019_data.csv")

prevalence <- read.csv("data/output/prevalence_of_seroreactivity.csv") %>% 
  mutate(date_begin = as.Date(date_begin,"%Y-%m-%d"))
prevalence_4_weeks <- read.csv("Data/output/prevalence_of_seroreactivity_per_four_weeks.csv") %>% 
  mutate(date_begin = as.Date(date_begin,"%Y-%m-%d"))

incidence <- read.csv("data/output/incidence_of_seroreactivity.csv") %>% mutate(date_begin = as.Date(date_begin,"%Y-%m-%d"))
incidence_among_susceptibles <- read.csv("data/output/incidence_of_seroreactivity_among_susceptibles.csv") %>% mutate(date_begin = as.Date(date_begin,"%Y-%m-%d"))

react2 <- read.csv("Data/input/REACT2_estimates.csv") %>% mutate(date_begin = as.Date(date_begin,"%d/%m/%Y"))

## analysis

# sample size

nrow(data) #11256

# Table 1: Classification and sample sizes of observations from October 2019 to September 2020 as seroreactive and non-seroreactive according to different thresholds of binding ratio (BR) values. 

#how many above 1.2
sum(data$seroreactive_high_threshold) #665
#how many between 1 and 1.2
sum(data$seroreactive) - sum(data$seroreactive_high_threshold) #35
#how many above 1
sum(data$seroreactive) #700
#how many between 0.8 and 1
sum(data$seroreactive_low_threshold) - sum(data$seroreactive) #70
#how many above 0.8
sum(data$seroreactive_low_threshold) #770
#how many below 0.8
nrow(data) - sum(data$seroreactive_low_threshold) #10486

# overall prevalence of seroreactivity and 95% CI
binom.confint(x=sum(data$seroreactive),n=nrow(data),method="exact")

# sample sizes of demographic variables
data %>% dplyr::select(age_group) %>% summary()
data %>% dplyr::select(ethnicity_group) %>% summary()
data %>% dplyr::select(IMD_group) %>% summary()


# Fisher's test for change in odds of observing seropositive sample

data_with_june <- rbind(data %>% dplyr::select(time_fortnight_label,seroreactive_label) %>% 
                          group_by(time_fortnight_label) %>% 
                          mutate(total=length(time_fortnight_label)),
                        june_2019_data %>% mutate(seroreactive_label=ifelse(label=="Negative",
                                                                            "Non-seroreactive","Seroreactive"),
                                                  total = length(seroreactive_label)) %>% 
                          dplyr::select(time_fortnight_label,seroreactive_label,total)) %>% data.frame()

fortnights <- unique(data$time_fortnight_label)
fisher_test_results_june <- data.frame()

for(i in 1:(length((unique(data$time_fortnight_label))))){
  fn1 <- fortnights[i]
  data_it <- rbind(data_with_june %>% filter(time_fortnight_label==fn1),
                   data_with_june %>% filter(time_fortnight_label=="June-19")) %>% droplevels()
  table_it <- table(data_it$seroreactive_label,data_it$time_fortnight_label)
  test_it <- fisher.test(table_it,alternative = "two.sided")
  fisher_test_results_it <- data.frame("fortnight_comp"=fn1,
                                       "pval"=test_it$p.value,"or_est"=test_it$estimate,
                                       "or_lower"=test_it$conf.int[1],"or_upper"=test_it$conf.int[2])
  fisher_test_results_june <- rbind(fisher_test_results_june,fisher_test_results_it)
}


# Figure 1: Binding ratio values among non-seroreactive (left) and seroreactive (right) observations per fortnight (October 2019 – September 2020) and per month for June 2019. 

# seroreactive observations

pos_plot_data <- rbind(data %>% 
                         filter(seroreactive_label=="Seroreactive") %>%
                         select(result,time_fortnight_label),
                       june_2019_data %>% 
                         filter(label=="Positive") %>% 
                         select(result,time_fortnight_label))
pos_plot_data$seroreactive_label <- factor("Seroreactive")

pos <- ggplot(pos_plot_data %>% filter(time_fortnight_label!="June-19"),aes(x=result,fill=seroreactive_label))+
  geom_boxplot(aes(group=time_fortnight_label),outlier.shape=NA)+
  facet_grid(time_fortnight_label~seroreactive_label)+
  theme_bw()+
  theme(panel.border = element_rect(colour="grey"),
        strip.background = element_rect(fill = "white", color = "grey"),
        axis.line = element_line(color="grey"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        strip.text.y = element_blank(),
        plot.margin=margin(l=-0.5,unit="cm"),
        axis.title.x = element_text(hjust=0.3))+
  scale_fill_manual(values="#4d5198")+
  labs(x="Binding ratio")+scale_x_continuous(expand=c(0,0),breaks=c(1,10,20,30),limits=c(1,33))

pos_june <- ggplot(pos_plot_data %>% filter(time_fortnight_label=="June-19"),aes(x=result,fill=seroreactive_label))+
  geom_boxplot(aes(group=time_fortnight_label),outlier.shape=NA)+
  facet_grid(time_fortnight_label~seroreactive_label)+
  theme_bw()+
  theme(panel.border = element_rect(colour="grey"),
        strip.background = element_rect(fill = "white", color = "grey"),
        axis.line = element_line(color="grey"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        strip.text.y = element_blank(),
        plot.margin=margin(l=-0.5,unit="cm"),
        axis.title.x = element_text(hjust=0.1))+
  scale_fill_manual(values="#4d5198")+
  labs(x=" ")+scale_x_continuous(expand=c(0,0),breaks=c(1,10,20,30),limits=c(1,33))

# non-seroreactive observations 

neg_plot_data <- rbind(data %>% 
                         filter(seroreactive_label=="Non-seroreactive") %>%
                         select(result,time_fortnight_label),
                       june_2019_data %>% 
                         filter(label=="Negative") %>% 
                         select(result,time_fortnight_label)) %>% filter(result>=0)
neg_plot_data <- rbind(neg_plot_data,
                       #to help with plot alignment 
                       setNames(data.frame(1,"June-19"),
                                names(neg_plot_data)))
neg_plot_data$seroreactive_label <- factor("Non-seroreactive")

neg <- ggplot(neg_plot_data %>% filter(time_fortnight_label!="June-19"),aes(x=result,fill=seroreactive_label))+
  geom_boxplot(aes(group=time_fortnight_label),outlier.shape=NA)+
  facet_grid(time_fortnight_label~seroreactive_label,switch="y")+
  theme_bw()+
  theme(panel.border = element_rect(colour="grey"),
        strip.background = element_rect(fill = "white", color = "grey"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_text(hjust=0.1))+
  scale_fill_manual(values="#81b7d2")+
  labs(x="")+scale_x_continuous(expand=c(0,0),breaks = c(0,0.5,1))

neg_june <- ggplot(neg_plot_data %>% filter(time_fortnight_label=="June-19"),aes(x=result,fill=seroreactive_label))+
  geom_boxplot(aes(group=time_fortnight_label),outlier.shape=NA)+
  facet_grid(time_fortnight_label~seroreactive_label,switch="y")+
  theme_bw()+
  theme(panel.border = element_rect(colour="grey"),
        strip.background = element_rect(fill = "white", color = "grey"),
        axis.line = element_line(color="grey"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none",
        strip.text.y.left = element_text(angle = 0,size=7),
        axis.title.x = element_text(hjust=0.1))+
  scale_fill_manual(values="#81b7d2")+
  labs(x=" ")+scale_x_continuous(expand=c(0,0),breaks=c(0,0.5,1))


# combine for june
june_comb <- plot_grid(neg_june,pos_june,align="h",rel_widths=c(1,2))

# combine for main study period
rest_comb <- plot_grid(neg,pos,align="h",rel_widths = c(1,2))

# combine 
plot_grid(june_comb,
          rest_comb,
          nrow=2,rel_heights = c(1,8),align="hv")
ggsave("data/output/figures/figure1.tiff",height=7,width=6)



# in text statistical tests

# # fortnights 12 - 15
# binom.confint(x=prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13"|time_fortnight=="2020-14"|time_fortnight=="2020-15") %>% dplyr::select(reactive) %>% sum(),
#               n=prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13"|time_fortnight=="2020-14"|time_fortnight=="2020-15") %>% dplyr::select(total) %>% sum(),
#               method="exact")
# 
# # fortnights 16 - 19
# binom.confint(x=prevalence %>% filter(time_fortnight=="2020-16"|time_fortnight=="2020-17"|time_fortnight=="2020-18"|time_fortnight=="2020-19") %>% dplyr::select(reactive) %>% sum(),
#               n=prevalence %>% filter(time_fortnight=="2020-16"|time_fortnight=="2020-17"|time_fortnight=="2020-18"|time_fortnight=="2020-19") %>% dplyr::select(total) %>% sum(),
#               method="exact")
# 
# # chi-squared test 
# 
# prop.test(x = c(prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13"|time_fortnight=="2020-14"|time_fortnight=="2020-15") %>% dplyr::select(reactive) %>% sum(),
#                 prevalence %>% filter(time_fortnight=="2020-16"|time_fortnight=="2020-17"|time_fortnight=="2020-18"|time_fortnight=="2020-19") %>% dplyr::select(reactive) %>% sum()),
#           n = c(prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13"|time_fortnight=="2020-14"|time_fortnight=="2020-15") %>% dplyr::select(total) %>% sum(),
#                 prevalence %>% filter(time_fortnight=="2020-16"|time_fortnight=="2020-17"|time_fortnight=="2020-18"|time_fortnight=="2020-19") %>% dplyr::select(total) %>% sum()))



# fortnights 12 - 13
binom.confint(x=prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13") %>% 
                dplyr::select(reactive) %>% sum(),
              n=prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13") %>% 
                dplyr::select(total) %>% sum(),
              method="exact")


# fortnights 18 - 19
binom.confint(x=prevalence %>% filter(time_fortnight=="2020-18"|time_fortnight=="2020-19") %>%
                dplyr::select(reactive) %>% sum(),
              n=prevalence %>% filter(time_fortnight=="2020-18"|time_fortnight=="2020-19") %>% 
                dplyr::select(total) %>% sum(),
              method="exact")


# chi-squared test

prop.test(x = c(prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13") %>% 
                  dplyr::select(reactive) %>% sum(),
                prevalence %>% filter(time_fortnight=="2020-18"|time_fortnight=="2020-19") %>% 
                  dplyr::select(reactive) %>% sum()),
          n = c(prevalence %>% filter(time_fortnight=="2020-12"|time_fortnight=="2020-13") %>% 
                  dplyr::select(total) %>% sum(),
                prevalence %>% filter(time_fortnight=="2020-18"|time_fortnight=="2020-19") %>%
                  dplyr::select(total) %>% sum()))



# Figure 2: Prevalence of seroreactivity over time

prevalence$col_width <- 13
prevalence$col_width[which(prevalence$time=="2019-26")] <- 13.25
prevalence$col_width[which(prevalence$time=="2020-1")] <- 13.25


plot_ests <- ggplot(prevalence,aes(x=date_begin+7,y=mean*100,ymin=lower*100,ymax=upper*100))+
  geom_point()+
  geom_errorbar(aes(width=col_width))+
  theme_bw()+
  scale_x_date(date_breaks = "1.5 month", date_labels = "%b-%y",expand=c(0.01,0.01))+
  coord_cartesian(ylim=c(0, 30))+
  labs(x=" ",y="Estimated \nprevalence of \nseroreactivity (%) \nwith 95% CIs",tag="A")+
  theme(axis.text.x = element_text(angle=0))

plot_sample <- ggplot(prevalence,aes(x=date_begin+7,y=total))+
  geom_col(width=prevalence$col_width,col="white",fill="black")+
  theme_bw()+
  labs(x=" ",y="Number of \nobservations",tag="B")+
  scale_x_date(date_breaks = "1.5 month", date_labels = "%b-%y",expand=c(0.01,0.01))+
  theme(axis.text.x = element_text(angle=0))

plot_positives <- ggplot(prevalence,aes(x=date_begin+7,y=reactive))+
  geom_col(width=prevalence$col_width,col="white",fill="black")+
  theme_bw()+
  labs(x="Date",y="Number of \nseroreactive \nobservations",tag="C")+
  scale_x_date(date_breaks = "1.5 month", date_labels = "%b-%y",expand=c(0.01,0.01))+
  theme(axis.text.x = element_text(angle=0))

cowplot::plot_grid(plot_ests, plot_sample, plot_positives, nrow = 3,align="v")
ggsave("data/output/figures/figure2.tiff",width=6,height=5)


# Table 2: The total number of observations and the number of seroreactive observations per fortnight and per four-week period using a seroreactive threshold of BR value ≥ 1 in the primary study period (and per month in June 2019). 


# prevalence in June 2019
binom.confint(x=length(which(june_2019_data$label=="Positive")),n=nrow(june_2019_data),method="exact")

# fortnightly prevalence
prevalence %>% mutate(estimate_tidy = paste0(round(mean*100,1),
                                                      " (",round(lower*100,1)," - ",round(upper*100,1),")"))

# four weekly prevalence
prevalence_4_weeks %>% mutate(estimate_tidy = paste0(round(mean*100,1),
                                                     " (",round(lower*100,1)," - ",round(upper*100,1),")"))


# Table 3: The number of total observations and the number of seroreactive observations per group within each demography variable across the primary study period (October 2019 – September 2020).

# age
data %>% group_by(age_group) %>% mutate(age_reactive=sum(seroreactive),age_total=length(seroreactive)) %>% 
  select(age_group,age_reactive,age_total) %>% 
  unique() %>% 
  mutate(binom.confint(x=age_reactive,n=age_total,method="exact")) %>% data.frame() %>% select(-x,-n) %>%
  mutate(estimate_tidy = paste0(round(mean*100,1), " (",round(lower*100,1)," - ",round(upper*100,1),")")) %>% 
  arrange(age_group)


# ethnicity
data %>% group_by(ethnicity_group) %>% mutate(ethnicity_reactive=sum(seroreactive),ethnicity_total=length(seroreactive)) %>% 
  select(ethnicity_group,ethnicity_reactive,ethnicity_total) %>% 
  unique() %>% 
  mutate(binom.confint(x=ethnicity_reactive,n=ethnicity_total,method="exact")) %>% data.frame() %>% select(-x,-n) %>%
  mutate(estimate_tidy = paste0(round(mean*100,1), " (",round(lower*100,1)," - ",round(upper*100,1),")")) %>% 
  arrange(ethnicity_group)


# IMD 
data %>% group_by(IMD_group) %>% mutate(IMD_reactive=sum(seroreactive),IMD_total=length(seroreactive)) %>% 
  select(IMD_group,IMD_reactive,IMD_total) %>% 
  unique() %>% 
  mutate(binom.confint(x=IMD_reactive,n=IMD_total,method="exact")) %>% data.frame() %>% select(-x,-n) %>%
  mutate(estimate_tidy = paste0(round(mean*100,1), " (",round(lower*100,1)," - ",round(upper*100,1),")")) %>% 
  arrange(IMD_group)


# multifactorial regression

# prepare data for regression

data_regression <- data %>% 
  mutate(age_group = relevel(age_group,ref="35 - 44"),
         ethnicity_group = relevel(ethnicity_group,ref="All white"),
         IMD_group = relevel(IMD_group,ref="IMD decile group 1 - 2"),
         time_fortnight = factor(time_fortnight),
         time_fortnight = fct_reorder(time_fortnight,date_begin,min),
         time_fortnight = relevel(time_fortnight,ref="2020-6"))


regression_model <- glm(formula = seroreactive ~ age_group + ethnicity_group + IMD_group + time_fortnight,
                        family = "binomial",
                        data = data_regression)

regression_model_summary <- summary(regression_model)$coefficient %>% data.frame()

# Table 4: Coefficient estimates, standard errors and p-values for the logistic regression model with age group, ethnicity, IMD and fortnight of sample as predictors of seroreactivity. 

cbind("Variable" = rownames(regression_model_summary),
      "Estimate" = round(regression_model_summary$Estimate,2),
      "Standard error" = round(regression_model_summary$Std..Error,2),
      "Odds ratio" = round(exp(regression_model_summary$Estimate),2),
      "p-value" = round(regression_model_summary$Pr...z..,3)) %>% data.frame()

# model with interaction between ethnicity and IMD 

regression_model_interaction <- glm(formula = seroreactive ~ age_group + ethnicity_group + IMD_group + ethnicity_group*IMD_group + time_fortnight,
                                    family = "binomial",
                                    data = data_regression)

# overall significance of interaction term
anova(regression_model_interaction,test="Chisq")


# Figure 3: Estimated incidence of seroreactivity over time with 95% bootstrapped confidence intervals. 

# this can be plotted using script 03_fit_incidence_to_hosps_and_deaths


# Supplementary Figure 1: Shape constrained P-spline (red) fit to prevalence of seroreactivity estimated using entire cleaned data set (black). 


prevalence$smooth <- scam(mean ~ s(as.numeric(date_begin), bs="mpi", k=24), data=prevalence)$fitted.values

ggplot(prevalence,aes(x=date_begin + 7,y=mean*100,ymin=lower*100,ymax=upper*100))+
  geom_line(aes(y=smooth*100,col="Shape constrained P-spline"),lwd=1.5)+
  geom_point(aes(fill="Estimated seroprevalence"))+
  geom_errorbar()+
  theme_bw()+
  scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+
  labs(x="Date",y="Estimated prevalence \nof seroreactivity (%) \nwith 95% CIs",col="",fill="")+
  theme(legend.position="bottom")+
  scale_colour_manual(values="red")
ggsave("data/output/figures/supplementaryfigure1.tiff",width=6,height=3)


# Supplementary Figure 2: Estimated incidence of seroreactivity among susceptible persons over time with 95% bootstrapped confidence intervals. 

ggplot(incidence_among_susceptibles,aes(x=date_begin,y=median*100,ymin=lower*100,ymax=upper*100))+
  geom_point()+
  geom_errorbar()+
  theme_bw()+
  scale_x_date(date_breaks = "2 month", date_labels = "%b-%y")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=as.Date("2020-03-23"),linetype="dashed",col="red")+
  labs(x="Date",
       y="Estimated incidence of \nseroreactivity (%) among susceptible \npersons with bootstrapped 95% CIs")
ggsave("data/output/figures/supplementaryfigure2.tiff",width=6,height=3)


# Supplementary Figure 3: Estimates of prevalence of seroreactivity using historic antenatal samples in this study (black) compared to 18 – 44-year-old women across our ethnicity and IMD groups from REACT-2 at a national-level (dark green) and only within London (light green). 

prevalence_comparison <- prevalence %>% filter(date_begin >= as.Date("2020-06-01")) %>% 
  mutate(data = "Historic antenatal samples") %>% select(-smooth)

comparison <- rbind(react2,prevalence_comparison)

ggplot(comparison,
       aes(x=date_begin+7,y=mean*100,ymin=lower*100,ymax=upper*100,col=data))+
  geom_point()+
  geom_errorbar()+
  theme_bw()+
  scale_x_date(date_breaks="1 month",date_labels="%b-%Y")+
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15,0.2,0.25)*100,limits=c(0,27.5),expand=c(0,0))+
  scale_colour_manual(values=c("#000000","#a9c0a6","#097770"))+
  labs(x="Date",y="Estimated prevalence of \nseroreactivity (%) with 95% CIs",col="Data underlying estimates")
ggsave("data/output/figures/supplementaryfigure3.tiff",height=3,width=6)

