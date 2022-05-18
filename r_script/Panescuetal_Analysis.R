##Created by Marvin Browne May 2022
##
##Here we perform statistical analyses(anova and repeated measures anova) for the publication: 
##Effects of Trehalose and Polyacrylate-Based Hydrogels on Tomato Growth Under Drought
##Authors: Priera H. Panescu, Marvin Browne, Kathleen K. Chen, Lawren Sack, and Heather D. Maynard
##
##Notes: The data may be accessed from the excel file within the repository but also from the 
##included Rds file and loaded at line 50. 

# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(readxl)
library(multcomp)
library(lme4)
library(lmerTest)
      

se<-function(x){sd(x)/sqrt(length(x))}

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Data --------------------------------------------------------------------

#load growth data
(tre_rgr_data<-readxl::read_xlsx(here("hydro_growth.xlsx"),
                             sheet = "rgr_trehalose")) 

#load repeated physio measurements
(tre_physio_data<-read_xlsx(here("hydro_growth.xlsx"), 
                      sheet = "physio_trehalose"))

tre_physio_data<-tre_physio_data%>%
        mutate(group=substring(individual, 1,2))

#load repeated physio measurements
(com_rgr_data<-read_xlsx(here("hydro_growth.xlsx"), 
                                  sheet = "rgr_commercial"))
#load repeated physio measurements
(com_physio_data<-read_xlsx(here("hydro_growth.xlsx"), 
                                  sheet = "physio_commercial"))

com_physio_data<-com_physio_data %>%
          mutate(group=substring(individual, 1,2))

#load data from Rdata

load(here("Panescuetal_Rdata.rds"))

com_rgr_data<-com_rgr_data %>%
  mutate(group=substring(individual, 1,2))

tre_rgr_data<-tre_rgr_data %>%
  mutate(group=substring(individual, 1,2))

##summarize means for RGR 
tre_rgr_mean<-tre_rgr_data%>%
  dplyr::select(4:14)%>%
  group_by(group)%>%
  summarise_at(.vars=1:10,
               #group=unique(group),
               .funs=mean)

com_rgr_mean<-com_rgr_data%>%
  dplyr::select(4:14)%>%
  group_by(group)%>%
  summarise_at(.vars=1:10,
               #group=unique(group),
               .funs=mean)
#SE
tre_rgr_se<-tre_rgr_data%>%
  dplyr::select(4:14)%>%
  group_by(group)%>%
  summarise_at(.vars=1:10,
               #group=unique(group),
               .funs=se)

com_rgr_se<-com_rgr_data%>%
  dplyr::select(4:14)%>%
  group_by(group)%>%
  summarise_at(.vars=1:10,
               #group=unique(group),
               .funs=se)


# ANOVA analysis(with RGR, and with all the others) ----------------------------------------------------------

#convert gel to factor
com_rgr_data$gel_nogel<-as.factor(com_rgr_data$gel_nogel)

tre_rgr_data$gel_nogel<-as.factor(tre_rgr_data$gel_nogel)

com_rgr_data$treatment<-as.factor(com_rgr_data$treatment)

tre_rgr_data$treatment<-as.factor(tre_rgr_data$treatment)

#RGR
rgr_treat_gel_c<-aov(formula = rgr~treatment+gel_nogel, 
                     data = com_rgr_data)

rgr_treat_gel_t<-aov(formula = rgr~treatment+gel_nogel, 
                     data = tre_rgr_data)

summary(rgr_treat_gel_c)
summary(rgr_treat_gel_t)

#L_mass
l_treat_gel_c<-aov(formula = leaf_mass~treatment+gel_nogel, 
                     data = com_rgr_data)

l_treat_gel_t<-aov(formula = leaf_mass~treatment+gel_nogel, 
                     data = tre_rgr_data)

summary(l_treat_gel_c)
summary(l_treat_gel_t)

#R_mass
r_treat_gel_c<-aov(formula = root_mass~treatment+gel_nogel, 
                     data = com_rgr_data)

r_treat_gel_t<-aov(formula = root_mass~treatment+gel_nogel, 
                     data = tre_rgr_data)

summary(r_treat_gel_c)
summary(r_treat_gel_t)

#Re_mass
re_treat_gel_c<-aov(formula = repro_mass~treatment+gel_nogel, 
                      data = com_rgr_data)

re_treat_gel_t<-aov(formula = repro_mass~treatment+gel_nogel, 
                      data = tre_rgr_data)

summary(re_treat_gel_c)
summary(re_treat_gel_t)

#S_mass
s_treat_gel_c<-aov(formula = stem_mass~treatment+gel_nogel, 
                     data = com_rgr_data)

s_treat_gel_t<-aov(formula = stem_mass~treatment+gel_nogel, 
                     data = tre_rgr_data)

summary(s_treat_gel_c)
summary(s_treat_gel_t)


# Repeated measures ANOVA -------------------------------------------------
##trehalose repeated measures
psi_leaf_tre<- lmer(psi_leaf ~ treatment*time + gel_nogel*time + (1|individual), data=tre_physio_data)
spad_tre<- lmer(spad ~ treatment*time + gel_nogel*time + (1|individual), data=tre_physio_data)
gs_tre<-lmer(gs ~ treatment*time + gel_nogel*time + (1|individual), data=tre_physio_data)

##commercial repeated measures
psi_leaf_com<- lmer(psi_leaf ~ treatment*time + gel_nogel*time + (1|individual), data=com_physio_data)
spad_com<- lmer(spad ~treatment*time + gel_nogel*time + (1|individual), data=com_physio_data)


#seperate aov for commerical gel treatment gs
gs_com<-aov(gs ~ treatment+ gel_nogel , data=com_physio_data)
TukeyHSD(gs_com)


rm_anovas<-c(psi_leaf_com,psi_leaf_tre, spad_com,spad_tre, gs_com, gs_tre)

capture.output(summary(psi_leaf_com),file = here("rm_anovas", "psi_leaf_com.txt"))

capture.output(summary(spad_com),file = here("rm_anovas", "spad_com.csv"))

capture.output(summary(gs_com),file = here("rm_anovas", "gs_com.txt"))

capture.output(summary(psi_leaf_tre),file = here("rm_anovas", "psi_leaf_tre.txt"))

capture.output(summary(spad_tre),file = here("rm_anovas", "spad_tre.txt"))

capture.output(summary(gs_tre),file = here("rm_anovas", "gs_tre.txt"))


# Post-hoc ----------------------------------------------------------------
growth_objects<-ls(pattern = "treat_gel")

growth_tukey<-lapply(growth_objects,
                     function(object) {print(object); TukeyHSD(eval(parse(text= object)))})
