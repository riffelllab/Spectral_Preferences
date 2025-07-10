#Imports necessary packages
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(stringr)
library(forcats)
library(readr)
library(emmeans)
#sets dplyr version of select to be called with "select()"
select<-dplyr::select

#Loads custom R functions
source("R code/wind tunnel functions.R")

#Sets the output directory
output_dir<-"../figures/"


#------------------------------------------------------
#-------- initial data import and processing ----------
#------------------------------------------------------

#Sets the subdirectory to search for "resp_tracks.csv" track files
sub_dir<-"output files/color vs gray 0.50, air, CO2-01, CO2-05, CO2-10/"

#---- Track level files
#Creates a list of "tracks data (resp).csv" track files
track.files<-
  list.files(sub_dir, pattern="tracks data \\(resp\\).csv", recursive = TRUE)
#Filters out track files moved into the "~exclude" folder
track.files<-
  grep("~exclude", track.files, value=T, invert=T)

#Creates emtpy dataframe
tracks<-NULL
#Reads in and rbinds each csv file
for(i in seq(1, length(track.files), 1)){
  tracks_i<-
    read_csv(paste0(sub_dir,track.files[i]), 
             #keeps all columns as characters
             col_types = cols(.default = "c"))
  tracks<-rbind(tracks, tracks_i)
  remove(tracks_i)
}

#Displays the total number of trajectories for each event
#across all the runs
table(tracks$event)
table(tracks$stim_series)


#Creates a subset of data and processes it for analysis and display
tracks<-
  #starts with "tracks"
  tracks %>%
  #Filters out the pre and post CO2 times
  filter(!str_detect(event, "CO2Time")) %>%
  #Creates CO2 factor
  mutate(CO2 = case_when(
    str_detect(stim_series, "clean air") ~ "0%",
    str_detect(stim_series, "CO2 01") ~ "1%",
    str_detect(stim_series, "CO2 05") ~ "5%",
    str_detect(stim_series, "CO2 10") ~ "10%",
    TRUE ~ NA), 
    .after = stim_series) %>%
  #Turns "trt_color" into a factor 
  mutate(trt_color = factor(trt_color, levels = c("gray","430","525","625","all")),
         CO2 = factor(CO2, levels=c("0%","1%", "5%", "10%"))) %>%
  #turns rep into a factor
  mutate(rep = as.factor(rep)) %>%
  #turns these variables numeric
  mutate_at(c("pref_index", "duration", "trt_time", 
              "ctr_time", "stim_time"), as.numeric)


#Displays the total number of trajectories recruited for 
#each experimental treatment across all the runs
(table.summary<-with(
      tracks, table(paste(stim_series, trt_color, sep=" - "), rep) ) )
#Displays the row and column means
as.data.frame(round(rowMeans(table.summary),0))
as.data.frame(round(colMeans(table.summary),0))



#---- Event level files
#Creates a list of "tracks data (resp).csv" track files
event.files<-
  list.files(sub_dir, pattern="event data.csv", recursive = TRUE)
#Filters out track files moved into the "~exclude" folder
event.files<-
  grep("~exclude", event.files, value=T, invert=T)

#Creates emtpy dataframe
events<-NULL
#Reads in and rbinds each csv file
for(i in seq(1, length(event.files), 1)){
  event_i<-
    read_csv(paste0(sub_dir,event.files[i]), 
             #keeps all columns as characters
             col_types = cols(.default = "c"))
  #row binds rows even if the columns don't match exactly
  events<-bind_rows(events, event_i)
  remove(event_i)
}


#Creates a subset of data and processes it for analysis and display
events<-
  #starts with "events"
  events %>%
  #Filters out the pre and post CO2 times
  filter(!str_detect(event, "CO2Time")) %>%
  #Creates CO2 factor
  mutate(CO2 = case_when(
    str_detect(stim_series, "clean air") ~ "0%",
    str_detect(stim_series, "CO2 01") ~ "1%",
    str_detect(stim_series, "CO2 05") ~ "5%",
    str_detect(stim_series, "CO2 10") ~ "10%",
    TRUE ~ NA), 
    .after = stim_series) %>%
  #Turns "trt_color" into a factor 
  mutate(trt_color = factor(trt_color, levels = c("gray","430","525","625","all")),
         CO2 = factor(CO2, levels=c("0%","1%", "5%", "10%"))) %>%
  #turns rep into a factor
  mutate(rep = as.factor(rep)) %>%
  #turns these variables numeric
  mutate_at(c("activity_index", "mean_duration", "mean_track_speed",
              "n_tracks", "n_tracks_resp", "p_tracks_resp", 
              "mean_pref_index",  "conf.low", "conf.high", "p_value"), 
            as.numeric)

#display the total number of reps for each treatment
with(events, table(trt_color, CO2) )




#------------------------------------------------------
#-------- mixed model analysis of preference ----------
#------------------------------------------------------

#Create two column vector of treatment and control time
#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(round(tracks$trt_time*1000,0), 
            round(tracks$ctr_time*1000,0))

#model with individual intercepts for each treatment combination
pref.model <- glmmTMB(resp ~ 0 + trt_color:CO2 + (1|rep:trt_color:CO2), 
                      data = tracks, family = betabinomial(link = "logit"))
summary(pref.model)

#Plot residuals to see if model fits
pref.model_simres<-simulateResiduals(pref.model)
plot(pref.model_simres)


#Null Model
pref.null <- glmmTMB(resp ~ 0 + (1|rep:trt_color:CO2), 
                     data = tracks, family = betabinomial(link = "logit"))
summary(pref.null)

#Run Anova for treatment effect for females
anova(pref.model, pref.null)


#---- Analysis without clean air stimuli
tracks.sub<-
  tracks %>%
  #Filters to runs with at least some CO2
  filter(!CO2=="0%")

#Create two column vector of treatment and control time
#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(round(tracks.sub$trt_time*1000,0), 
            round(tracks.sub$ctr_time*1000,0))

#model with effect of treatment color and CO2
pref.sub.model <- glmmTMB(resp ~ 0 + trt_color + CO2 + trt_color:CO2 + (1|rep:trt_color:CO2), 
                          data = tracks.sub, family = betabinomial(link = "logit"))
summary(pref.sub.model)

#model with effect of treatment color and CO2 but no interaction
pref.sub.model.no.inter <- glmmTMB(resp ~ 0 + trt_color + CO2 + (1|rep:trt_color:CO2), 
                          data = tracks.sub, family = betabinomial(link = "logit"))
summary(pref.sub.model.no.inter)

#model with effect of treatment color
pref.sub.model.no.CO2 <- glmmTMB(resp ~ 0 + trt_color + (1|rep:trt_color:CO2), 
                          data = tracks.sub, family = betabinomial(link = "logit"))
summary(pref.sub.model.no.CO2)

#Null Model
pref.null <- glmmTMB(resp ~ (1|rep:trt_color:CO2), 
                     data = tracks.sub, family = betabinomial(link = "logit"))
summary(pref.null)

#Plot residuals to see if model fits
pref.sub.model_simres<-simulateResiduals(pref.sub.model)
plot(pref.sub.model_simres)

#Run Anova for effect of CO2
anova(pref.sub.model, pref.sub.model.no.CO2)

#Run Anova for overall model fit
anova(pref.sub.model.no.CO2, pref.null)




#---------------------------------------------------------------
#------------ Extract preference model predictions -------------
#---------------------------------------------------------------

#Extracts predictions for individual treatment:rep combinations
pref.point.pred<-
  predict_response(pref.model, terms=c("CO2", "trt_color","rep"), 
            type="random", ci_level=NA) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(CO2=x, trt_color=group, rep=facet, point_pred=predicted) %>%
  #relocate rep column
  #rearrange column order
  relocate(trt_color, .after = CO2) %>%
  relocate(rep, .after= trt_color) %>%
  #converts from probabilities to preference index
  mutate(point_pred=2*point_pred-1)


#Extracts model predictions, SE and confidence intervals
#There is an expected error about focal terms included as random effects
#but they are only included as an interaction with "rep"
pref.model.pred<-
  predict_response(pref.model, terms=c("CO2", "trt_color"),
                   type="fixed", ci_level = 0.95) %>%
  as.data.frame() %>%
  #rename variables
  rename(CO2=x, trt_color=group) %>%
  #rearrange column order
  relocate(trt_color, .after = CO2) %>%
  #extracts and adds the p-value from model summary
  mutate(p_value=summary(pref.model)$coefficients$cond[,4])



  
#------------------------------------------------------
#---- Collate display info for the preference graph ---
#------------------------------------------------------

#sets space between bars in the graph
sp_btw_bars<-1

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info.pref<-
  #expands these factors to all possible combinations
  expand.grid(CO2 = levels(tracks$CO2),
              trt_color = levels(tracks$trt_color)) %>%
  #reorders the rows
  arrange(CO2, trt_color) %>%
  #adds the other columns based on "trt_color"
  mutate(
    trt_intensity = case_when(
      trt_color == "all" ~ "0.00",
      TRUE ~ "1.00"),
    ctr_color =  "gray",
    ctr_intensity = case_when(
      trt_color == "gray" ~ "1.00",
      TRUE ~ "0.50")) %>%
  #----This section joins hex colors from "color.table"
  #adds a treatment color column for graph display
  mutate(trt.disp.col = color_lookup(trt_color, trt_intensity), 
         .after=trt_intensity) %>%
  #turns the color to NA to omit these treatment levels in the graph
  #as the number of recruited trajectories is too low
  mutate(trt.disp.col = case_when(
    CO2 == "0%" ~ NA,
    TRUE ~ trt.disp.col
  )) %>%
  #adds an alpha (40%) to the color
  mutate(trt.disp.col.alpha = color.add.alpha(trt.disp.col, 40), 
         .after=trt.disp.col) %>%
  #turns the color to NA to omit these treatment levels in the graph
  #as the number of recruited trajectories is too low
  mutate(trt.disp.col.alpha = case_when(
    CO2 == "0%" ~ NA,
    TRUE ~ trt.disp.col.alpha
  )) %>%
  #adds a control color column for graph display
  mutate(ctr.disp.col = color_lookup(ctr_color, ctr_intensity), 
         .after=ctr_intensity) %>%
  #adds an alpha (40%) to the color
  mutate(ctr.disp.col.alpha = color.add.alpha(ctr.disp.col, 40), 
         .after=ctr.disp.col) %>%
  #Set x positions for the graph
  #Sets the order of the color stimuli (gray, blue, green, red, black)
  mutate(x.pos=case_when(
    trt_color == "gray" ~ 1*sp_btw_bars,
    trt_color == "430"  ~ 2*sp_btw_bars,
    trt_color == "525"  ~ 3*sp_btw_bars,
    trt_color == "625"  ~ 4*sp_btw_bars,
    trt_color == "all"  ~ 5*sp_btw_bars)) %>%
  #Further orders based on "CO2"
  mutate(x.pos=case_when(
    CO2 == "0%"  ~ x.pos + 0*sp_btw_bars,
    CO2 == "1%"  ~ x.pos + 6.5*sp_btw_bars,
    CO2 == "5%"  ~ x.pos + 13.0*sp_btw_bars,
    CO2 == "10%" ~ x.pos + 19.5*sp_btw_bars)) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .before=CO2) %>%
  #sorts by "x.pos"
  arrange(x.pos) %>%
  #----this section joins model predictions
  left_join(pref.model.pred, by = join_by(CO2, trt_color)) %>%
  #drops unnecessary columns
  select(!c(std.error)) %>%
  #converts these variables from proportion to preference index
  mutate_at(c("predicted", "conf.low", "conf.high"), function(x) x*2-1) %>%
  #Determines significance stars for display
  mutate(sig_stars=case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) )) %>%
  #Sets the stars to NA in the treatment levels in the graph
  #as the number of recruited trajectories is too low
  mutate(sig_stars = case_when(
    CO2 == "0%" ~ NA,
    TRUE ~ sig_stars
  ))




#------------------------------------------------------
#--- mixed model analysis of experiment recruitment ---
#------------------------------------------------------

#Create two column vector of trajectories recruited to the stimuli and total trajectories
resp<-cbind(events$n_tracks_resp, 
            events$n_tracks)

#model with individual intercept estimation for all
recruit.model <- glmmTMB(resp ~ 0 + trt_color:CO2 + (1|rep),
                      data = events, family = betabinomial(link = "logit"))
summary(recruit.model)

#Plot residuals to see if model fits
recruit.model_simres<-simulateResiduals(recruit.model)
plot(recruit.model_simres)


#Null Model
recruit.null <- glmmTMB(resp ~ (1|rep), data = events, 
                        family = betabinomial(link = "logit"))
summary(recruit.null)

#likelhood ratio test for overall model fit
anova(recruit.model, recruit.null)


#no CO2 model
recruit.model.1 <- glmmTMB(resp ~ 0 + trt_color + (1|rep),
                      data = events, family = betabinomial(link = "logit"))
summary(recruit.model.1)

#likelhood ratio test to test significance of CO2 effect
anova(recruit.model, recruit.model.1)



#Creates a dummy variable combining "trt_color" and "CO2"
dummy<-with(events, paste0("CO2:",CO2,"_trt_color:",trt_color))
dummy<-factor(dummy, levels = 
                levels(as.factor(dummy))[c(5,1,2,3,4,10,6,7,8,9,15,11,12,13,14,20,16,17,18,19)])

#uses dummy variable to contrast all other treatment levels with no CO2, grey
recruit.model.sig <- glmmTMB(resp ~ dummy + (1|rep),
                          data = events, family = betabinomial(link = "logit"))
summary(recruit.model.sig)


#Calculate emmeans for all of the "trt_color:CO2:odor" treatment combinations
recruit.model.emm<-emmeans(recruit.model, specs = ~ trt_color:CO2)

#specify different combination of treatment   
#combinations to compare with vectors of 1s and 0s
CO2_00 <- c(rep(1,5), rep(0,5), rep(0,5), rep(0,5))
CO2_01 <- c(rep(0,5), rep(1,5), rep(0,5), rep(0,5))
CO2_05 <- c(rep(0,5), rep(0,5), rep(1,5), rep(0,5))
CO2_10 <- c(rep(0,5), rep(0,5), rep(0,5), rep(1,5))

#calculate the significance of the contrast between the different levels
contrast(recruit.model.emm, adjust = "sidak", method = 
           list("CO2_00 - CO2_01" = CO2_00 - CO2_01,
                "CO2_00 - CO2_05" = CO2_00 - CO2_05,
                "CO2_00 - CO2_10" = CO2_00 - CO2_10,
                "CO2_01 - CO2_05" = CO2_01 - CO2_05,
                "CO2_01 - CO2_10" = CO2_01 - CO2_10,
                "CO2_05 - CO2_10" = CO2_05 - CO2_10))



#Extracts predictions for individual treatment:rep combinations
recruit.point.pred<-
  predict_response(recruit.model, terms=c("CO2", "trt_color","rep"), 
                   type="random", ci_level=NA) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(CO2=x, trt_color=group, rep=facet, point_pred=predicted) %>%
  #rearrange column order
  relocate(trt_color, .after = CO2) %>%
  relocate(rep, .after= trt_color)


#Extracts model predictions, SE and confidence intervals
recruit.model.pred<-
  predict_response(recruit.model, terms=c("CO2", "trt_color"),
                   type="fixed", ci_level = 0.95) %>%
  as.data.frame() %>%
  #rename variables
  rename(CO2=x, trt_color=group) %>%
  #rearrange column order
  relocate(trt_color, .after = CO2) %>%
  #extracts and adds the p-value from model summary
  mutate(p_value=summary(recruit.model.sig)$coefficients$cond[,4])




#------------------------------------------------------
#--- Collate display info for the recruitment graph ---
#------------------------------------------------------

#sets space between bars in the graph
sp_btw_bars<-1

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info.recruit<-
  #uses events as input
  events %>%
  #----this section summarizes by treatment_factor
  group_by(CO2, trt_color) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups = "drop") %>%
  #----This section joins hex colors from "color.table"
  #adds a treatment color column for graph display
  mutate(trt.disp.col = color_lookup(trt_color, trt_intensity), 
         .after=trt_intensity) %>%
  #adds an alpha (70%) to the color
  mutate(trt.disp.col.alpha = color.add.alpha(trt.disp.col, 40), 
         .after=trt.disp.col) %>%
  #adds a control color column for graph display
  mutate(ctr.disp.col = color_lookup(ctr_color, ctr_intensity), 
         .after=ctr_intensity) %>%
  #adds an alpha (70%) to the color
  mutate(ctr.disp.col.alpha = color.add.alpha(ctr.disp.col, 40), 
         .after=ctr.disp.col) %>%
  #Set x positions for the graph
  #Sets the order of the color stimuli (gray, blue, green, red, black)
  mutate(x.pos=case_when(
    trt_color == "gray" ~ 1*sp_btw_bars,
    trt_color == "430"  ~ 2*sp_btw_bars,
    trt_color == "525"  ~ 3*sp_btw_bars,
    trt_color == "625"  ~ 4*sp_btw_bars,
    trt_color == "all"  ~ 5*sp_btw_bars)) %>%
  #Further orders based on "CO2"
  mutate(x.pos=case_when(
    CO2 == "0%"  ~ x.pos + 0*sp_btw_bars,
    CO2 == "1%"  ~ x.pos + 6.5*sp_btw_bars,
    CO2 == "5%"  ~ x.pos + 13.0*sp_btw_bars,
    CO2 == "10%" ~ x.pos + 19.5*sp_btw_bars)) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .before=CO2) %>%
  #sorts by "x.pos"
  arrange(x.pos) %>%
  #----this section joins model predictions
  left_join(recruit.model.pred, by = join_by(CO2, trt_color)) %>%
  #drops unnecessary columns
  select(!c(std.error)) %>%
  #Determines significance stars for display
  mutate(sig_stars=case_when(
    CO2=="0%" & trt_color=="gray" ~ as.character(NA),
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))
 


#------------------------------------------------------
#---- mixed model analysis of experiment activation ---
#------------------------------------------------------

#model with individual intercept estimation for all
activ.model <- glmmTMB(n_tracks ~ 0 + trt_color:CO2 + (1|rep),
                      data = events, family = nbinom1(link = "log"))
summary(activ.model)

#Plot residuals to see if model fits
activ.model_simres<-simulateResiduals(activ.model)
plot(activ.model_simres)


#no CO2 model
activ.model.1 <- glmmTMB(n_tracks ~ 0 + trt_color + (1|rep),
                         data = events, family = nbinom1(link = "log"))
summary(activ.model.1)

#likelhood ratio test to test significance of CO2 effect
anova(activ.model, activ.model.1)


#Null Model
activ.null <- glmmTMB(n_tracks ~ 0 + (1|rep),
                      data = events, family = nbinom1(link = "log"))
summary(activ.null)

#Run Anova for overall model fit
anova(activ.model, activ.null)


#Creates a dummy variable combining "trt_color" and "CO2"
dummy<-with(events, paste0("CO2:",CO2,"_trt_color:",trt_color))
dummy<-factor(dummy, levels = 
                levels(as.factor(dummy))[c(5,1,2,3,4,10,6,7,8,9,15,11,12,13,14,20,16,17,18,19)])

#uses dummy variable to contrast all other treatment levels with no CO2, grey
activ.model.sig <- glmmTMB(n_tracks ~ dummy + (1|rep),
                          data = events, family = nbinom1(link = "log"))
summary(activ.model.sig)


#Calculate emmeans for all of the "trt_color:CO2:odor" treatment combinations
activ.model.emm<-emmeans(activ.model, specs = ~ trt_color:CO2)

#specify different combination of treatment   
#combinations to compare with vectors of 1s and 0s
CO2_00 <- c(rep(1,5), rep(0,5), rep(0,5), rep(0,5))
CO2_01 <- c(rep(0,5), rep(1,5), rep(0,5), rep(0,5))
CO2_05 <- c(rep(0,5), rep(0,5), rep(1,5), rep(0,5))
CO2_10 <- c(rep(0,5), rep(0,5), rep(0,5), rep(1,5))

#calculate the significance of the contrast between the different levels
contrast(activ.model.emm, adjust = "sidak", method = 
           list("CO2_00 - CO2_01" = CO2_00 - CO2_01,
                "CO2_00 - CO2_05" = CO2_00 - CO2_05,
                "CO2_00 - CO2_10" = CO2_00 - CO2_10,
                "CO2_01 - CO2_05" = CO2_01 - CO2_05,
                "CO2_01 - CO2_10" = CO2_01 - CO2_10,
                "CO2_05 - CO2_10" = CO2_05 - CO2_10))



#Extracts predictions for individual treatment:rep combinations
activ.point.pred<-
  predict_response(activ.model, terms=c("CO2", "trt_color","rep"), 
                   type="random", ci_level=NA) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(CO2=x, trt_color=group, rep=facet, point_pred=predicted) %>%
  #rearrange column order
  relocate(trt_color, .after = CO2) %>%
  relocate(rep, .after= trt_color)


#Extracts model predictions, SE and confidence intervals
activ.model.pred<-
  predict_response(activ.model, terms=c("CO2", "trt_color"),
                   type="fixed", ci_level = 0.95) %>%
  as.data.frame() %>%
  #rename variables
  rename(CO2=x, trt_color=group) %>%
  #rearrange column order
  relocate(trt_color, .after = CO2) %>%
  #extracts and adds the p-value from model summary
  mutate(p_value=summary(activ.model.sig)$coefficients$cond[,4])




#------------------------------------------------------
#--- Collate display info for the activation graph ----
#------------------------------------------------------

#sets space between bars in the graph
sp_btw_bars<-1

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info.activ<-
  #uses events as input
  events %>%
  #----this section summarizes by treatment_factor
  group_by(CO2, trt_color) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups = "drop") %>%
  #----This section joins hex colors from "color.table"
  #adds a treatment color column for graph display
  mutate(trt.disp.col = color_lookup(trt_color, trt_intensity), 
         .after=trt_intensity) %>%
  #adds an alpha (70%) to the color
  mutate(trt.disp.col.alpha = color.add.alpha(trt.disp.col, 40), 
         .after=trt.disp.col) %>%
  #adds a control color column for graph display
  mutate(ctr.disp.col = color_lookup(ctr_color, ctr_intensity), 
         .after=ctr_intensity) %>%
  #adds an alpha (70%) to the color
  mutate(ctr.disp.col.alpha = color.add.alpha(ctr.disp.col, 40), 
         .after=ctr.disp.col) %>%
  #Set x positions for the graph
  #Sets the order of the color stimuli (gray, blue, green, red, black)
  mutate(x.pos=case_when(
    trt_color == "gray" ~ 1*sp_btw_bars,
    trt_color == "430"  ~ 2*sp_btw_bars,
    trt_color == "525"  ~ 3*sp_btw_bars,
    trt_color == "625"  ~ 4*sp_btw_bars,
    trt_color == "all"  ~ 5*sp_btw_bars)) %>%
  #Further orders based on "CO2"
  mutate(x.pos=case_when(
    CO2 == "0%"  ~ x.pos + 0*sp_btw_bars,
    CO2 == "1%"  ~ x.pos + 6.5*sp_btw_bars,
    CO2 == "5%"  ~ x.pos + 13.0*sp_btw_bars,
    CO2 == "10%" ~ x.pos + 19.5*sp_btw_bars)) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .before=CO2) %>%
  #sorts by "x.pos"
  arrange(x.pos) %>%
  #----this section joins model predictions
  left_join(activ.model.pred, by = join_by(CO2, trt_color)) %>%
  #drops unnecessary columns
  select(!c(std.error)) %>%
  #Determines significance stars for display
  mutate(sig_stars=case_when(
    CO2=="0%" & trt_color=="gray" ~ as.character(NA),
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))




#--------------------------------------------------------------
#----------------- code to generate graph ---------------------
#--------------------------------------------------------------

#------ this section has variables used in the graph call below

#Set figure width, margins and calculate plot sizes and figure height.
#All inputs in inches
f_width <- 5.5
f_height <- 8.116667
m_0 <- 0.0
m_bottom <- 0.3
m_bottom_outer <- 0.65
m_left <- 0.60
m_top <- 0.35
m_top_outer <- 0.35
m_right <- 0.1
plot_width <- f_width-m_left-m_right
f_height_no_margin <- f_height - 2*m_0 - 2*m_top - m_top_outer - m_bottom_outer
subplot_1_height <- (2/5) * f_height_no_margin
subplot_2_height <- (2/5) * f_height_no_margin
subplot_3_height <- (1/5) * f_height_no_margin
subplot_ypos<- c(1, 
                 1-(subplot_1_height + m_0 + m_top_outer)/f_height,
                 1-(subplot_1_height + 
                    subplot_2_height + 2*m_0 + 2*m_top_outer)/f_height, 
                 0)

#X and Y limits subplot 1
xlim.sp1<-c(0, max(display.info.recruit$x.pos)+1)
ylim.sp1<-c(0, 900)
yrange.sp1<-ylim.sp1[2]-ylim.sp1[1]
xaxis_ticks.sp1<-display.info.recruit$x.pos
yaxis_ticks.sp1<-seq(ylim.sp1[1], ylim.sp1[2], 200)
#size of the spacer between the bar and sigificance stars
sig_star_spacer.sp1 <- 0.04*yrange.sp1 #space between the bars and the star

#X and Y limits subplot 2
xlim.sp2<-c(0, max(display.info.recruit$x.pos)+1)
ylim.sp2<-c(0.00, 0.20)
yrange.sp2<-ylim.sp2[2]-ylim.sp2[1]
xaxis_ticks.sp2<-display.info.recruit$x.pos
yaxis_ticks.sp2<-seq(ylim.sp2[1], ylim.sp2[2], 0.05)
#size of the spacer between the bar and sigificance stars
sig_star_spacer.sp2 <- 0.04*yrange.sp2 #space between the bars and the star

#X and Y limits subplot 3
xlim.sp3<-c(0, max(display.info.pref$x.pos)+1)
ylim.sp3<-c(-0.20, 0.85)
yrange.sp3<-ylim.sp3[2]-ylim.sp3[1]
xaxis_ticks.sp3<-display.info.pref$x.pos
yaxis_ticks.sp3<-seq(0, ylim.sp3[2], 0.25)
#size of the spacer between the bar and sigificance stars
sig_star_spacer.sp3 <- 0.04*yrange.sp3 #space between the bars and the stars


axis_mgp <- c(3,0.45,0) #2nd value controls position of tick labels
tck <- -0.015 #Value to control tick mark length default is -0.01
CI_width <- 0.35 #half width of confidence interval rectangles
jitter <- 0.195 #half width of the jitter cloud
mean_lwd <- 2 #width of the mean line
sb.lwd <- 1.25 #width of the significance bars
#alternating ypos of sample size text
sample_size_ypos.sp1 <- c(ylim.sp1[2]+0.02*yrange.sp1, ylim.sp1[2]+0.08*yrange.sp1)
sample_size_ypos.sp2 <- c(ylim.sp2[2]+0.02*yrange.sp2, ylim.sp2[2]+0.08*yrange.sp2)
sample_size_cex <- 0.75 #size of sample size text
sample_size_col <- "gray70"
subplot_label_cex <- 1.5
subplot_label_xpos <- 2.8 #x position of subplot labels in lines (0.2"/line)
subplot_label_ypos.sp1 <- ylim.sp1[2] #y position of subplot labels in plot units
subplot_label_ypos.sp2 <- ylim.sp2[2] #y position of subplot labels in plot units
subplot_label_ypos.sp3 <- ylim.sp3[2] #y position of subplot labels in plot units

#Variables for circle size and positioning
#incorporating yrange and ylim[1] keeps the circle position constant 
#when ylims change
trt_circle_y <- ylim.sp3[1]-0.14*yrange.sp3 
ctr_circle_y <- trt_circle_y-0.13*yrange.sp3 #expressed relative to "trt_circle_y"
circle_text_x_spacer <- 1.3 #space between the circles and their labels
circle_size<-3 #character expansion factor default=1

#Variables for the positioning of the boxes around the circles
box_lwd <- 1.25
box_upper <- trt_circle_y+0.08*yrange.sp3
box_lower <- ctr_circle_y-0.08*yrange.sp3
box_horiz_sp<-0.5*sp_btw_bars
box_1_left <-display.info.pref$x.pos[1]-box_horiz_sp
box_1_right <- display.info.pref$x.pos[5]+box_horiz_sp
box_2_left <- display.info.pref$x.pos[6]-box_horiz_sp
box_2_right <- display.info.pref$x.pos[10]+box_horiz_sp
box_3_left <- display.info.pref$x.pos[11]-box_horiz_sp
box_3_right <- display.info.pref$x.pos[15]+box_horiz_sp
box_4_left <- display.info.pref$x.pos[16]-box_horiz_sp
box_4_right <- display.info.pref$x.pos[20]+box_horiz_sp
box_col <- "black"

#Calculates the x positions of the CO2 odor labels
CO2_odor_labels_x<-
  display.info.recruit %>%
  group_by(CO2) %>%
  summarise(x.pos = mean(x.pos), .groups = "drop") %>%
  pull(x.pos)
#Calculates a y-postion of the label
CO2_odor_labels_y<-box_lower-0.04*yrange.sp3
CO2_odor_labels_y<-rep(CO2_odor_labels_y, 
                       length.out=length(CO2_odor_labels_x))
CO2_levels<-levels(display.info.recruit$CO2)
#Specifies "CO2_odor_labels"
CO2_odor_labels<-c("clean air",
                   bquote( .(CO2_levels[2]) ~ CO[2] ),
                   bquote( .(CO2_levels[3]) ~ CO[2] ),
                   bquote( .(CO2_levels[4]) ~ CO[2] ) ) %>%
  as.expression()



#-------- Execute code from here to end to generate graph

# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure 2.pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)

#Export figure to a png with the width and height from above
#also sets font type and size
# png(file=paste0(output_dir,"figure 2.png"),
#     width=f_width, height=f_height, units="in", res=150,
#     family="sans", pointsize=8)

par(  #specifies the outer (whole plot) margin in inches
  omi=c(0,0,0,0),
  #allows drawing outside of the plotting area
  xpd = TRUE)



#---- Subplot 1

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_0,m_left,m_top_outer,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[2],subplot_ypos[1]),
  cex=1) #prevents font scaling in multi-panel figures

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim.sp1, ylim=ylim.sp1, xaxs="i", yaxs="i")

#axes
axis(side=2, las=1, at=yaxis_ticks.sp1, mgp=axis_mgp, tck=tck)
box(bty="l")

#axis label
mtext("activation (number of trajectories)", side=2, line=2.8)

#subtext label
mtext(bquote(bold("A")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos.sp1, las=1, cex=subplot_label_cex)

#significance bars
bar.1.y<-0.9*ylim.sp1[2]
sig.bar.spacer<-0.075*yrange.sp1
lines( x = c( rep(display.info.pref$x.pos[3], 2), rep(display.info.pref$x.pos[18], 2) ),
       y = c( -0.5*sig_star_spacer.sp1, 0, 0, -0.5*sig_star_spacer.sp1 ) + bar.1.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(3,18)] ),
      y = bar.1.y  + 0.5*sig_star_spacer.sp1,
      labels = "***", adj = c(0.5,0.5) )
bar.2.y<-bar.1.y - 1*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[8], 2), rep(display.info.pref$x.pos[18], 2) ),
       y = c( -0.5*sig_star_spacer.sp1, 0, 0, -0.5*sig_star_spacer.sp1 ) + bar.2.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(8,18)] ),
      y = bar.2.y  + 0.5*sig_star_spacer.sp1,
      labels = "***", adj = c(0.5,0.5) )
bar.3.y<-bar.1.y - 2*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[13], 2), rep(display.info.pref$x.pos[18], 2) ),
       y = c( -0.5*sig_star_spacer.sp1, 0, 0, -0.5*sig_star_spacer.sp1 ) + bar.3.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(13,18)] ),
      y = bar.3.y  + 0.5*sig_star_spacer.sp1,
      labels = "**", adj = c(0.5,0.5) )
bar.4.y<-bar.1.y - 3*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[3], 2), rep(display.info.pref$x.pos[13], 2) ),
       y = c( -0.5*sig_star_spacer.sp1, 0, 0, -0.5*sig_star_spacer.sp1 ) + bar.4.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(3,13)] ),
      y = bar.4.y  + 0.5*sig_star_spacer.sp1,
      labels = "***", adj = c(0.5,0.5) )
bar.5.y<-bar.1.y - 4*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[8], 2), rep(display.info.pref$x.pos[13], 2) ),
       y = c( -0.5*sig_star_spacer.sp1, 0, 0, -0.5*sig_star_spacer.sp1 ) + bar.5.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(8,13)] ),
      y = bar.5.y  + 0.75*sig_star_spacer.sp1,
      labels = "n.s.", adj = c(0.5,0.5) )
bar.6.y<-bar.1.y - 5*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[3], 2), rep(display.info.pref$x.pos[8], 2) ),
       y = c( -0.5*sig_star_spacer.sp1, 0, 0, -0.5*sig_star_spacer.sp1 ) + bar.6.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(3,8)] ),
      y = bar.6.y  + 0.5*sig_star_spacer.sp1,
      labels = "***", adj = c(0.5,0.5) )

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.activ$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width) 
          + display.info.activ$x.pos[i], 
          y=c(rep(display.info.activ$conf.low[i],2), 
              rep(display.info.activ$conf.high[i],2)),
          col=display.info.activ$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.activ$x.pos[i], 
        y=rep(display.info.activ$predicted[i],2),
        col=display.info.activ$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#individual points
stripchart(point_pred~trt_color+CO2, data=activ.point.pred, at=display.info.activ$x.pos,
           method = "jitter", jitter=jitter, vertical=TRUE, add=TRUE,
           pch=16, col=display.info.activ$trt.disp.col.alpha)

#significance stars
text(x=display.info.activ$x.pos,
     y=display.info.activ$conf.high+sig_star_spacer.sp1,
     labels=display.info.activ$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <- events %>% 
  mutate(odor = factor(odor, levels=c("air","CO2 01","CO2 05","CO2 10"))) %>% 
  group_by(odor, trt_color) %>% 
  summarise(n_tracks = sum(n_tracks), .groups="drop") %>%
  mutate(n_tracks = paste0( "(", n_tracks, ")" )) %>%
  pull(n_tracks)
text(x=display.info.activ$x.pos, y=rep(sample_size_ypos.sp1,length.out=length(display.info.activ$x.pos)), 
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)



#---- Subplot 2

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_0,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[3],subplot_ypos[2]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim.sp2, ylim=ylim.sp2, xaxs="i", yaxs="i")

#axes
axis(side=2, las=1, at=yaxis_ticks.sp2, mgp=axis_mgp, tck=tck)
box(bty="l")

#axis label
mtext("proportion recruited", side=2, line=2.8)

#subtext label
mtext(bquote(bold("B")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos.sp2, las=1, cex=subplot_label_cex)

#significance bars
bar.7.y<-0.9*ylim.sp2[2]
sig.bar.spacer<-0.075*yrange.sp2
lines( x = c( rep(display.info.pref$x.pos[3], 2), rep(display.info.pref$x.pos[18], 2) ),
       y = c( -0.5*sig_star_spacer.sp2, 0, 0, -0.5*sig_star_spacer.sp2 ) + bar.7.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(3,18)] ),
      y = bar.7.y  + 0.5*sig_star_spacer.sp2,
      labels = "***", adj = c(0.5,0.5) )
bar.8.y<-bar.7.y - 1*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[8], 2), rep(display.info.pref$x.pos[18], 2) ),
       y = c( -0.5*sig_star_spacer.sp2, 0, 0, -0.5*sig_star_spacer.sp2 ) + bar.8.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(8,18)] ),
      y = bar.8.y  + 0.5*sig_star_spacer.sp2,
      labels = "***", adj = c(0.5,0.5) )
bar.9.y<-bar.7.y - 2*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[13], 2), rep(display.info.pref$x.pos[18], 2) ),
       y = c( -0.5*sig_star_spacer.sp2, 0, 0, -0.5*sig_star_spacer.sp2 ) + bar.9.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(13,18)] ),
      y = bar.9.y  + 0.5*sig_star_spacer.sp2,
      labels = "*", adj = c(0.5,0.5) )
bar.10.y<-bar.7.y - 3*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[3], 2), rep(display.info.pref$x.pos[13], 2) ),
       y = c( -0.5*sig_star_spacer.sp2, 0, 0, -0.5*sig_star_spacer.sp2 ) + bar.10.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(3,13)] ),
      y = bar.10.y  + 0.5*sig_star_spacer.sp2,
      labels = "***", adj = c(0.5,0.5) )
bar.11.y<-bar.7.y - 4*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[8], 2), rep(display.info.pref$x.pos[13], 2) ),
       y = c( -0.5*sig_star_spacer.sp2, 0, 0, -0.5*sig_star_spacer.sp2 ) + bar.11.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(8,13)] ),
      y = bar.11.y  + 0.5*sig_star_spacer.sp2,
      labels = "***", adj = c(0.5,0.5) )
bar.12.y<-bar.7.y - 5*sig.bar.spacer
lines( x = c( rep(display.info.pref$x.pos[3], 2), rep(display.info.pref$x.pos[8], 2) ),
       y = c( -0.5*sig_star_spacer.sp2, 0, 0, -0.5*sig_star_spacer.sp2 ) + bar.12.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(3,8)] ),
      y = bar.12.y  + 0.5*sig_star_spacer.sp2,
      labels = "***", adj = c(0.5,0.5) )


#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.recruit$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width) 
          + display.info.recruit$x.pos[i], 
          y=c(rep(display.info.recruit$conf.low[i],2), 
              rep(display.info.recruit$conf.high[i],2)),
          col=display.info.recruit$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.recruit$x.pos[i], 
        y=rep(display.info.recruit$predicted[i],2),
        col=display.info.recruit$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#individual points
stripchart(point_pred~trt_color+CO2, data=recruit.point.pred, at=display.info.recruit$x.pos,
           method = "jitter", jitter=jitter, vertical=TRUE, add=TRUE,
           pch=16, col=display.info.recruit$trt.disp.col.alpha)

#significance stars
text(x=display.info.recruit$x.pos,
     y=display.info.recruit$conf.high+sig_star_spacer.sp2,
     labels=display.info.recruit$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels<-paste0("(", with(tracks, table(trt_color, CO2)), ")")
text(x=display.info.recruit$x.pos, y=rep(sample_size_ypos.sp2,length.out=length(display.info.recruit$x.pos)), 
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)



#---- Subplot 3

par(
  #prevents font scaling in multi-panel figures
  mai=c(m_bottom_outer,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[4],subplot_ypos[3]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim.sp3, ylim=ylim.sp3, xaxs="i", yaxs="i")

#axes
axis(side=2, las=1, at=yaxis_ticks.sp3, mgp=axis_mgp, tck=tck*1.5)
lines(x=rep(xlim.sp3[1],2), y=ylim.sp3)
lines(x=xlim.sp3, y=rep(0,2), col="grey70", lty=2, lwd=1.5, lend=1)

#axis label
mtext("preference index", side=2, line=2.8)

#subtext label
mtext(bquote(bold("C")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos.sp3, las=1, cex=subplot_label_cex)

#significance bar
bar.13.y<-0.80
lines( x = c( rep(display.info.pref$x.pos[8], 2), rep(display.info.pref$x.pos[13], 3), 
              rep(display.info.pref$x.pos[18], 2) ),
       y = c( -sig_star_spacer.sp3, 0, 0, -sig_star_spacer.sp3, 
              0, 0, -sig_star_spacer.sp3) + bar.13.y, 
       lwd = sb.lwd )
text( x = mean( display.info.pref$x.pos[c(8,18)] ),
      y = bar.13.y  + 2.25*sig_star_spacer.sp3,
      labels = "n.s.", adj = c(0.5,0.5) )

#insufficient recruitment
text( x = mean( display.info.pref$x.pos[1:5] ),
      y = 0.18, labels = "insufficient\nrecruitment", adj = c(0.5,0.5) )

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.pref$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width) 
          + display.info.pref$x.pos[i], 
          y=c(rep(display.info.pref$conf.low[i],2), 
              rep(display.info.pref$conf.high[i],2)),
          col=display.info.pref$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.pref$x.pos[i], 
        y=rep(display.info.pref$predicted[i],2),
        col=display.info.pref$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#individual points
stripchart(point_pred~trt_color+CO2, data=pref.point.pred, at=display.info.pref$x.pos,
           method = "jitter", jitter=jitter, vertical=TRUE, add=TRUE,
           pch=16, col=display.info.pref$trt.disp.col.alpha)

#significance stars
text(x=display.info.pref$x.pos, 
     y=display.info.pref$conf.high + 2*sig_star_spacer.sp3,
     labels=display.info.pref$sig_stars,
     adj=c(0.5,1))

#Circle text labels
text(x=min(display.info.recruit$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.recruit$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.recruit$x.pos, 
       y=rep(trt_circle_y, length(display.info.recruit$x.pos)), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.recruit$trt.disp.col,
       #circle outline color
       col=NA)
points(x=display.info.recruit$x.pos, 
       y=rep(ctr_circle_y, length(display.info.recruit$x.pos)), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.recruit$ctr.disp.col,
       #circle outline color
       col=NA)

#Boxes around the circles
polygon(x=c(box_1_left, box_1_right, box_1_right, box_1_left),
        y=c(box_upper, box_upper, box_lower, box_lower),
        col=NA, border=box_col, lwd=box_lwd)
polygon(x=c(box_2_left, box_2_right, box_2_right, box_2_left),
        y=c(box_upper, box_upper, box_lower, box_lower),
        col=NA, border=box_col, lwd=box_lwd)
polygon(x=c(box_3_left, box_3_right, box_3_right, box_3_left),
        y=c(box_upper, box_upper, box_lower, box_lower),
        col=NA, border=box_col, lwd=box_lwd)
polygon(x=c(box_4_left, box_4_right, box_4_right, box_4_left),
        y=c(box_upper, box_upper, box_lower, box_lower),
        col=NA, border=box_col, lwd=box_lwd)

#odor and CO2 box labels
text(x=CO2_odor_labels_x, y=CO2_odor_labels_y, labels=CO2_odor_labels, adj=c(0.5,1))


#Close Graph "Device"
dev.off()

