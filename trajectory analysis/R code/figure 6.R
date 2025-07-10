#Imports necessary packages
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(readr)
library(grImport2)
library(grid)
library(drc)
library(medrc)
library(multcomp)
#sets dplyr version of select to be called with "select()"
select<-dplyr::select

#Loads custom R functions
source("R code/wind tunnel functions.R")

#Sets the output directory
output_dir<-"../figures/"


#------------------------------------------------------
#-------- initial data import and processing ----------
#------------------------------------------------------

#Read in SVG images to use in graph
odor_source_dir<-"../odor source diagrams/"
CO2_SVG <- file.path(odor_source_dir, "CO2-cairo.svg") %>% readPicture()
tansy_SVG <- file.path(odor_source_dir, "tansy-cairo v2.svg") %>% readPicture()
foot_SVG <- file.path(odor_source_dir, "foot-cairo.svg") %>% readPicture()
alfalfa_SVG <- file.path(odor_source_dir, "plant infusion-cairo.svg") %>% readPicture()

#Read in relative intensities of each color ramp.
#While the amount of light from the LEDs scales linearly
#with the setting (IE 0.50 or 1.00) this is displayed on top
#of the background illumination which is not spectrally flat
#for this reason the same LED intensity settings modifying
#isoquantal stimuli will give different relative intensities
#when background illumination is taken into account.
rel_intensity<-read.csv("R code/relative intensities.csv", colClasses="character")

#Sets the sub directories to search for "resp_tracks.csv" track files
sub_dir_CO2<-"output files/color intensity ramps, CO2/"
sub_dir_CO2_VGrY<-"output files/color intensity ramps, CO2 (VGrY)/"
sub_dir_CO2_tansy<-"output files/color intensity ramps, CO2 + tansy/"
sub_dir_CO2_foot<-"output files/color intensity ramps, CO2 + foot/"
sub_dir_CO2_alfalfa<-"output files/color intensity ramps, CO2 + alfalfa/"

#Extracts track level data from each groups of runs
tracks<-track_extract(sub_dir_CO2)
tracks<-track_extract(sub_dir_CO2_VGrY, tracks)
tracks<-track_extract(sub_dir_CO2_tansy, tracks)
tracks<-track_extract(sub_dir_CO2_foot, tracks)
tracks<-track_extract(sub_dir_CO2_alfalfa, tracks)


#Creates a subset of data and processes it for analysis and display
tracks<-
  #starts with "tracks"
  tracks %>%
  #limits to intensity ramp stimuli
  filter(str_detect(stim_series, "intensity ramp")) %>%
  #The "trt_color" color of black vs grey stimuli are labeled according to their
  #stim series so the appropriate black vs gray sorts with the right color
  #intensity ramp
  mutate(trt_color = case_when(
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="430 nm intensity ramp" ~ "430",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="470 nm intensity ramp" ~ "470",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="525 nm intensity ramp" ~ "525",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="545 nm intensity ramp" ~ "545",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="625 nm intensity ramp" ~ "625",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="gray intensity ramp" ~ "gray",
    #all trt_colors are left unmodified
    TRUE ~ trt_color)) %>%
  #turns these variables into factors
  mutate(odor = factor(odor, levels=c("CO2","CO2 + tansy",
                                      "CO2 + foot","CO2 + alfalfa")),
         trt_color = as.factor(trt_color),
         trt_intensity = as.factor(trt_intensity),
         rep = as.factor(rep) ) %>%
  #turns these variables numeric
  mutate_at(c("pref_index", "duration", "trt_time", 
              "ctr_time", "stim_time"), as.numeric)

#Displays the total number of trajectories for each color, intensity & odors
table(list(tracks$trt_color, tracks$trt_intensity, tracks$odor))


#----------------------------------------------------------
#---------- mixed model analysis of preference ------------
#----------------------------------------------------------

#---- Gray intensity ramp only model
#Analysed separately as "trt_intensity" levels differ and unlike the other
#colors the gray ramp wasn't performed with different odors

#Filters "tracks" to gray intensity ramps only
tracks.gray<-
  tracks %>%
  filter(stim_series == "gray intensity ramp")

#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(round(tracks.gray$trt_time*1000,0), 
            round(tracks.gray$ctr_time*1000,0))

#The pref.model.gray is used to generate individual point and group predictions
#for the gray intensity ramp and to test if each group has a preference 
#differing from 50:50
pref.model.gray <- glmmTMB(resp ~ 0 + trt_intensity + (1|rep), 
                      data = tracks.gray, family = betabinomial(link = "logit"))

summary(pref.model.gray)

#Plot residuals to see if model fits
pref.model.gray_simres<-simulateResiduals(pref.model.gray)
plot(pref.model.gray_simres)


#Null Model
pref.null.gray <- glmmTMB(resp ~ 0 + (1|rep),  
                     data = tracks.gray, family = betabinomial(link = "logit"))
summary(pref.null.gray)

#Run Anova for overall model fit
anova(pref.model.gray, pref.null.gray)



#---- color intensity ramps  model

#Filters "tracks" to remove gray intensity ramps
tracks.sub<-
  tracks %>%
  filter(stim_series != "gray intensity ramp") %>%
  mutate(trt_color = droplevels(trt_color))


#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(round(tracks.sub$trt_time*1000,0), 
            round(tracks.sub$ctr_time*1000,0))

#The pref.model is used to generate individual point and group predictions
#and to test if each group has a preference differing from 50:50

#all stimuli model used for predictions many combinations of
#trt_color:trt_intensity:odor were not performed so model is rank deficient 
#and will give multiple warnings. This model and the ones below will 
#drop trt_color:trt_intensity:odor combinations with no data
pref.model <- glmmTMB(resp ~ 0 + trt_color:trt_intensity:odor + 
                        (1|trt_color:trt_intensity:rep), 
                      data = tracks.sub, family = betabinomial(link = "logit"))

#Processes model summary to remove dropped 
#trt_color:trt_intensity:odor combinations
pref.model.summary <- 
  #extracts coefficients
  summary(pref.model)$coefficients$cond %>% 
  as.data.frame() %>% 
  #filters the dropped trt_color:trt_intensity:odor levels
  filter(!is.na(Estimate)) %>% 
  #turns the row names into a column
  tibble::rownames_to_column() %>% 
  #splits the column to extract factor info
  separate_wider_delim(rowname, delim = ":",  
                       names = c("trt_color", "trt_intensity", "odor"), 
                       too_few="align_start") %>%
  #removes factor labels from the levels
  mutate( trt_color = gsub("trt_color", "", trt_color),
          trt_intensity = gsub("trt_intensity", "", trt_intensity),
          odor = gsub("odor", "", odor)) %>%
  #turns the new columns into factors using levels from "events
  mutate(trt_color = 
           factor(trt_color, levels=levels(tracks.sub$trt_color)),
         trt_intensity = 
           factor(trt_intensity, levels=levels(tracks.sub$trt_intensity)),
         odor = 
           factor(odor, levels=levels(tracks.sub$odor))) %>%
  rename(p_value = "Pr(>|z|)") %>%
  #sorts the table by odor color and then intensity
  arrange(odor, trt_color, trt_intensity)

#Plot residuals to see if model fits
pref.model_simres<-simulateResiduals(pref.model)
plot(pref.model_simres)
testOutliers(pref.model_simres, type="bootstrap")

#Null Model
pref.null <- glmmTMB(resp ~ 0 + (1|trt_color:trt_intensity:rep),  
                     data = tracks.sub, family = betabinomial(link = "logit"))
summary(pref.null)

#Run Anova for overall model fit
anova(pref.model, pref.null)


#---- Analysis testing individual effects

#Full Model
pref.model.1 <- glmmTMB(resp ~ trt_color + trt_intensity + odor + 
                          trt_color:trt_intensity + 
                          trt_color:odor + 
                          trt_intensity:odor + 
                          trt_color:trt_intensity:odor + 
                          (1|trt_color:trt_intensity:rep), 
                        data = tracks.sub, 
                        family = betabinomial(link = "logit"))
summary(pref.model.1)


#no 3-way interaction model
pref.model.2 <- glmmTMB(resp ~ trt_color + trt_intensity + odor + 
                          trt_color:trt_intensity + 
                          trt_color:odor + 
                          trt_intensity:odor +
                          (1|trt_color:trt_intensity:rep), 
                        data = tracks.sub, 
                        family = betabinomial(link = "logit"))
summary(pref.model.2)

#anova testing the 3-way interaction
anova(pref.model.1, pref.model.2)


#no trt_color:odor interaction 
pref.model.2.5 <- glmmTMB(resp ~ trt_color + trt_intensity + odor + 
                          trt_color:trt_intensity + 
                          trt_intensity:odor +
                          (1|trt_color:trt_intensity:rep), 
                        data = tracks.sub, 
                        family = betabinomial(link = "logit"))
summary(pref.model.2.5)

#anova testing the trt_color:odor interaction
anova(pref.model.2, pref.model.2.5)


#no odor model
pref.model.3 <- glmmTMB(resp ~ trt_color + trt_intensity +
                          trt_color:trt_intensity +
                          (1|trt_color:trt_intensity:rep), 
                        data = tracks.sub, 
                        family = betabinomial(link = "logit"))
summary(pref.model.3)

#anova testing the effect of odor
anova(pref.model.1, pref.model.3)

#anova testing the effect of odor (excluding the 3-way interaction)
anova(pref.model.2, pref.model.3)

#anova testing the effect of odor (excluding the 3-way & trt_color:odor interaction)
anova(pref.model.2.5, pref.model.3)


#no color model
pref.model.4 <- glmmTMB(resp ~ trt_intensity + odor +
                          trt_intensity:odor +
                          (1|trt_color:trt_intensity:rep), 
                        data = tracks.sub, 
                        family = betabinomial(link = "logit"))
summary(pref.model.4)

#anova testing the effect of color
anova(pref.model.1, pref.model.4)


#no intensity model
pref.model.5 <- glmmTMB(resp ~ trt_color + odor + 
                          trt_color:odor +
                          (1|trt_color:trt_intensity:rep), 
                        data = tracks.sub, 
                        family = betabinomial(link = "logit"))
summary(pref.model.5)

#anova testing the effect of intensity
anova(pref.model.1, pref.model.5)




#---------------------------------------------------------------
#------------ Extract preference model predictions -------------
#---------------------------------------------------------------

#---- Gray intensity ramp model predictions

#creates list of unique combinations of "trt_color",  
#"trt_intensity", "odor", and "rep"
newdata <- 
  tracks.gray %>% select(trt_color, trt_intensity, rep) %>% distinct()

#point predictions
point.pred.gray <-
  predict(pref.model.gray, re.form = NULL, type="response", newdata=newdata) %>%
  #combine with "newdata" frame
  cbind(newdata, .) %>%
  #rename prediction column
  rename(prop_pred = ".") %>%
  #converts from proportion to preference index
  mutate(point_pred=2*prop_pred-1) %>%
  #Adds relative intensity column from "rel_intensity" read in above
  left_join(rel_intensity, by = join_by(trt_intensity == int_setting, 
                                        trt_color == color)) %>%
  relocate(rel_intensity, .after=trt_intensity) %>%
  #correct data types
  mutate(trt_color = as.factor(trt_color),
         trt_intensity = as.factor(trt_intensity),
         rel_intensity = as.numeric(rel_intensity))

#Extracts overall model predictions, SE and confidence intervals
#There is an expected error about focal terms included as random effects
#but they are only included as an interaction with "rep"
model.pred.gray<-
  predict_response(pref.model.gray, terms=c("trt_intensity"), 
                   type="fixed", ci_level = 0.95) %>%
  #convert to dataframe
  as.data.frame() %>%
  #drops the group column
  select(!group) %>%
  #renames calculated columns to original terms
  rename(trt_intensity=x) %>%
  #drops SE column
  select(!std.error) %>%
  #relocate factor columns
  relocate(c(trt_intensity), .after=NULL) %>%
  #converts from proportion to preference index
  mutate(across(!trt_intensity, ~ 2 * .x - 1)) %>%
  #extracts and adds the p-value from model summary
  mutate(p_value=summary(pref.model.gray)$coefficients$cond[,4])


#---- Color intensity ramp model predictions

#creates list of unique combinations of "trt_color",  
#"trt_intensity", "odor", and "rep"
newdata <- 
  tracks.sub %>% select(trt_color, trt_intensity, odor, rep) %>% distinct()

#point predictions
point.pred <-
  predict(pref.model, re.form = NULL, type="response", newdata=newdata) %>%
  #combine with "newdata" frame
  cbind(newdata, .) %>%
  #rename prediction column
  rename(prop_pred = ".") %>%
  #converts from proportion to preference index
  mutate(point_pred=2*prop_pred-1) %>%
  #Adds relative intensity column from "rel_intensity" read in above
  left_join(rel_intensity, by = join_by(trt_intensity == int_setting, 
                                        trt_color == color)) %>%
  relocate(rel_intensity, .after=trt_intensity) %>%
  #correct data types
  mutate(trt_color = as.factor(trt_color),
         trt_intensity = as.factor(trt_intensity),
         rel_intensity = as.numeric(rel_intensity))




#--------------------------------------------------------------
#------------ sigmoid intensity model fitting -----------------
#--------------------------------------------------------------

#---- sigmoid model fit of gray intensity ramp

#model fitting 430, 470, 525, & 545 nm intensity ramps
sig.pref.model.gray <- 
  medrm( prop_pred ~ rel_intensity, 
         fct=LL.4(),
         start=c(3,0.4,0.7,1.5),
         random = b + e ~ 1|rep,
         data = point.pred.gray,
         control= nlmeControl(msMaxIter=200, pnlsTol=0.1))

summary(sig.pref.model.gray)

#check model fit
plot(sig.pref.model.gray)
plot(residuals(sig.pref.model.gray) ~ fitted(sig.pref.model.gray))

#writes gray coefficients to file
coef(sig.pref.model.gray) %>% as.data.frame() %>%
  rename(coef_value = 1) %>%
  write.csv(file="R Code/gray sigmoid coefficients.csv")

#The upper (d) and lower (c) asymptotes are taken from the gray model
#These variables are used as fixed terms the sigmoid fitting below.
#The coefficients must be rounded before including them in the medrm call
#not sure why.
coef_c <- coef(sig.pref.model.gray)[["c"]] %>% round(7)
coef_d <- coef(sig.pref.model.gray)[["d"]] %>% round(7)


#---- sigmoid model fit of all non red curves

#Creates and processed a subset of predicted points for the sigmoid modeling
point.pred.subset.no_reds<-
  point.pred %>%
  filter(trt_color!="625") %>%
  #Creates a combination factor to avoid missing color odor combinations
  mutate(odor_color = fct_cross(odor, trt_color, sep=" | ") )

#individual curves are fitted for each of these combined factor levels
#in that order. Starting values below follow this odor
levels(point.pred.subset.no_reds$odor_color)

#model fitting 430, 470, 525, & 545 nm intensity ramps
sig.pref.model.no_reds <- 
  medrm( prop_pred ~ rel_intensity, 
         fct=LL.4( fixed=c(NA, coef_c, coef_d, NA) ),
         start=list( fixed = c( rep(3, 11), rep(1.5, 11) ) ),
         curveid = b + e ~ odor_color,
         random = b + e ~ 1|rep,
         data = point.pred.subset.no_reds,
         control= nlmeControl(msMaxIter=100, pnlsTol=0.1))

summary(sig.pref.model.no_reds)

#check model fit
plot(sig.pref.model.no_reds)
plot(residuals(sig.pref.model.no_reds) ~ fitted(sig.pref.model.no_reds))


#---- sigmoid model fit of just red curves

point.pred.subset.reds<-
  point.pred %>%
  filter(trt_color=="625") %>%
  #filters out the anomalously low response to the black stimuli
  #in the CO2 625 nm intensity ramp series. These values caused fitting issues
  filter( !(trt_intensity=="0.00" & odor=="CO2") ) %>%
  #Creates a combination factor to avoid missing color odor combinations
  mutate(odor_color = fct_cross(odor, trt_color, sep=" | ") )

#individual curves are fitted for each of these combined factor levels
#in that order. Starting values below follow this odor
levels(point.pred.subset.reds$odor_color)

#model fitting 625 nm intensity ramps, fitted separately as the inclusion of 
#the 625 intensity ramps interfered with model fit of the other intensity ramps
#as the preference was relatively unaffected by the intensity of 625 nm LEDs
sig.pref.model.reds <- 
  medrm( prop_pred ~ rel_intensity, 
         fct=LL.4( fixed=c(NA, coef_c, coef_d, NA) ),
         start=list( fixed = c( rep(0.75, 3), rep(40, 3) ) ), 
         curveid = b + e ~ odor_color,
         random = b + e ~ 1|rep,
         data = point.pred.subset.reds,
         control= nlmeControl(msMaxIter=100, pnlsTol=0.3))

summary(sig.pref.model.reds)

#check model fit
plot(sig.pref.model.reds)
plot(residuals(sig.pref.model.reds) ~ fitted(sig.pref.model.reds))


#---- sigmoid model fit of all CO2 curves

point.pred.CO2.subset<-
  point.pred %>%
  filter( odor=="CO2" ) %>%
  filter( !(trt_intensity=="0.00" & trt_color=="625") )

sig.pref.model.CO2 <- 
  medrm( prop_pred ~ rel_intensity, 
         fct=LL.4( fixed=c(NA, coef_c, coef_d, NA) ),
         start=list( fixed = c(3, 3, 3, 3, 0.75,
                               1.5, 1.5, 1.5, 1.5, 40) ),
         curveid = b + e ~ trt_color,
         random = b + e ~ 1|rep,
         data = point.pred.CO2.subset,
         control= nlmeControl(msMaxIter=400, pnlsTol=0.1) )

summary(sig.pref.model.CO2)

#check model fit
plot(sig.pref.model.CO2)
plot(residuals(sig.pref.model.CO2) ~ fitted(sig.pref.model.CO2))




#--------------------------------------------------------------
#--------------- sigmoid ED50 comparisons ---------------------
#--------------------------------------------------------------

#sigmoid comparisons between different odors but with the same "trt_color"
sigmoid_comps.odors<-
  rbind(
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 430", "CO2 + foot | 430")) ,
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 430", "CO2 + alfalfa | 430")) ,
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 470", "CO2 + tansy | 470")) ,
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 525", "CO2 + tansy | 525")) ,
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 525", "CO2 + foot | 525")) ,
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 525", "CO2 + alfalfa | 525")) ,
    EDcomp(sig.pref.model.no_reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 545", "CO2 + foot | 545")) ,
    EDcomp(sig.pref.model.reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 625", "CO2 + tansy | 625")) ,
    EDcomp(sig.pref.model.reds, percVec=c(50,50), display = FALSE,
           compMatch = c("CO2 | 625", "CO2 + alfalfa | 625")) 
  ) %>% 
  as.data.frame() %>% rename(p_value = "p-value") %>%
  mutate(sig_stars=case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))

#sigmoid comparisons between different "trt_color"s but with CO2 alone
sigmoid_comps.CO2 <- 
  rbind(
    EDcomp(sig.pref.model.CO2, percVec=c(50,50),
           #reversed=TRUE so all values are compared to 525 nm
           display = FALSE, reverse = TRUE, 
           compMatch = c("525", "430")) ,
    EDcomp(sig.pref.model.CO2, percVec=c(50,50), 
           display = FALSE, reverse = TRUE,
           compMatch = c("525", "470")) ,
    EDcomp(sig.pref.model.CO2, percVec=c(50,50), 
           display = FALSE, reverse = FALSE,
           compMatch = c("525", "545")) ,
    EDcomp(sig.pref.model.CO2, percVec=c(50,50), 
           display = FALSE, reverse = FALSE,
           compMatch = c("525", "625"))
  ) %>% 
  as.data.frame() %>% rename(p_value = "p-value") %>%
  mutate(sig_stars=case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))




#--------------------------------------------------------------
#------------- sigmoid model line predictions -----------------
#--------------------------------------------------------------

#---- sets variables for line predictions below

#color max values for prediction below based on the max intensity rounded to
#the nearest tick mark division
max_rel_int.gray <- 0.5 * ceiling(max(point.pred.gray$rel_intensity) / 0.5 )
max_rel_int.reds <- 0.2 * ceiling(max(point.pred$rel_intensity) / 0.2 )
max_rel_int <- 1.8


#---- CO2 line predictions

line.CO2.gray <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.gray,0.01) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.gray, 
                             type="marginal", newdata=., se.fit=TRUE) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)

line.CO2.430 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 430", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.CO2.430.sub <- subset(line.CO2.430, rel_intensity<=max_rel_int)

line.CO2.430.CI<-
  #predicts 95% confidence interval for "effective doses" from 1 to 99% of
  #interval between the upper and lower asymptotes (coef_c, coef_d)
  ED(sig.pref.model.CO2, clevel="430", c(1:99), 
     interval="delta", level=0.95, display=FALSE) %>%
  as.data.frame() %>%
  #calculates the proportion values for all percentage values 
  #between "coef_c" and "coef_d"
  mutate(pref_index = coef_d - seq(0.01,0.99,0.01)*(coef_d-coef_c), .before=everything()) %>%
  #converts from proportion to preference index
  mutate(pref_index = 2*pref_index-1) %>%
  #filters to values within the xlims of the graph
  filter(Upper >= 1 & Lower <= max_rel_int) %>%
  #clips ends of confidence interval so nothing extends outside the 1 to max_rel_int range
  mutate( Lower = case_when( Lower < 1 ~ 1, TRUE ~ Lower),
          Upper = case_when( Upper > max_rel_int ~ max_rel_int, TRUE ~ Upper))

line.CO2.470 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 470", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.CO2.470.sub <- subset(line.CO2.470, rel_intensity<=max_rel_int)

line.CO2.470.CI<-
  #predicts 95% confidence interval for "effective doses" from 1 to 99% of
  #interval between the upper and lower asymptotes (coef_c, coef_d)
  ED(sig.pref.model.CO2, clevel="470", c(1:99), 
     interval="delta", level=0.95, display=FALSE) %>%
  as.data.frame() %>%
  #calculates the proportion values for all percentage values 
  #between "coef_c" and "coef_d"
  mutate(pref_index = coef_d - seq(0.01,0.99,0.01)*(coef_d-coef_c), .before=everything()) %>%
  #converts from proportion to preference index
  mutate(pref_index = 2*pref_index-1) %>%
  #filters to values within the xlims of the graph
  filter(Upper >= 1 & Lower <= max_rel_int) %>%
  #clips ends of confidence interval so nothing extends outside the 1 to max_rel_int range
  mutate( Lower = case_when( Lower < 1 ~ 1, TRUE ~ Lower),
          Upper = case_when( Upper > max_rel_int ~ max_rel_int, TRUE ~ Upper))

line.CO2.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.CO2.525.sub <- subset(line.CO2.525, rel_intensity<=max_rel_int)

#confidence interval polygon for "line.CO2.525"
line.CO2.525.CI<-
  #predicts 95% confidence interval for "effective doses" from 1 to 99% of
  #interval between the upper and lower asymptotes (coef_c, coef_d)
  ED(sig.pref.model.CO2, clevel="525", c(1:99), 
     interval="delta", level=0.95, display=FALSE) %>%
  as.data.frame() %>%
  #calculates the proportion values for all percentage values 
  #between "coef_c" and "coef_d"
  mutate(pref_index = coef_d - seq(0.01,0.99,0.01)*(coef_d-coef_c), .before=everything()) %>%
  #converts from proportion to preference index
  mutate(pref_index = 2*pref_index-1) %>%
  #filters to values within the xlims of the graph
  filter(Upper >= 1 & Lower <= max_rel_int.reds) %>%
  #clips ends of confidence interval so nothing extends outside the 1 to max_rel_int.reds range
  mutate( Lower = case_when( Lower < 1 ~ 1, TRUE ~ Lower),
          Upper = case_when( Upper > max_rel_int.reds ~ max_rel_int.reds, TRUE ~ Upper))
#subsets to the lower "max_rel_int"
line.CO2.525.CI.sub <- line.CO2.525.CI %>%
  filter(Lower <= max_rel_int) %>%
  #clips ends of confidence interval so nothing extends outside the 1 to max_rel_ints range
  mutate( Upper = case_when( Upper > max_rel_int ~ max_rel_int, TRUE ~ Upper))
  

line.CO2.545 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 545", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.CO2.545.sub <- subset(line.CO2.545, rel_intensity<=max_rel_int)

line.CO2.545.CI<-
  #predicts 95% confidence interval for "effective doses" from 1 to 99% of
  #interval between the upper and lower asymptotes (coef_c, coef_d)
  ED(sig.pref.model.CO2, clevel="545", c(1:99), 
     interval="delta", level=0.95, display=FALSE) %>%
  as.data.frame() %>%
  #calculates the proportion values for all percentage values 
  #between "coef_c" and "coef_d"
  mutate(pref_index = coef_d - seq(0.01,0.99,0.01)*(coef_d-coef_c), .before=everything()) %>%
  #converts from proportion to preference index
  mutate(pref_index = 2*pref_index-1) %>%
  #filters to values within the xlims of the graph
  filter(Upper >= 1 & Lower <= max_rel_int) %>%
  #clips ends of confidence interval so nothing extends outside the 1 to max_rel_int range
  mutate( Lower = case_when( Lower < 1 ~ 1, TRUE ~ Lower),
          Upper = case_when( Upper > max_rel_int ~ max_rel_int, TRUE ~ Upper))

line.CO2.625 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 625", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)

line.CO2.625.CI<-
  #predicts 95% confidence interval for "effective doses" from 1 to 99% of
  #interval between the upper and lower asymptotes (coef_c, coef_d)
  ED(sig.pref.model.CO2, clevel="625", seq(0.1,99,0.1), 
     interval="delta", level=0.95, display=FALSE) %>%
  as.data.frame() %>%
  #calculates the proportion values for all percentage values 
  #between "coef_c" and "coef_d"
  mutate(pref_index = coef_d - seq(0.001,0.99,0.001)*(coef_d-coef_c), .before=everything()) %>%
  #converts from proportion to preference index
  mutate(pref_index = 2*pref_index-1) %>%
  #filters to values within the xlims of the graph
  filter(Upper >= 1 & Estimate <= max_rel_int.reds) %>%
  #clips ends of confidence interval so nothing extends outside the 1 to max_rel_int.reds range
  mutate( Upper = case_when( Upper > max_rel_int.reds ~ max_rel_int.reds, TRUE ~ Upper))


#---- tansy line predictions

line.tansy.470 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + tansy | 470", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.tansy.470.sub <- subset(line.tansy.470, rel_intensity<=max_rel_int)

line.tansy.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + tansy | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.tansy.525.sub <- subset(line.tansy.525, rel_intensity<=max_rel_int)

line.tansy.625 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + tansy | 625", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)


#---- foot line predictions

line.foot.430 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + foot | 430", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.foot.430.sub <- subset(line.foot.430, rel_intensity<=max_rel_int)

line.foot.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + foot | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.foot.525.sub <- subset(line.foot.525, rel_intensity<=max_rel_int)

line.foot.545 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + foot | 545", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.foot.545.sub <- subset(line.foot.545, rel_intensity<=max_rel_int)


#---- alfalfa prediction lines

line.alfalfa.430 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + alfalfa | 430", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.alfalfa.430.sub <- subset(line.alfalfa.430, rel_intensity<=max_rel_int)

line.alfalfa.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + alfalfa | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)
line.alfalfa.525.sub <- subset(line.alfalfa.525, rel_intensity<=max_rel_int)

line.alfalfa.625 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_rel_int.reds,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + alfalfa | 625", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1)




#------------------------------------------------------
#-------- Collate display info for the graph ----------
#------------------------------------------------------
#sets space between bars in the graph
sp_btw_bars<-0.075

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info.gray<-
  #uses tracks.gray as input
  tracks.gray %>%
  #----this section summarizes by trt_color and odor
  group_by(trt_color, trt_intensity, odor) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_color=first(trt_color),
            trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups="drop") %>%
  #moves "odor" to the first column 
  relocate(odor, .after=NULL) %>%
  #----This section joins hex colors from "color.table"
  #adds a treatment color column for graph display
  mutate(trt.disp.col = color_lookup(trt_color, trt_intensity), 
         .after=trt_intensity) %>%
  #adds an alpha (40%) to the color
  mutate(trt.disp.col.alpha = color.add.alpha(trt.disp.col, 40), 
         .after=trt.disp.col) %>%
  #adds a control color column for graph display
  mutate(ctr.disp.col = color_lookup(ctr_color, ctr_intensity), 
         .after=ctr_intensity) %>%
  #adds an alpha (40%) to the color
  mutate(ctr.disp.col.alpha = color.add.alpha(ctr.disp.col, 40), 
         .after=ctr.disp.col) %>%
  #Adds relative intensity column from "rel_intensity" read in above
  left_join(rel_intensity, by = join_by(trt_intensity == int_setting, 
                                        trt_color == color)) %>%
  mutate(rel_intensity = log10(as.numeric(rel_intensity))) %>%
  #moves new intensity column after "trt_intensity"
  relocate(rel_intensity, .after=trt_intensity) %>%
  #----This section determines x positions to display data
  mutate(x.pos=rel_intensity) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=trt_intensity) %>%
  #arrange by "odor" and then "x.pos"
  arrange(odor, x.pos) %>%
  #----this section joins model predictions
  left_join(model.pred.gray, by="trt_intensity") %>%
  #Determines significance stars for display
  mutate(sig_stars=case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))


#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info<-
  #uses tracks.gray as input
  tracks %>%
  #filters out gray intensity ramp
  filter(trt_color!="gray") %>%
  #----this section summarizes by trt_color and odor
  group_by(trt_color, trt_intensity) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_color=first(trt_color),
            trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups="drop") %>%
  #----This section joins hex colors from "color.table"
  #adds a treatment color column for graph display
  mutate(trt.disp.col = color_lookup(trt_color, trt_intensity), 
         .after=trt_intensity) %>%
  #adds a control color column for graph display
  mutate(ctr.disp.col = color_lookup(ctr_color, ctr_intensity), 
         .after=ctr_intensity) %>%
  #Adds relative intensity column from "rel_intensity" read in above
  left_join(rel_intensity, by = join_by(trt_intensity == int_setting, 
                                        trt_color == color)) %>%
  mutate(rel_intensity = log10(as.numeric(rel_intensity))) %>%
  #moves new intensity column after "trt_intensity"
  relocate(rel_intensity, .after=trt_intensity) %>%
  #----This section determines x positions to display data
  mutate(x.pos=rel_intensity) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=trt_intensity)

#subsets for each wavelength
display.info.430<-subset(display.info, trt_color=="430")
display.info.470<-subset(display.info, trt_color=="470")
display.info.525<-subset(display.info, trt_color=="525")
display.info.545<-subset(display.info, trt_color=="545")
display.info.625<-subset(display.info, trt_color=="625")




#--------------------------------------------------------------
#----------------- code to generate graph ---------------------
#--------------------------------------------------------------

#------ this section has variables used in the graph call below

#Set figure width, margins and calculate plot sizes and figure height.
#All inputs in inches
f_width <- 7
m_bottom_row1 <- 0.7
m_bottom_row2 <- 0.7
m_bottom_row3 <- 0.9
m_left <- 0.03
m_left_outer <- 0.5
m_top <- 0.1
m_top_outer <- 0.3
m_right <- 0.12
m_right_outer <- 0.1
subplot_height <- 1.81
subplot_row1_height <- subplot_height + m_top_outer + m_bottom_row1
subplot_row2_height <- subplot_height + m_top + m_bottom_row2
subplot_row3_height <- subplot_height + m_top + m_bottom_row3
f_height <- sum(subplot_row1_height, subplot_row2_height, subplot_row3_height)
subplot_width_2col <- (f_width - m_left_outer-m_right - m_left-m_right_outer)/2
subplot_width_4col <- (f_width - m_left_outer-m_right - 2*(m_left+m_right)
                       - m_left-m_right_outer) / 4
subplot_col1_width <- subplot_width_4col + m_left_outer + m_right
subplot_col2_width <- subplot_width_4col + m_left + m_right
subplot_col3_width <- subplot_width_4col + m_left + m_right
subplot_col4_width <- subplot_width_4col + m_left + m_right_outer

#calculates subplot x positions
subplot_xpos_col <- c(0, 
                       (subplot_col1_width) / f_width, 
                       (subplot_col1_width + 
                              subplot_col2_width) / f_width,
                       (subplot_col1_width + 
                              subplot_col2_width + 
                              subplot_col3_width ) / f_width,
                       1)

#calculates subplot y positions
subplot_ypos <- c(1, 
                  1 - (subplot_row1_height) / f_height, 
                  1 - (subplot_row1_height + subplot_row2_height) / f_height, 
                  0)

#x limits and tick labels
xlim.gray<-c(-0.4*sp_btw_bars, log10(max_rel_int.gray) )
xlim.reds<-c(-0.36*sp_btw_bars, log10(max_rel_int.reds) )
xlim<-c(-0.24*sp_btw_bars, log10(max_rel_int) )
xrange.gray<-xlim.gray[2]-xlim.gray[1]
xrange.reds<-xlim.reds[2]-xlim.reds[1]
xrange<-xlim[2]-xlim[1]
#calculates tick labels and positions from "rel_intensity"
xaxis_tick_labels.gray<-
  seq(1, 10^xlim.gray[2], 0.5) %>% format(digits = 2)
xaxis_ticks.gray <- xaxis_tick_labels.gray %>% as.numeric() %>% log10()
xaxis_tick_labels.reds<-
  seq(1, 10^xlim.reds[2], 0.2) %>% format(digits = 2)
xaxis_ticks.reds <- xaxis_tick_labels.reds %>% as.numeric() %>% log10() 
xaxis_tick_labels<-
  seq(1, 10^xlim[2], 0.2) %>% format(digits = 2)
xaxis_ticks <- xaxis_tick_labels %>% as.numeric() %>% log10() 

#Y limits and tick labels
ylim.gray<-c(-0.35, 0.60)
yrange.gray<-ylim.gray[2]-ylim.gray[1]
yaxis_ticks.gray<-seq(-0.25, ylim.gray[2], 0.25)
ylim<-c(-0.15, 0.50)
yrange<-ylim[2]-ylim[1]
yaxis_ticks<-seq(0, ylim[2], 0.25)

axis_mgp.x <- c(3,0.05,0) #2nd value controls position of tick labels
axis_mgp.y <- c(3,0.35,0) #2nd value controls position of tick labels
tck <- -0.015 #Value to control tick mark length default is -0.01
#number of treatments
n_trts <- length(levels(point.pred.gray$trt_intensity))
CI_width <- 0.18*sp_btw_bars #half width of confidence interval rectangles in plot units
jitter <- 0.557*CI_width #half width of the jitter cloud
mean_lwd <- 2 #width of the mean line
#alternating ypos of sample size text
sample_size_ypos <- c(ylim.gray[2]+0.10*yrange, ylim.gray[2]+0.02*yrange)
sample_size_cex <- 5/8 #size of sample size text
subplot_label_cex <- 12/8 #size of subplot levels
subplot_label_xpos1 <- 0.6 #x position of subplot labels in lines (0.2"/line)
subplot_label_xpos2 <- -1.6 #x position of subplot labels in lines (0.2"/line)
subplot_label_ypos1 <- ylim.gray[2] #y position of subplot labels in plot units
subplot_label_ypos2 <- ylim[2] #y position of subplot labels in plot units
cex.axis <- 6/8 #tick label font size
sample_size_col <- "gray70"
#space between the bars and the stars
sig_star_spacer <- 0.05*yrange

#line styles
np_line_col <- "gray70"
np_line.lty <- 2 #no preference line type
np_line.lwd <- 1.5 #no preference line width
np_line.lend <- 1 #line end type "butt"

#gray divider line positions
div.line.red <- rep(0 - 0.055*xrange.reds, 2)
div.line.1col <- rep(0 - 0.080*xrange, 2)
div.line.2col <- rep(0 - 0.055*xrange, 2)

#Variables for circle size and positioning
#incorporating yrange and ylim[1] keeps the circle position constant 
#when ylims change
trt_circle_row1_y <- ylim.gray[1]-0.150*yrange.gray 
ctr_circle_row1_y <- trt_circle_row1_y-0.090*yrange.gray #expressed relative to "trt_circle_row1_y"
trt_circle_row2_y <- ylim[1]-0.150*yrange 
ctr_circle_row2_y <- trt_circle_row2_y-0.090*yrange #expressed relative to "trt_circle_row2_y"
trt_circle_row3_y <- ylim[1]-0.250*yrange 
ctr_circle_row3_y <- trt_circle_row3_y-0.090*yrange #expressed relative to "trt_circle_row3_y"
circle_text_x_spacer <- 0.05 #space between the circles and their labels in inches
circle_size<-3 #character expansion factor default=1
trt_circle_outline <- NA
ctr_circle_outline <- NA

#Variables for the positioning of the boxes around the circles
box_lwd <- 1.25
int_box_row1_upper <- trt_circle_row1_y + 0.050*yrange.gray
int_box_row1_lower <- ctr_circle_row1_y - 0.050*yrange.gray
int_box_row2_upper <- trt_circle_row2_y + 0.050*yrange
int_box_row2_lower <- ctr_circle_row2_y - 0.050*yrange
int_box_row3_upper <- trt_circle_row3_y + 0.050*yrange
int_box_row3_lower <- ctr_circle_row3_y - 0.050*yrange
box_label_row1_y <- int_box_row1_lower-0.02*yrange.gray
box_label_row2_y <- int_box_row2_lower-0.02*yrange
box_label_row3_y <- int_box_row3_lower-0.02*yrange
int_box_gray_left <- xlim.gray[1]
int_box_gray_right <- xlim.gray[2]
int_box_430_left <- xlim[1]
int_box_430_right <- xlim[2]
int_box_470_left <- div.line.1col[1]
int_box_470_right <- xlim[2] + 0.06*xrange
int_box_525_left <- div.line.2col[1] + 0.02*xrange
int_box_525_right <- xlim[2]
int_box_545_left <- xlim[1] + 0.027*xrange
int_box_545_right <- xlim[2]
int_box_625_left <- div.line.red[1] + 0.02*xrange.reds
int_box_625_right <- xlim.reds[2] + 0.01*xrange.reds
#box colors
box_col_gray <-  
  subset(display.info.gray, trt_intensity=="0.60", select=trt.disp.col) %>% 
  unique() %>%  as.character()
box_col_430 <-
  subset(display.info.430, trt_intensity=="1.00", select=trt.disp.col) %>% 
  unique() %>%  as.character()
box_col_470 <-
  subset(display.info.470, trt_intensity=="1.00", select=trt.disp.col) %>% 
  unique() %>%  as.character()
box_col_525 <-
  subset(display.info.525, trt_intensity=="1.00", select=trt.disp.col) %>%
  unique() %>%  as.character()
box_col_545 <-
  subset(display.info.545, trt_intensity=="1.00", select=trt.disp.col) %>%
  unique() %>%  as.character()
box_col_625 <-
  subset(display.info.625, trt_intensity=="1.00", select=trt.disp.col) %>%
  unique() %>%  as.character()

#sigmoid line styling and colors
sig.lwd <- 2 #width of the sigmoid lines
CO2.lwd <- 2.5 #width of the sigmoid lines
CO2.lty     <- "32" #CO2 line type
od_label.cex1 <- 8/8 #odor label size subplots a,b
od_label.cex2 <- 6/8 #odor label size subplots c-g
tansy_col   <- rgb(255, 240,  10, alpha=128, maxColorValue = 255)
foot_col    <- rgb(200, 145, 105, alpha=128, maxColorValue = 255)
alfalfa_col <- rgb(100, 155,  60, alpha=128, maxColorValue = 255)
line_col_430 <- color.add.alpha(box_col_430, 50)
line_col_470 <- color.add.alpha(box_col_470, 50)
line_col_525 <- color.add.alpha(box_col_525, 50)
line_col_545 <- color.add.alpha(box_col_545, 50)
line_col_625 <- color.add.alpha(box_col_625, 50)
line_col_430_CI <- color.add.alpha(box_col_430, 15)
line_col_470_CI <- color.add.alpha(box_col_470, 15)
line_col_525_CI <- color.add.alpha(box_col_525, 15)
line_col_545_CI <- color.add.alpha(box_col_545, 15)
line_col_625_CI <- color.add.alpha(box_col_625, 7)



#-------- Execute code from here to end to generate graph

#Export figure to a pdf with the width and height from above
#also sets font type and size
pdf(file=paste0(output_dir,"figure 6 (org).pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)


par(
  #specifies the outer (whole plot) margin in inches
  omi=c(0,0,0,0),
  #allows drawing outside of the plotting area
  xpd = TRUE)



#---- Subplot a

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row1,m_left_outer,m_top_outer,m_right),
  #specifies subplot position
  fig=c(subplot_xpos_col[1], subplot_xpos_col[3],
        subplot_ypos[2],subplot_ypos[1]),
  cex=1) #prevents font scaling in multi-panel figures

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim.gray, ylim=ylim.gray, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("A")), side=2, line=subplot_label_xpos1, 
      at=subplot_label_ypos1, las=1, cex=subplot_label_cex)

#axis
axis(side=1, las=1, at=xaxis_ticks.gray, labels=xaxis_tick_labels.gray,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks.gray, mgp=axis_mgp.y, tck=tck, cex.axis=cex.axis)
box(bty="l")

#no preference lines
lines(x=xlim.gray, y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.gray$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width)
          + display.info.gray$x.pos[i],
          y=c(rep(display.info.gray$conf.low[i],2),
              rep(display.info.gray$conf.high[i],2)),
          col=display.info.gray$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.gray$x.pos[i],
        y=rep(display.info.gray$predicted[i],2),
        col=display.info.gray$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#set seed for jitter
set.seed(123456)

#stripchart
stripchart(point_pred~trt_intensity*trt_color, data=point.pred.gray, 
           at=display.info.gray$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.gray$trt.disp.col.alpha)

#draws sigmoid model lines
lines(x = log10(line.CO2.gray$rel_intensity), y = line.CO2.gray$pred_line, 
      col=box_col_gray, lwd=sig.lwd)

#significance stars
text(x=display.info.gray$x.pos,
     y=display.info.gray$conf.high+sig_star_spacer,
     labels=display.info.gray$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <-
  tracks %>% filter(odor == "CO2" & trt_color== "gray") %>% 
  group_by(trt_color, trt_intensity, odor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = length(trt_color), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.gray$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #odor labels
# text(x=log10(1), y=-0.25, adj=c(0,0.5), cex=od_label.cex1, 
#      labels=bquote(bold(CO[2]~"alone")))

#Circle text labels
text(x=xlim.gray[1]-circle_text_x_spacer*(xrange.gray/subplot_width_2col), 
     y=trt_circle_row1_y, adj=c(1,0.5), labels="test")
text(x=xlim.gray[1]-circle_text_x_spacer*(xrange.gray/subplot_width_2col), 
     y=ctr_circle_row1_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.gray$x.pos, y=rep(trt_circle_row1_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.gray$trt.disp.col,
       col=NA)
points(x=display.info.gray$x.pos, y=rep(ctr_circle_row1_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.gray$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_gray_left, int_box_gray_right, int_box_gray_right, int_box_gray_left),
        y=c(int_box_row1_upper, int_box_row1_upper, int_box_row1_lower, int_box_row1_lower),
        col=NA, border=box_col_gray, lwd=box_lwd)

#box labels
text(x=int_box_gray_right, y=box_label_row1_y, adj=c(1,1),
     labels=c("gray intensity ramp"), col=box_col_gray)



#---- Subplot b

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row1,m_left,m_top_outer,m_right_outer),
  #specifies subplot position
  fig=c(subplot_xpos_col[3], subplot_xpos_col[5],
        subplot_ypos[2],subplot_ypos[1]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim.reds, ylim=ylim.gray, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("B")), side=2, line=subplot_label_xpos2,
      at=subplot_label_ypos1, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks.reds, labels=xaxis_tick_labels.reds,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0)

#no preference lines
lines(x=c(0,xlim.reds[2]), y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=div.line.red, y=ylim.gray, col="grey90", lwd=1.5)

#confidence interval around 525 line
polygon(x=log10( c(line.CO2.525.CI$Upper, rev(line.CO2.525.CI$Lower) ) ),
        y=c(line.CO2.525.CI$pref_index, rev(line.CO2.525.CI$pref_index)),
        col=line_col_525_CI, border=NA)

#draws sigmoid model lines
lines(x = log10(line.CO2.625$rel_intensity), y = line.CO2.625$pred_line,
      col=line_col_625, lwd=CO2.lwd, lty=CO2.lty)
lines(x = log10(line.CO2.430$rel_intensity), y = line.CO2.430$pred_line,
      col=line_col_430, lwd=CO2.lwd, lty=CO2.lty)
lines(x = log10(line.CO2.545$rel_intensity), y = line.CO2.545$pred_line,
      col=line_col_545, lwd=CO2.lwd, lty=CO2.lty)
lines(x = log10(line.CO2.470$rel_intensity), y = line.CO2.470$pred_line,
      col=line_col_470, lwd=CO2.lwd, lty=CO2.lty)
lines(x = log10(line.CO2.525$rel_intensity), y = line.CO2.525$pred_line,
      col=line_col_525, lwd=CO2.lwd, lty=CO2.lty)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #sigmoid line labels
# text(x=log10(2.2), y=0.48, adj=c(0.5,0.5), srt=0,
#      labels="621 nm", col=box_col_625)
# text(x=log10(2.2), y=0.085, adj=c(0.5,0.5), srt=-21,
#      labels="435 nm", col=box_col_430)
# text(x=log10(1.385), y=0.27, adj=c(0.5,0.5), srt=-27,
#      labels="552 nm", col=box_col_545)
# text(x=log10(1.585), y=0.138, adj=c(0.5,0.5), srt=-26,
#      labels="478 nm", col=box_col_470)
# text(x=log10(2.2), y=-0.24, adj=c(0.5,0.5), srt=0,
#      labels="527 nm", col=box_col_525)
# 
# #sigmoid line significance stars
# text(x=log10(2.3), y=0.10, adj=c(0.5,0.5), srt=-21, 
#      labels=sigmoid_comps.CO2$sig_stars[1])
# text(x=log10(1.65), y=0.13, adj=c(0.5,0.5), srt=-26, 
#      labels=sigmoid_comps.CO2$sig_stars[2])
# text(x=log10(1.45), y=0.26, adj=c(0.5,0.5), srt=-27, 
#      labels=sigmoid_comps.CO2$sig_stars[3])
# text(x=log10(2.3), y=0.52, adj=c(0.5,0.5), srt=0, 
#      labels=sigmoid_comps.CO2$sig_stars[4])
# 
# #odor labels
# text(x=log10(1), y=-0.25, adj=c(0,0.5), cex=od_label.cex1, 
#      labels=bquote(bold(CO[2]~"alone")))

#Treatment and control circles
points(x=display.info.625$x.pos, y=rep(trt_circle_row1_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.625$trt.disp.col,
       col=NA)
points(x=display.info.625$x.pos, y=rep(ctr_circle_row1_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.625$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_625_left, int_box_625_right, int_box_625_right, int_box_625_left),
        y=c(int_box_row1_upper, int_box_row1_upper, int_box_row1_lower, int_box_row1_lower),
        col=NA, border=box_col_gray, lwd=box_lwd)

#box labels
text(x=int_box_625_right, y=box_label_row1_y, adj=c(1,1),
     labels=c("color intensity ramps"), col=box_col_gray)



#---- Subplot c

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row2,m_left_outer,m_top,m_right),
  #specifies subplot position
  fig=c(subplot_xpos_col[1], subplot_xpos_col[2],
        subplot_ypos[3],subplot_ypos[2]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("C")), side=2, line=subplot_label_xpos2,
      at=subplot_label_ypos2, las=1, cex=subplot_label_cex)

#axis
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp.y, tck=tck, cex.axis=cex.axis)
box(bty="l")

#no preference lines
lines(x=xlim, y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#confidence interval around 430 line
polygon(x=log10( c(line.CO2.430.CI$Upper, rev(line.CO2.430.CI$Lower) ) ),
        y=c(line.CO2.430.CI$pref_index, rev(line.CO2.430.CI$pref_index)),
        col=line_col_430_CI, border=NA)

#draws sigmoid model lines
lines(x = log10(line.foot.430.sub$rel_intensity), y = line.foot.430.sub$pred_line,
      col=foot_col, lwd=sig.lwd)
lines(x = log10(line.alfalfa.430.sub$rel_intensity), y = line.alfalfa.430.sub$pred_line,
      col=alfalfa_col, lwd=sig.lwd)
lines(x = log10(line.CO2.430.sub$rel_intensity), y = line.CO2.430.sub$pred_line,
      col=line_col_430, lwd=CO2.lwd, lty=CO2.lty)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #sigmoid line labels
# text(x=log10(1), y=0.32, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold(CO[2]~"alone")))
# text(x=log10(1), y=0.24, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("host odor")))
# text(x=log10(1), y=0.16, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("plant infusion")))
# #sigmoid line significance stars
# text(x=log10(1.2), y=0.26, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[1])
# text(x=log10(1.2), y=0.18, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[2])

#Circle text labels
text(x=xlim[1]-circle_text_x_spacer*(xrange/subplot_width_4col), 
     y=trt_circle_row2_y, adj=c(1,0.5), labels="test")
text(x=xlim[1]-circle_text_x_spacer*(xrange/subplot_width_4col), 
     y=ctr_circle_row2_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.430$x.pos, y=rep(trt_circle_row2_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.430$trt.disp.col,
       col=NA)
points(x=display.info.430$x.pos, y=rep(ctr_circle_row2_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.430$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_430_left, int_box_430_right, int_box_430_right, int_box_430_left),
        y=c(int_box_row2_upper, int_box_row2_upper, int_box_row2_lower, int_box_row2_lower),
        col=NA, border=box_col_430, lwd=box_lwd)

#box labels
text(x=int_box_430_right, y=box_label_row2_y, adj=c(1,1),
     labels=c("435 nm intensity ramps"), col=box_col_430)



#---- Subplot d

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row2,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(subplot_xpos_col[2], subplot_xpos_col[3],
        subplot_ypos[3],subplot_ypos[2]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("D")), side=2, line=subplot_label_xpos2,
      at=subplot_label_ypos2, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0)

#no preference lines
lines(x=c(0,xlim[2]), y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=div.line.1col, y=ylim, col="grey90", lwd=1.5)

#confidence interval around 525 line
polygon(x=log10( c(line.CO2.470.CI$Upper, rev(line.CO2.470.CI$Lower) ) ),
        y=c(line.CO2.470.CI$pref_index, rev(line.CO2.470.CI$pref_index)),
        col=line_col_470_CI, border=NA)

#draws sigmoid model lines
lines(x = log10(line.tansy.470.sub$rel_intensity), y = line.tansy.470.sub$pred_line,
      col=tansy_col, lwd=sig.lwd)
lines(x = log10(line.CO2.470.sub$rel_intensity), y = line.CO2.470.sub$pred_line,
      col=line_col_470, lwd=CO2.lwd, lty=CO2.lty)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #sigmoid line labels
# text(x=log10(1), y=0.32, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold(CO[2]~"alone")))
# text(x=log10(1), y=0.24, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("floral odor")))
# #sigmoid line significance stars
# text(x=log10(1.2), y=0.26, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[3])

#Treatment and control circles
points(x=display.info.470$x.pos, y=rep(trt_circle_row2_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.470$trt.disp.col,
       col=NA)
points(x=display.info.470$x.pos, y=rep(ctr_circle_row2_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.470$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_470_left, int_box_470_right, int_box_470_right, int_box_470_left),
        y=c(int_box_row2_upper, int_box_row2_upper, int_box_row2_lower, int_box_row2_lower),
        col=NA, border=box_col_470, lwd=box_lwd)

#box labels
text(x=int_box_470_right, y=box_label_row2_y, adj=c(1,1),
     labels=c("478 nm intensity ramps"), col=box_col_470)



#---- Subplot e

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row2,m_left,m_top,m_right_outer),
  #specifies subplot position
  fig=c(subplot_xpos_col[3], subplot_xpos_col[5],
        subplot_ypos[3],subplot_ypos[2]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("E")), side=2, line=subplot_label_xpos2,
      at=subplot_label_ypos2, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0)

#no preference lines
lines(x=c(0,xlim[2]), y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=div.line.2col, y=ylim, col="grey90", lwd=1.5)

#confidence interval around 525 line
polygon(x=log10( c(line.CO2.525.CI.sub$Upper, rev(line.CO2.525.CI.sub$Lower) ) ),
        y=c(line.CO2.525.CI.sub$pref_index, rev(line.CO2.525.CI.sub$pref_index)),
        col=line_col_525_CI, border=NA)

#draws sigmoid model lines
lines(x = log10(line.tansy.525.sub$rel_intensity), y = line.tansy.525.sub$pred_line,
      col=tansy_col, lwd=sig.lwd)
lines(x = log10(line.foot.525.sub$rel_intensity), y = line.foot.525.sub$pred_line,
      col=foot_col, lwd=sig.lwd)
lines(x = log10(line.alfalfa.525.sub$rel_intensity), y = line.alfalfa.525.sub$pred_line,
      col=alfalfa_col, lwd=sig.lwd)
lines(x = log10(line.CO2.525.sub$rel_intensity), y = line.CO2.525.sub$pred_line,
      col=line_col_525, lwd=CO2.lwd, lty=CO2.lty)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #sigmoid line labels
# text(x=log10(1), y=0.32, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold(CO[2]~"alone")))
# text(x=log10(1), y=0.24, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("floral odor")))
# text(x=log10(1), y=0.16, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("host odor")))
# text(x=log10(1), y=0.08, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("plant infusion")))
# #sigmoid line significance stars
# text(x=log10(1.1), y=0.26, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[4])
# text(x=log10(1.1), y=0.18, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[5])
# text(x=log10(1.1), y=0.10, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[6])

#Treatment and control circles
points(x=display.info.525$x.pos, y=rep(trt_circle_row2_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.525$trt.disp.col,
       col=NA)
points(x=display.info.525$x.pos, y=rep(ctr_circle_row2_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.525$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_525_left, int_box_525_right, int_box_525_right, int_box_525_left),
        y=c(int_box_row2_upper, int_box_row2_upper, int_box_row2_lower, int_box_row2_lower),
        col=NA, border=box_col_525, lwd=box_lwd)

#box labels
text(x=int_box_525_right, y=box_label_row2_y, adj=c(1,1),
     labels=c("527 nm intensity ramps"), col=box_col_525)



#---- Subplot f

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row3,m_left_outer,m_top,m_right),
  #specifies subplot position
  fig=c(subplot_xpos_col[1], subplot_xpos_col[3],
        subplot_ypos[4],subplot_ypos[3]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=c(xlim[1]+0.1*sp_btw_bars,xlim[2]), ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("F")), side=2, line=subplot_label_xpos2,
      at=subplot_label_ypos2, las=1, cex=subplot_label_cex)

#axis
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp.y, tck=tck, cex.axis=cex.axis)
box(bty="l")
#no preference lines
lines(x=c(xlim[1]+0.1*sp_btw_bars,xlim[2]), y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#axis labels
mtext("preference index", side=2, line=2.2, at=1.065)
mtext("relative intensity", side=1, line=1.2, at=log10(1.85))

#confidence interval around 545 line
polygon(x=log10( c(line.CO2.545.CI$Upper, rev(line.CO2.545.CI$Lower) ) ),
        y=c(line.CO2.545.CI$pref_index, rev(line.CO2.545.CI$pref_index)),
        col=line_col_545_CI, border=NA)

#draws sigmoid model lines
lines(x = log10(line.foot.545.sub$rel_intensity), y = line.foot.545.sub$pred_line,
      col=foot_col, lwd=sig.lwd)
lines(x = log10(line.CO2.545.sub$rel_intensity), y = line.CO2.545.sub$pred_line,
      col=line_col_545, lwd=CO2.lwd, lty=CO2.lty)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #sigmoid line labels
# text(x=log10(1), y=0.32, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold(CO[2]~"alone")))
# text(x=log10(1), y=0.24, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("host odor")))
# #sigmoid line significance stars
# text(x=log10(1.1), y=0.26, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[7])

#Circle text labels
text(x=xlim[1]+0.1*sp_btw_bars-circle_text_x_spacer*(xrange/subplot_width_2col), 
     y=trt_circle_row3_y, adj=c(1,0.5), labels="test")
text(x=xlim[1]+0.1*sp_btw_bars-circle_text_x_spacer*(xrange/subplot_width_2col), 
     y=ctr_circle_row3_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.545$x.pos, y=rep(trt_circle_row3_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.545$trt.disp.col,
       col=NA)
points(x=display.info.545$x.pos, y=rep(ctr_circle_row3_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.545$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_545_left, int_box_545_right, int_box_545_right, int_box_545_left),
        y=c(int_box_row3_upper, int_box_row3_upper, int_box_row3_lower, int_box_row3_lower),
        col=NA, border=box_col_545, lwd=box_lwd)

#box labels
text(x=int_box_545_right, y=box_label_row3_y, adj=c(1,1),
     labels=c("552 nm intensity ramps"), col=box_col_545)



#---- Subplot g

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_row3,m_left,m_top,m_right_outer),
  #specifies subplot position
  fig=c(subplot_xpos_col[3], subplot_xpos_col[5],
        subplot_ypos[4],subplot_ypos[3]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim.reds, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("G")), side=2, line=subplot_label_xpos2,
      at=subplot_label_ypos2, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks.reds, labels=xaxis_tick_labels.reds,
     mgp=axis_mgp.x, tck=tck, cex.axis=cex.axis, gap.axis = 0)

#no preference lines
lines(x=c(0,xlim.reds[2]), y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=div.line.red, y=ylim, col="grey90", lwd=1.5)

#confidence interval around 625 line
#The lower interval becomes very large and is difficult to 
#estimate and is omitted for this panel
# polygon(x=log10( c(line.CO2.625.CI$Upper, 10^xlim.reds[2], 1) ),
#         y=c(line.CO2.625.CI$pref_index, ylim[1], ylim[1]),
#         col=line_col_625_CI, border=NA)

#draws sigmoid model lines
lines(x = log10(line.tansy.625$rel_intensity), y = line.tansy.625$pred_line,
      col=tansy_col, lwd=sig.lwd)
lines(x = log10(line.alfalfa.625$rel_intensity), y = line.alfalfa.625$pred_line,
      col=alfalfa_col, lwd=sig.lwd)
lines(x = log10(line.CO2.625$rel_intensity), y = line.CO2.625$pred_line,
      col=line_col_625, lwd=CO2.lwd, lty=CO2.lty)

#--- These annotations were used to generate the "figure 5 annotations.svg"
#--- file but are omitted for subsequent revisions
# #sigmoid line labels
# text(x=log10(1), y=0.32, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold(CO[2]~"alone")))
# text(x=log10(1), y=0.24, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("floral odor")))
# text(x=log10(1), y=0.16, adj=c(0,0.5), cex=od_label.cex2, 
#      labels=bquote(bold("plant infusion")))
# #sigmoid line significance stars
# text(x=log10(1.1), y=0.26, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[8])
# text(x=log10(1.1), y=0.18, adj=c(0,0.5), labels=sigmoid_comps.odors$sig_stars[9])

#Treatment and control circles
points(x=display.info.625$x.pos, y=rep(trt_circle_row3_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.625$trt.disp.col,
       col=NA)
points(x=display.info.625$x.pos, y=rep(ctr_circle_row3_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.625$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_625_left, int_box_625_right, int_box_625_right, int_box_625_left),
        y=c(int_box_row3_upper, int_box_row3_upper, int_box_row3_lower, int_box_row3_lower),
        col=NA, border=box_col_625, lwd=box_lwd)

#box labels
text(x=int_box_625_right, y=box_label_row3_y, adj=c(1,1),
     labels=c("621 nm intensity ramps"), col=box_col_625)



#Close Graph "Device"
dev.off()
