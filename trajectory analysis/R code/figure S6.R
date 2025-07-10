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
  #omits gray intensity ramp, that data covered in (Fig. 5)
  filter(stim_series != "gray intensity ramp") %>%
  #black vs grey stimuli are labeled according to their stim series
  #so the appropriate black vs gray sorts with the right color
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

#Create two column vector of treatment and control time
#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(round(tracks$trt_time*1000,0), 
            round(tracks$ctr_time*1000,0))

#The pref.model is used to generate individual point and group predictions
#and to test if each group has a preference differing from 50:50

#all stimuli model used for predictions many combinations of
#trt_color:trt_intensity:odor were not performed so model is rank deficient 
#and will give multiple warnings. This model and the ones below will 
#drop trt_color:trt_intensity:odor combinations with no data
pref.model <- glmmTMB(resp ~ 0 + trt_color:trt_intensity:odor + 
                        (1|trt_color:trt_intensity:rep), 
                      data = tracks, family = betabinomial(link = "logit"))

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
           factor(trt_color, levels=levels(tracks$trt_color)),
         trt_intensity = 
           factor(trt_intensity, levels=levels(tracks$trt_intensity)),
         odor = 
           factor(odor, levels=levels(tracks$odor))) %>%
  rename(p_value = "Pr(>|z|)") %>%
  #sorts the table by odor color and then intensity
  arrange(odor, trt_color, trt_intensity)

#Plot residuals to see if model fits
pref.model_simres<-simulateResiduals(pref.model)
plot(pref.model_simres)
testOutliers(pref.model_simres, type="bootstrap")

#Null Model
pref.null <- glmmTMB(resp ~ 0 + (1|trt_color:trt_intensity:rep),  
                     data = tracks, family = betabinomial(link = "logit"))
summary(pref.null)

#Run Anova for overall model fit
anova(pref.model, pref.null)




#---------------------------------------------------------------
#---------------- preference model predictions -----------------
#---------------------------------------------------------------

#creates list of unique combinations of "trt_color",  
#"trt_intensity", "odor", and "rep"
newdata <- 
  tracks %>% select(trt_color, trt_intensity, odor, rep) %>% distinct()

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
  
#creates subsets for each odor
point.pred.CO2 <-
  point.pred %>%
  filter(odor=="CO2") %>%
  mutate(trt_color = droplevels(trt_color))
point.pred.tansy <-
  point.pred %>%
  filter(odor=="CO2 + tansy") %>%
  mutate(trt_color = droplevels(trt_color))
point.pred.foot <-
  point.pred %>%
  filter(odor=="CO2 + foot") %>%
  mutate(trt_color = droplevels(trt_color))
point.pred.alfalfa <-
  point.pred %>%
  filter(odor=="CO2 + alfalfa") %>%
  mutate(trt_color = droplevels(trt_color))


#Extracts overall model predictions, SE and confidence intervals
model.pred<-
  predict_response(pref.model, terms=c("odor","trt_color", "trt_intensity"), 
                   type="fixed", ci_level = 0.95) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(odor=x, trt_color=group, trt_intensity=facet) %>%
  #drops SE column
  select(!std.error) %>%
  #relocate rep column
  relocate(c(trt_color, trt_intensity, odor), .after=NULL) %>%
  #converts from proportion to preference index
  mutate(across(!c(trt_color, trt_intensity, odor), ~ 2 * .x - 1)) %>%
  #joins the p-values from from the "pref.model"
  left_join(pref.model.summary %>% 
              select(odor, trt_color, trt_intensity, p_value), 
            join_by(odor, trt_color, trt_intensity)) %>%
  #this filters out unmeasured factor combinations
  filter(!is.na(p_value))




#--------------------------------------------------------------
#------------ sigmoid intensity model fitting -----------------
#--------------------------------------------------------------

#Read in the coefficients for the sigmoid curve fitted to the gray intensity ramp
sig.pref.gray.coef<-read.csv(file="R Code/gray sigmoid coefficients.csv", row.names = 1)

#The upper (d) and lower (c) asymptotes are taken from the gray model
#These variables are used as fixed terms the sigmoid fitting below.
#The coefficients must be rounded before including them in the medrm call
#not sure why.
coef_c <- sig.pref.gray.coef["c",] %>% round(7)
coef_d <- sig.pref.gray.coef["d",] %>% round(7)


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


#---- sets variables for line predictions below

#sets space between bars in the graph
sp_btw_bars<-0.2125

#color max values for prediction below
max_430<-subset(point.pred, trt_color=="430", select=rel_intensity) %>% max()
max_470<-subset(point.pred, trt_color=="470", select=rel_intensity) %>% max()
max_525<-subset(point.pred, trt_color=="525", select=rel_intensity) %>% max()
max_545<-subset(point.pred, trt_color=="545", select=rel_intensity) %>% max()
max_625<-subset(point.pred, trt_color=="625", select=rel_intensity) %>% max()


#---- CO2 line predictions

line.CO2.430 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_430,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 430", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 0.0*sp_btw_bars + 0 + log10(rel_intensity) / log10(max_430), 
         .before=rel_intensity)

line.CO2.470 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_470,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 470", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 1.0*sp_btw_bars + 1 + log10(rel_intensity) / log10(max_470), 
         .before=rel_intensity)

line.CO2.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_525,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 2.0*sp_btw_bars + 2 + log10(rel_intensity) / log10(max_525), 
         .before=rel_intensity)

line.CO2.545 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_545,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 545", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 3.0*sp_btw_bars + 3 + log10(rel_intensity) / log10(max_545), 
         .before=rel_intensity)

line.CO2.625 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_625,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 | 625", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 4.0*sp_btw_bars + 4 + log10(rel_intensity) / log10(max_625), 
         .before=rel_intensity)


#---- tansy line predictions

line.tansy.470 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_470,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + tansy | 470", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 1.0*sp_btw_bars + 1 + log10(rel_intensity) / log10(max_470), 
         .before=rel_intensity)

line.tansy.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_525,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + tansy | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 2.0*sp_btw_bars + 2 + log10(rel_intensity) / log10(max_525), 
         .before=rel_intensity)

line.tansy.625 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_625,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + tansy | 625", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 4.0*sp_btw_bars + 4 + log10(rel_intensity) / log10(max_625), 
         .before=rel_intensity)


#---- foot line predictions

line.foot.430 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_430,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + foot | 430", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 0.0*sp_btw_bars + 0 + log10(rel_intensity) / log10(max_430), 
         .before=rel_intensity)

line.foot.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_525,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + foot | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 2.0*sp_btw_bars + 2 + log10(rel_intensity) / log10(max_525), 
         .before=rel_intensity)

line.foot.545 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_545,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + foot | 545", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 3.0*sp_btw_bars + 3 + log10(rel_intensity) / log10(max_545), 
         .before=rel_intensity)


#---- alfalfa prediction lines

line.alfalfa.430 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_430,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + alfalfa | 430", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 0.0*sp_btw_bars + 0 + log10(rel_intensity) / log10(max_430), 
         .before=rel_intensity)

line.alfalfa.525 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_525,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + alfalfa | 525", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.no_reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 2.0*sp_btw_bars + 2 + log10(rel_intensity) / log10(max_525), 
         .before=rel_intensity)

line.alfalfa.625 <- 
  #generates a sequence of intensity
  data.frame(rel_intensity = seq(1,max_625,0.01) ) %>%
  #adds columns for "odor_color"
  mutate(odor_color = as.factor( rep("CO2 + alfalfa | 625", length(rel_intensity)) ) ) %>%
  #generates a prediction for each intensity from model
  mutate(pred_line = predict(sig.pref.model.reds, 
                             type="marginal", newdata=.) ) %>%
  #converts from proportion to preference index
  mutate(pred_line = 2*pred_line - 1) %>%
  #converts the relative intensity into an x position
  mutate(xpos = 4.0*sp_btw_bars + 4 + log10(rel_intensity) / log10(max_625), 
         .before=rel_intensity)




#------------------------------------------------------
#-------- Collate display info for the graph ----------
#------------------------------------------------------

#sets space between bars in the graph
sp_btw_bars<-0.2125

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info<-
  #uses tracks as input
  tracks %>%
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
  mutate(x.pos=case_when(
    trt_color == "430" ~ 
      0.0*sp_btw_bars + 
      x.pos / max(subset(., trt_color=="430", select=rel_intensity)),
    trt_color == "470" ~ 
      1.0*sp_btw_bars + 1 +
      x.pos / max(subset(., trt_color=="470", select=rel_intensity)),
    trt_color == "525" ~  
      2.0*sp_btw_bars + 2 +
      x.pos / max(subset(., trt_color=="525", select=rel_intensity)),
    trt_color == "545" ~  
      3.0*sp_btw_bars + 3 + 
      x.pos / max(subset(., trt_color=="545", select=rel_intensity)),
    trt_color == "625" ~  
      4.0*sp_btw_bars + 4 +
      x.pos / max(subset(., trt_color=="625", select=rel_intensity)),
    TRUE ~ x.pos
  )) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=trt_intensity) %>%
  #arrange by "odor" and then "x.pos"
  arrange(odor, x.pos) %>%
  #----this section joins model predictions
  left_join(model.pred, by=c("trt_color", "trt_intensity", "odor")) %>%
  #Determines significance stars for display
  mutate(sig_stars=case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))


#subsets for individual odors
display.info.CO2<-subset(display.info, odor=="CO2")
display.info.tansy<-subset(display.info, odor=="CO2 + tansy")
display.info.foot<-subset(display.info, odor=="CO2 + foot")
display.info.alfalfa<-subset(display.info, odor=="CO2 + alfalfa")




#--------------------------------------------------------------
#----------------- code to generate graph ---------------------
#--------------------------------------------------------------

#------ this section has variables used in the graph call below

#Set figure width, margins and calculate plot sizes and figure height.
#All inputs in inches
f_width <- 11.5
m_bottom <- 0.01
m_bottom_outer <- 1.0
m_left <- 0.5
m_top <- 0.49
m_top_outer <- 0.3
m_right <- 0.1
subplot_height <- 1.81
subplot_width <- f_width-m_left-m_right
subplot_a_height <- subplot_height + m_top_outer + m_bottom
subplot_b_height <- subplot_height + m_top + m_bottom
subplot_c_height <- subplot_height + m_top + m_bottom
subplot_d_height <- subplot_height + m_top + m_bottom_outer
f_height <- sum(subplot_a_height, subplot_b_height, 
                subplot_c_height, subplot_d_height)
#calculates subplot y position
subplot_ypos <- c(1, 
                  1-(subplot_a_height) / f_height, 
                  1-(subplot_a_height + subplot_b_height) / f_height, 
                  1-(subplot_a_height + subplot_b_height + 
                       subplot_c_height) / f_height,
                  0)

#x limits and tick labels
xlim<-c(-0.5*sp_btw_bars, max(display.info$x.pos) + 0.5*sp_btw_bars)
xrange<-xlim[2]-xlim[1]
#calculates tick labels and positions from "rel_intensity"
rel_max_430<-10^max(subset(display.info, trt_color=="430", select=rel_intensity))
rel_max_470<-10^max(subset(display.info, trt_color=="470", select=rel_intensity))
rel_max_525<-10^max(subset(display.info, trt_color=="525", select=rel_intensity))
rel_max_545<-10^max(subset(display.info, trt_color=="545", select=rel_intensity))
rel_max_625<-10^max(subset(display.info, trt_color=="625", select=rel_intensity))
xaxis_ticks<-
  c( 0.0*sp_btw_bars + 0 + log10(seq(1, 1.05*rel_max_430, 0.2)) / log10(rel_max_430),
     1.0*sp_btw_bars + 1 + log10(seq(1, 1.05*rel_max_470, 0.2)) / log10(rel_max_470),
     2.0*sp_btw_bars + 2 + log10(seq(1, 1.05*rel_max_525, 0.2)) / log10(rel_max_525),
     3.0*sp_btw_bars + 3 + log10(seq(1, 1.05*rel_max_545, 0.2)) / log10(rel_max_545),
     4.0*sp_btw_bars + 4 + log10(seq(1, 1.05*rel_max_625, 0.2)) / log10(rel_max_625))
xaxis_tick_labels<-
  c( seq(1, 1.05*rel_max_430, 0.2),
     seq(1, 1.05*rel_max_470, 0.2),
     seq(1, 1.05*rel_max_525, 0.2),
     seq(1, 1.05*rel_max_545, 0.2),
     seq(1, 1.05*rel_max_625, 0.2)) %>%
  format(digits = 2)

#Y limits and tick labels
ylim<-c(-0.25, 0.70)
yrange<-ylim[2]-ylim[1]
yaxis_ticks<-seq(ylim[1], ylim[2], 0.25)

axis_mgp <- c(3,0.35,0) #2nd value controls position of tick labels
tck <- -0.015 #Value to control tick mark length default is -0.01
#number of treatments
n_trts <- length(levels(point.pred.CO2$trt_color)) * 
  length(levels(point.pred.CO2$trt_intensity)) 
CI_width <- 0.14*sp_btw_bars #half width of confidence interval rectangles in plot units
jitter <- 0.078*sp_btw_bars #half width of the jitter cloud
mean_lwd <- 2 #width of the mean line
sig.lwd <- 1.5 #width of the sigmoid lines
sig.col <- "grey60" # color of the sigmoid model predicted lines
#xaxis line positions
xaxis_pos <- c(xlim[1],                0.25*sp_btw_bars + 1,
               0.75*sp_btw_bars + 1,   1.25*sp_btw_bars + 2,
               1.75*sp_btw_bars + 2,   2.25*sp_btw_bars + 3,
               2.75*sp_btw_bars + 3,   3.25*sp_btw_bars + 4,
               3.75*sp_btw_bars + 4,   xlim[2])
#alternating ypos of sample size text
sample_size_ypos <- c(ylim[2]+0.04*yrange, ylim[2]+0.10*yrange)
sample_size_cex <- 5/8 #size of sample size text
subplot_label_cex <- 12/8 #size of subplot levels
subplot_label_xpos <- 0.6 #x position of subplot labels in lines (0.2"/line)
subplot_label_ypos <- ylim[2] #y position of subplot labels in plot units
cex.axis <- 6/8 #tick label font size
sample_size_col <- "gray70"
#space between the bars and the stars
sig_star_spacer <- 0.10*yrange
sig_star_cex <- 1 #9.5/8 #size of significance stars

#line styles
np_line_col <- "gray70"
np_line.lty <- 2 #no preference line type
np_line.lwd <- 1.5 #no preference line width
np_line.lend <- 1 #line end type "butt"

#odor diagram variables
od_xcenter <- 0.19
od_xcenter_npc <- 
  (od_xcenter - (xlim[1]-m_left*xrange/subplot_width) ) / 
  (f_width*xrange/subplot_width)
od.x1 <- unit(od_xcenter_npc, "npc")
od.x2 <- unit(od_xcenter_npc-0.013, "npc")
od.x3 <- unit(od_xcenter_npc+0.021, "npc")
od.x4 <- unit(od_xcenter_npc+0.020, "npc")
od.x5 <- unit(od_xcenter_npc+0.016, "npc")
od.y.spacer <- -0.215
od.y1 <- unit(subplot_ypos[1] + od.y.spacer + 0.019, "npc")
od.y2 <- unit(subplot_ypos[2] + od.y.spacer + 0.000, "npc")
od.y3 <- unit(subplot_ypos[2] + od.y.spacer + 0.005, "npc")
od.y4 <- unit(subplot_ypos[3] + od.y.spacer + 0.000, "npc")
od.y5 <- unit(subplot_ypos[3] + od.y.spacer + 0.000, "npc")
od.y6 <- unit(subplot_ypos[4] + od.y.spacer + 0.003, "npc")
od.y7 <- unit(subplot_ypos[4] + od.y.spacer + 0.000, "npc")
od.width <- unit(0.75, "inch")
od_label.x <- od_xcenter
od_label.y <- -0.18
od_label.cex <-1

#Variables for circle size and positioning
#incorporating yrange and ylim[1] keeps the circle position constant 
#when ylims change
trt_circle_y <- ylim[1]-0.300*yrange 
ctr_circle_y <- trt_circle_y-0.090*yrange #expressed relative to "trt_circle_y"
circle_text_x_spacer <- 0.12 #space between the circles and their labels plot units
circle_size<-3 #character expansion factor default=1
trt_circle_outline <- NA
ctr_circle_outline <- NA

#Variables for the positioning of the boxes around the circles
box_lwd <- 1.25
int_box_upper <- trt_circle_y + 0.05*yrange
int_box_lower <- ctr_circle_y - 0.05*yrange
int_box_1_left <- display.info.CO2$x.pos[1] - 0.4*sp_btw_bars
int_box_1_right <- display.info.CO2$x.pos[8] + 0.4*sp_btw_bars
int_box_2_left <- display.info.CO2$x.pos[9] - 0.4*sp_btw_bars
int_box_2_right <- display.info.CO2$x.pos[16] + 0.4*sp_btw_bars
int_box_3_left <- display.info.CO2$x.pos[17] - 0.4*sp_btw_bars
int_box_3_right <- display.info.CO2$x.pos[24] + 0.4*sp_btw_bars
int_box_4_left <- display.info.CO2$x.pos[25] - 0.4*sp_btw_bars
int_box_4_right <- display.info.CO2$x.pos[32] + 0.4*sp_btw_bars
int_box_5_left <- display.info.CO2$x.pos[33] - 0.4*sp_btw_bars
int_box_5_right <- display.info.CO2$x.pos[40] + 0.4*sp_btw_bars
box_col_430 <-  subset(display.info, trt_color=="430" & trt_intensity=="1.00", 
                       select=trt.disp.col) %>% unique() %>%  as.character()
box_col_470 <-  subset(display.info, trt_color=="470" & trt_intensity=="1.00", 
                       select=trt.disp.col) %>% unique() %>%  as.character()
box_col_525 <-  subset(display.info, trt_color=="525" & trt_intensity=="1.00", 
                       select=trt.disp.col) %>% unique() %>%  as.character()
box_col_545 <-  subset(display.info, trt_color=="545" & trt_intensity=="1.00", 
                       select=trt.disp.col) %>% unique() %>%  as.character()
box_col_625 <-  subset(display.info, trt_color=="625" & trt_intensity=="1.00", 
                       select=trt.disp.col) %>% unique() %>%  as.character()
odor_box_col <- "gray70"
box_label_y <- int_box_lower-0.02*yrange



#-------- Execute code from here to end to generate graph

# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure S6.pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)

#Export figure to a png with the width and height from above
#also sets font type and size
# png(file=paste0(output_dir,"figure S6.png"),
#     width=f_width, height=f_height, units="in", res=150,
#     family="sans", pointsize=8)


par(
  #specifies the outer (whole plot) margin in inches
  omi=c(0,0,0,0),
  #allows drawing outside of the plotting area
  xpd = TRUE)



#---- Subplot a

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom,m_left,m_top_outer,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[2],subplot_ypos[1]),
  cex=1) #prevents font scaling in multi-panel figures

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("A")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos, las=1, cex=subplot_label_cex)

#axis
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)

#no preference lines
lines(x=xaxis_pos[1:2], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[3:4], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[5:6], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[7:8], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[9:10], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim, col="grey90", lwd=1.5)

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.CO2$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width)
          + display.info.CO2$x.pos[i],
          y=c(rep(display.info.CO2$conf.low[i],2),
              rep(display.info.CO2$conf.high[i],2)),
          col=display.info.CO2$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.CO2$x.pos[i],
        y=rep(display.info.CO2$predicted[i],2),
        col=display.info.CO2$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#set seed for jitter
set.seed(123456)

#stripchart
stripchart(point_pred~trt_intensity*trt_color, data=point.pred.CO2, 
           at=display.info.CO2$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.CO2$trt.disp.col.alpha)

#draws sigmoid model lines
lines(x = line.CO2.430$xpos, y = line.CO2.430$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.CO2.470$xpos, y = line.CO2.470$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.CO2.525$xpos, y = line.CO2.525$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.CO2.545$xpos, y = line.CO2.545$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.CO2.625$xpos, y = line.CO2.625$pred_line, 
      col=sig.col, lwd=sig.lwd)

#significance stars
text(x=display.info.CO2$x.pos,
     y=display.info.CO2$conf.high+sig_star_spacer,
     labels=display.info.CO2$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  tracks %>% filter(odor == "CO2") %>% 
  group_by(trt_color, trt_intensity, odor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = length(trt_color), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.CO2$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#white line covering the no preference line under the odor diagrams
lines(x = c(0.00,0.38), y = c(0,0), lwd=np_line.lwd, col="white")
#odor diagrams
grid.picture(CO2_SVG, x=od.x1, y=od.y1, just=c(0.5,0), width=od.width)
text(x=od_label.x, y=od_label.y, labels=bquote(bold(CO[2]~"alone")),
     adj=c(0.5,0.5), cex=od_label.cex)



#---- Subplot b

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[3],subplot_ypos[2]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("B")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)

#no preference lines
lines(x=xaxis_pos[3:4], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[5:6], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[9:10], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim, col="grey90", lwd=1.5)

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.tansy$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width)
          + display.info.tansy$x.pos[i],
          y=c(rep(display.info.tansy$conf.low[i],2),
              rep(display.info.tansy$conf.high[i],2)),
          col=display.info.tansy$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.tansy$x.pos[i],
        y=rep(display.info.tansy$predicted[i],2),
        col=display.info.tansy$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#set seed for jitter
set.seed(123456)

#stripchart
stripchart(point_pred~trt_intensity*trt_color, data=point.pred.tansy, 
           at=display.info.tansy$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.tansy$trt.disp.col.alpha)

#draws sigmoid model lines
lines(x = line.tansy.470$xpos, y = line.tansy.470$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.tansy.525$xpos, y = line.tansy.525$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.tansy.625$xpos, y = line.tansy.625$pred_line, 
      col=sig.col, lwd=sig.lwd)

#significance stars
text(x=display.info.tansy$x.pos,
     y=display.info.tansy$conf.high+sig_star_spacer,
     labels=display.info.tansy$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  tracks %>% filter(odor == "CO2 + tansy") %>% 
  group_by(trt_color, trt_intensity, odor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = length(trt_color), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.tansy$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y2, just=c(0.5,0), width=od.width)
grid.picture(tansy_SVG, x=od.x3, y=od.y3, just=c(0.5,0), height=od.width*0.8)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("floral odor")),
     adj=c(0.5,0.5), cex=od_label.cex)



#---- Subplot c

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[4],subplot_ypos[3]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("C")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)

#no preference lines
lines(x=xaxis_pos[1:2], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[5:6], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[7:8], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim, col="grey90", lwd=1.5)

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.foot$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width)
          + display.info.foot$x.pos[i],
          y=c(rep(display.info.foot$conf.low[i],2),
              rep(display.info.foot$conf.high[i],2)),
          col=display.info.foot$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.foot$x.pos[i],
        y=rep(display.info.foot$predicted[i],2),
        col=display.info.foot$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#set seed for jitter
set.seed(123456)

#stripchart
stripchart(point_pred~trt_intensity*trt_color, data=point.pred.foot, 
           at=display.info.foot$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.foot$trt.disp.col.alpha)

#draws sigmoid model lines
lines(x = line.foot.430$xpos, y = line.foot.430$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.foot.525$xpos, y = line.foot.525$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.foot.545$xpos, y = line.foot.545$pred_line, 
      col=sig.col, lwd=sig.lwd)

#significance stars
text(x=display.info.foot$x.pos,
     y=display.info.foot$conf.high+sig_star_spacer,
     labels=display.info.foot$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  tracks %>% filter(odor == "CO2 + foot") %>% 
  group_by(trt_color, trt_intensity, odor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = length(trt_color), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.foot$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#white line covering the no preference line under the odor diagrams
lines(x = c(-0.05,0.38), y = c(0,0), lwd=np_line.lwd, col="white")
#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y4, just=c(0.5,0), width=od.width)
grid.picture(foot_SVG, x=od.x4, y=od.y5, just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("host odor")),
     adj=c(0.5,0.5), cex=od_label.cex)



#---- Subplot d

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom_outer,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[5],subplot_ypos[4]),
  #prevents font scaling in multi-panel figures
  cex=1, 
  new=T)

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("D")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks, labels=xaxis_tick_labels,
     mgp=axis_mgp, tck=tck, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
lines(x=xaxis_pos[1:2], y=rep(ylim[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim[1],2))
lines(x=xaxis_pos[7:8], y=rep(ylim[1],2))
lines(x=xaxis_pos[9:10], y=rep(ylim[1],2))
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)

#no preference lines
lines(x=xaxis_pos[1:2], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[5:6], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)
lines(x=xaxis_pos[9:10], y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim, col="grey90", lwd=1.5)

#axis labels
mtext("relative intensity", side=1, line=1.6)
mtext("preference index", side=2, line=2.2, at=2.0)

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.alfalfa$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width)
          + display.info.alfalfa$x.pos[i],
          y=c(rep(display.info.alfalfa$conf.low[i],2),
              rep(display.info.alfalfa$conf.high[i],2)),
          col=display.info.alfalfa$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.alfalfa$x.pos[i],
        y=rep(display.info.alfalfa$predicted[i],2),
        col=display.info.alfalfa$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#set seed for jitter
set.seed(123456)

# #stripchart
stripchart(point_pred~trt_intensity*trt_color, data=point.pred.alfalfa,
           at=display.info.alfalfa$x.pos, method = "jitter", jitter=jitter,
           vertical=TRUE, add=TRUE, pch=16,
           col=display.info.alfalfa$trt.disp.col.alpha)

#draws sigmoid model lines
lines(x = line.alfalfa.430$xpos, y = line.alfalfa.430$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.alfalfa.525$xpos, y = line.alfalfa.525$pred_line, 
      col=sig.col, lwd=sig.lwd)
lines(x = line.alfalfa.625$xpos, y = line.alfalfa.625$pred_line, 
      col=sig.col, lwd=sig.lwd)

#significance stars
text(x=display.info.alfalfa$x.pos,
     y=display.info.alfalfa$conf.high+sig_star_spacer,
     labels=display.info.alfalfa$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  tracks %>% filter(odor == "CO2 + alfalfa") %>% 
  group_by(trt_color, trt_intensity, odor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = length(trt_color), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.alfalfa$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#white line covering the no preference line under the odor diagrams
lines(x = c(-0.05,0.38), y = c(0,0), lwd=np_line.lwd, col="white")
#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y6, just=c(0.5,0), width=od.width)
grid.picture(alfalfa_SVG, x=od.x5, y=od.y7,just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("plant infusion")),
     adj=c(0.5,0.5), cex=od_label.cex)

#Circle text labels
text(x=min(display.info.CO2$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.CO2$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.CO2$x.pos, y=rep(trt_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.CO2$trt.disp.col,
       col=NA)
points(x=display.info.CO2$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.CO2$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_1_left, int_box_1_right, int_box_1_right, int_box_1_left),
        y=c(int_box_upper, int_box_upper, int_box_lower, int_box_lower),
        col=NA, border=box_col_430, lwd=box_lwd)
polygon(x=c(int_box_2_left, int_box_2_right, int_box_2_right, int_box_2_left),
        y=c(int_box_upper, int_box_upper, int_box_lower, int_box_lower),
        col=NA, border=box_col_470, lwd=box_lwd)
polygon(x=c(int_box_3_left, int_box_3_right, int_box_3_right, int_box_3_left),
        y=c(int_box_upper, int_box_upper, int_box_lower, int_box_lower),
        col=NA, border=box_col_525, lwd=box_lwd)
polygon(x=c(int_box_4_left, int_box_4_right, int_box_4_right, int_box_4_left),
        y=c(int_box_upper, int_box_upper, int_box_lower, int_box_lower),
        col=NA, border=box_col_545, lwd=box_lwd)
polygon(x=c(int_box_5_left, int_box_5_right, int_box_5_right, int_box_5_left),
        y=c(int_box_upper, int_box_upper, int_box_lower, int_box_lower),
        col=NA, border=box_col_625, lwd=box_lwd)
#box labels
text(x=c(int_box_1_right, int_box_2_right, int_box_3_right, 
         int_box_4_right, int_box_5_right), 
     y=box_label_y, adj=c(1,1),
     labels=c("435 nm intensity ramp","478 nm intensity ramp",
              "527 nm intensity ramp","552 nm intensity ramp",
              "621 nm intensity ramp"),
     col=c(box_col_430, box_col_470, box_col_525, 
           box_col_545, box_col_625))

#Close Graph "Device"
dev.off()



