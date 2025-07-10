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
tansy_SVG <- file.path(odor_source_dir, "tansy-cairo.svg") %>% readPicture()
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

#extracts the violet, grey and yellow CO2 runs and renames the odor
#so the control stimulus pairs can be be estimated separately from
#the blue, green, and red runs
events<-event_extract(sub_dir_CO2_VGrY) %>%
  mutate(odor = case_when(
    odor=="CO2" ~ "CO2_VGrY",
    TRUE ~ odor
  ))

#Extracts track level data from the remaining groups of runs
events<-event_extract(sub_dir_CO2, events)
events<-event_extract(sub_dir_CO2_tansy, events)
events<-event_extract(sub_dir_CO2_foot, events)
events<-event_extract(sub_dir_CO2_alfalfa, events)


#list of rep odors, used for relabeling the air controls
rep_vs_odor <- events %>% select(rep, odor) %>% filter(!odor=="air") %>% distinct()

#Creates a subset of data and processes it for analysis and display
events<-
  #starts with "events"
  events %>%
  #joins list of rep odors
  left_join(rep_vs_odor, by="rep") %>%
  #replace original odor with joined odor
  relocate(odor.y, .before=odor.x) %>%
  select(!odor.x) %>%
  rename(odor=odor.y) %>%
  #turns odor into a factor
  mutate(odor = factor(odor, levels=c("CO2","CO2_VGrY","CO2 + tansy",
                                      "CO2 + foot","CO2 + alfalfa"))) %>%
  #sets the treatment factor
  mutate(treatment_factor = case_when(
    #controls specified by individual stim_series values
    stim_series=="PreCO2Time" ~ "PreCO2Time",
    stim_series=="resp check, initial no odor" ~ "resp check, initial no odor",
    stim_series=="PreStim" ~ "PreStim",
    stim_series=="resp check, initial odor" ~ "resp check, initial odor",
    stim_series=="PostCO2Time" ~ "PostCO2Time",
    #black vs grey stimuli are labeled according to their stim series
    #so the appropriate black vs gray sorts with the right color
    #intensity ramp
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="430 nm intensity ramp" ~ "430-0.00",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="470 nm intensity ramp" ~ "470-0.00",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="525 nm intensity ramp" ~ "525-0.00",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="545 nm intensity ramp" ~ "545-0.00",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="625 nm intensity ramp" ~ "625-0.00",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      stim_series=="gray intensity ramp" ~ "gray-0.00",
    #all others are assigned color-intensity
    TRUE ~ paste0(trt_color,"-",trt_intensity)), 
    .after = stim_series) %>%
  #Turns "treatment_factor" and rep into a factor 
  mutate(treatment_factor = as.factor(treatment_factor), 
         rep = as.factor(rep)) %>%
  #Moves the non color factor levels before the others
  mutate(treatment_factor = fct_relevel(treatment_factor,
                                        "PreCO2Time",
                                        "resp check, initial no odor",
                                        "PreStim",
                                        "resp check, initial odor")) %>%
  #Move the color intensity blacks before the relevant color's 0.20 intensity
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, "430-0.00",
                       #the after formula here looks up the position of the treatment level
                       #and puts the above treatment level before the below treatment level
                       after=match("430-0.20", levels(treatment_factor))-2)) %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, "470-0.00",
                       after=match("470-0.20", levels(treatment_factor))-2)) %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, "525-0.00",
                       after=match("525-0.20", levels(treatment_factor))-2)) %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, "545-0.00",
                       after=match("545-0.20", levels(treatment_factor))-2)) %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, "625-0.00",
                       after=match("625-0.20", levels(treatment_factor))-2))  %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, "gray-0.00",
                       after=match("gray-0.20", levels(treatment_factor))-2))  %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, 
                       grep("gray", unique(treatment_factor), value=TRUE),
                       after=match("525-3.00", levels(treatment_factor))-2))  %>%
  #black vs grey stimuli are labeled according to their stim series
  #so the appropriate black vs gray sorts with the right color
  #intensity ramp
  mutate(trt_color = case_when(
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="1.00" &
      treatment_factor=="resp check, initial no odor" ~ "gray",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="1.00" &
      treatment_factor=="resp check, initial odor" ~ "gray",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      treatment_factor=="430-0.00" ~ "430",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      treatment_factor=="470-0.00" ~ "470",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      treatment_factor=="525-0.00" ~ "525",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      treatment_factor=="545-0.00" ~ "545",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      treatment_factor=="625-0.00" ~ "625",
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50" &
      treatment_factor=="gray-0.00" ~ "gray",
    #all trt_colors are left unmodified
    TRUE ~ trt_color)) %>%
  #turns these variables numeric
  mutate_at(c("activity_index", "mean_duration", "mean_track_speed",
              "n_tracks", "n_tracks_resp", "p_tracks_resp", 
              "mean_pref_index",  "conf.low", "conf.high", "p_value"), 
            as.numeric)

#Displays the total number of reps for each treatment_factor & odors
table(list(events$treatment_factor, events$odor))




#------------------------------------------------------
#-------- mixed model analysis of activation ----------
#------------------------------------------------------

#all stimuli model used for predictions many combinations of
#treatment_factor:odor were not performed so model is rank deficient 
#and will give multiple warnings. This model and the ones below will 
#drop treatment_factor:odor combinations with no data
activ.model <- 
  glmmTMB(n_tracks ~ 0 + odor + treatment_factor:odor + (1|rep), 
          data = events, family = nbinom1(link = "log"))
#Processes model summary to remove dropped treatment_factor:odor combinations
activ.model.summary <- 
  #extracts coefficients
  summary(activ.model)$coefficients$cond %>% 
  as.data.frame() %>% 
  #filters the dropped treatment_factor:odors
  filter(!is.na(Estimate)) %>% 
  #turns the row names into a column
  tibble::rownames_to_column() %>% 
  #splits the column to extract factor info
  separate_wider_delim(rowname, delim = ":",  
                       names = c("odor", "treatment_factor"), 
                       too_few="align_start") %>%
  #removes factor labels from the levels
  mutate(odor = gsub("odor", "", odor),
         treatment_factor = gsub("treatment_factor", "", treatment_factor)) %>%
  #adds "PreCO2Time" to the first factor level, missing in summary
  mutate(treatment_factor = case_when(
    is.na(treatment_factor) ~ "PreCO2Time",
    TRUE ~ treatment_factor)) %>%
  #turns the new columns into factors using levels from "events
  mutate(odor = factor(odor, levels=levels(events$odor)),
         treatment_factor = factor(treatment_factor, levels=levels(events$treatment_factor))) %>%
  #sorts the table by odor and then by treatment_factor
  arrange(odor, treatment_factor)


#Plot residuals to see if model fits
activ.model_simres<-simulateResiduals(activ.model)
plot(activ.model_simres)
testOutliers(activ.model_simres, type="bootstrap")

#Null Model
activ.null <- glmmTMB(resp ~ (1|rep), data = events, 
                     family = betabinomial(link = "logit"))
summary(activ.null)

#Run Anova for overall model fit
anova(activ.model, activ.null)


#---- Analysis with only intensity ramps
events.sub<-
  events %>%
  #Filters to only spectral sweeps
  filter(str_detect(stim_series, "intensity ramp")) %>%
  mutate(treatment_factor = droplevels(treatment_factor))

#Full Model
activ.model.1 <- glmmTMB(n_tracks ~ trt_color + trt_intensity + odor + 
                          trt_color:trt_intensity + trt_color:odor + trt_intensity:odor + 
                          trt_color:trt_intensity:odor +(1|rep), 
                        data = events.sub, family = nbinom1(link = "log"),
                        control = glmmTMBControl(optCtrl = list(iter.max=1000, eval.max=1100)))
summary(activ.model.1)

#Plot residuals to see if model fits
activ.model.1_simres<-simulateResiduals(activ.model.1)
plot(activ.model.1_simres)
testOutliers(activ.model.1_simres, type="bootstrap")


#no 3-way interaction model
activ.model.2 <- glmmTMB(n_tracks ~ trt_color + trt_intensity + odor + 
                          trt_color:trt_intensity + trt_color:odor + trt_intensity:odor + 
                          (1|rep), 
                         data = events.sub, family = nbinom1(link = "log"),
                         control = glmmTMBControl(optCtrl = list(iter.max=1000, eval.max=1100)))
summary(activ.model.2)

#anova testing the 3-way interaction
anova(activ.model.1, activ.model.2)

#no odor model
activ.model.3 <- glmmTMB(n_tracks ~ trt_color + trt_intensity + 
                          trt_color:trt_intensity + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.3)

#anova testing the effect of odor
anova(activ.model.1, activ.model.3)


#no color model
activ.model.4 <- glmmTMB(n_tracks ~ trt_intensity + odor + 
                          trt_intensity:odor + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.4)

#anova testing the effect of color
anova(activ.model.2, activ.model.4)


#no intensity model
activ.model.5 <- glmmTMB(n_tracks ~ trt_color + odor + 
                          trt_color:odor + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.5)

#anova testing the effect of intensity
anova(activ.model.2, activ.model.5)


#null model
activ.null <- glmmTMB(n_tracks ~ (1|rep), 
                     data = events.sub, 
                     family = nbinom1(link = "log"))
summary(activ.null)

#anova testing overall model fit
anova(activ.model.1, activ.null)


#minimally adequate model
activ.model.qq <- glmmTMB(n_tracks ~ trt_intensity + trt_color:odor + (1|rep), 
                          data = events.sub, 
                          family = nbinom1(link = "log"))
summary(activ.model.11)

#test of model simplification
anova(activ.model.2, activ.model.11)




#---------------------------------------------------------------
#----------- Extract activation model predictions --------------
#---------------------------------------------------------------

#creates list of unique combinations of "treatment_factor", "odor", "rep"
newdata <- events %>% select(treatment_factor, odor, rep) %>% distinct()

#point predictions
point.pred <-
  predict(activ.model, re.form = NULL, type="response", newdata=newdata) %>%
  #combine with "newdata" frame
  cbind(newdata, .) %>%
  #rename prediction column
  rename(point_pred = ".")

#creates subsets for each odor
point.pred.CO2 <-
  point.pred %>%
  filter(odor=="CO2") %>%
  mutate(treatment_factor = droplevels(treatment_factor))
point.pred.CO2_VGrY <-
  point.pred %>%
  filter(odor=="CO2_VGrY") %>%
  mutate(treatment_factor = droplevels(treatment_factor))
point.pred.tansy <-
  point.pred %>%
  filter(odor=="CO2 + tansy") %>%
  mutate(treatment_factor = droplevels(treatment_factor))
point.pred.foot <-
  point.pred %>%
  filter(odor=="CO2 + foot") %>%
  mutate(treatment_factor = droplevels(treatment_factor))
point.pred.alfalfa <-
  point.pred %>%
  filter(odor=="CO2 + alfalfa") %>%
  mutate(treatment_factor = droplevels(treatment_factor))


#Extracts overall model predictions, SE and confidence intervals
model.pred<-
  predict_response(activ.model, terms=c("odor","treatment_factor"), 
                   type="fixed", ci_level = 0.95) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(treatment_factor=group, odor=x) %>%
  #drops SE column
  select(!std.error) %>%
  #relocate rep column
  relocate(c(odor), .after=treatment_factor) %>%
  left_join(activ.model.summary %>% select(odor, treatment_factor, "Pr(>|z|)"), 
            join_by(odor, treatment_factor)) %>%
  rename(p_value = "Pr(>|z|)") %>%
  filter(!is.na(p_value))




#------------------------------------------------------
#-------- Collate display info for the graph ----------
#------------------------------------------------------

#sets space between bars in the graph
sp_btw_bars<-0.2125

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info<-
  #uses tracks as input
  events %>%
  #----this section summarizes by trt_color and odor
  group_by(treatment_factor, odor) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_color=first(trt_color),
            trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups="drop") %>%
  #sort by "odor" and then "trt_color"
  arrange(odor, treatment_factor) %>%
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
  #Adds relative intensity column from "rel_intensity" read in above
  left_join(rel_intensity, by = join_by(trt_intensity == int_setting, 
                                        trt_color == color)) %>%
  mutate(rel_intensity = log10(as.numeric(rel_intensity))) %>%
  #moves new intensity column after "trt_intensity"
  relocate(rel_intensity, .after=trt_intensity) %>%
  #----This section determines x positions to display data
  mutate(x.pos=rel_intensity) %>%
  mutate(x.pos=case_when(
    treatment_factor == "PreCO2Time" ~ 0,
    treatment_factor == "resp check, initial no odor" ~  0.5*sp_btw_bars,
    treatment_factor == "PreStim" ~  1.5*sp_btw_bars,
    treatment_factor == "resp check, initial odor" ~ 2.0*sp_btw_bars,
    str_detect(treatment_factor, "470")  ~ 
      3.0*sp_btw_bars + 
      x.pos / max(subset(., trt_color=="470", select=rel_intensity)),
    str_detect(treatment_factor, "430") ~ 
      3.0*sp_btw_bars + 
      x.pos / max(subset(., trt_color=="430", select=rel_intensity)),
    str_detect(treatment_factor, "525")  ~ 
      4.0*sp_btw_bars + 1 +
      x.pos / max(subset(., trt_color=="525", select=rel_intensity)),
    str_detect(treatment_factor, "gray")  ~ 
      4.0*sp_btw_bars + 1 + 
      x.pos / max(subset(., trt_color=="gray", select=rel_intensity)),
    str_detect(treatment_factor, "625")  ~ 
      5.0*sp_btw_bars + 2 +
      x.pos / max(subset(., trt_color=="625", select=rel_intensity)),
    str_detect(treatment_factor, "545") ~ 
      5.0*sp_btw_bars + 2 + 
      x.pos / max(subset(., trt_color=="545", select=rel_intensity)),
    treatment_factor == "PostCO2Time" ~ 
      6.0*sp_btw_bars + 3,
    TRUE ~ x.pos
  )) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=treatment_factor) %>%
  #arrange by "odor" and then "x.pos"
  arrange(odor, x.pos) %>%
  #----this section joins model predictions
  left_join(model.pred, by=c("treatment_factor", "odor")) %>%
  #Determines significance stars for display
  mutate(sig_stars=case_when(
    treatment_factor == "PreCO2Time" ~ as.character(NA),
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ as.character(NA) ))


#subsets for individual odors
display.info.CO2<-subset(display.info, odor=="CO2")
display.info.CO2_VGrY<-subset(display.info, odor=="CO2_VGrY")
display.info.tansy<-subset(display.info, odor=="CO2 + tansy")
display.info.foot<-subset(display.info, odor=="CO2 + foot")
display.info.alfalfa<-subset(display.info, odor=="CO2 + alfalfa")




#--------------------------------------------------------------
#----------------- code to generate graph ---------------------
#--------------------------------------------------------------

#------ this section has variables used in the graph call below

#x limits and tick labels
xlim<-c(-1.7*sp_btw_bars, max(display.info$x.pos) + 0.3*sp_btw_bars)
xrange<-xlim[2]-xlim[1]

#calculates tick labels and positions from "rel_intensity"
rel_max_430<-10^max(subset(display.info, trt_color=="430", select=rel_intensity))
rel_max_470<-10^max(subset(display.info, trt_color=="470", select=rel_intensity))
rel_max_525<-10^max(subset(display.info, trt_color=="525", select=rel_intensity))
rel_max_545<-10^max(subset(display.info, trt_color=="545", select=rel_intensity))
rel_max_625<-10^max(subset(display.info, trt_color=="625", select=rel_intensity))
rel_max_gray<-10^max(subset(display.info, trt_color=="gray", select=rel_intensity))
xaxis_ticks_CO2<-
  c( 3.0*sp_btw_bars + log10(seq(1, 1.05*rel_max_470, 0.2)) / log10(rel_max_470),
     4.0*sp_btw_bars + 1 + log10(seq(1, 1.05*rel_max_525, 0.2)) / log10(rel_max_525),
     5.0*sp_btw_bars + 2 + log10(seq(1, 1.05*rel_max_625, 0.2)) / log10(rel_max_625))
xaxis_ticks_CO2_VGrY<-
  c( 3.0*sp_btw_bars + log10(seq(1, 1.05*rel_max_430, 0.2)) / log10(rel_max_470),
     4.0*sp_btw_bars + 1 + log10(seq(1, 1.05*rel_max_gray, 1.0)) / log10(rel_max_gray),
     5.0*sp_btw_bars + 2 + log10(seq(1, 1.05*rel_max_545, 0.2)) / log10(rel_max_545))
xaxis_ticks_tansy<-
  c( 3.0*sp_btw_bars + log10(seq(1, 1.05*rel_max_470, 0.2)) / log10(rel_max_470),
     4.0*sp_btw_bars + 1 + log10(seq(1, 1.05*rel_max_525, 0.2)) / log10(rel_max_525),
     5.0*sp_btw_bars + 2 + log10(seq(1, 1.05*rel_max_625, 0.2)) / log10(rel_max_625))
xaxis_ticks_foot<-
  c( 3.0*sp_btw_bars + log10(seq(1, 1.05*rel_max_430, 0.2)) / log10(rel_max_430),
     4.0*sp_btw_bars + 1 + log10(seq(1, 1.05*rel_max_525, 0.2)) / log10(rel_max_525),
     5.0*sp_btw_bars + 2 + log10(seq(1, 1.05*rel_max_545, 0.2)) / log10(rel_max_545))
xaxis_ticks_alfalfa<-
  c( 3.0*sp_btw_bars + log10(seq(1, 1.05*rel_max_430, 0.2)) / log10(rel_max_430),
     4.0*sp_btw_bars + 1 + log10(seq(1, 1.05*rel_max_525, 0.2)) / log10(rel_max_525),
     5.0*sp_btw_bars + 2 + log10(seq(1, 1.05*rel_max_625, 0.2)) / log10(rel_max_625))
xaxis_tick_labels_CO2<-
  c( seq(1, 1.05*rel_max_470, 0.2),
     seq(1, 1.05*rel_max_525, 0.2),
     seq(1, 1.05*rel_max_625, 0.2)) %>%
  format(digits = 2)
xaxis_tick_labels_CO2_VGrY<-
  c( seq(1, 1.05*rel_max_430, 0.2),
     seq(1, 1.05*rel_max_gray, 1.0),
     seq(1, 1.05*rel_max_545, 0.2)) %>%
  format(digits = 2)
xaxis_tick_labels_tansy<-
  c( seq(1, 1.05*rel_max_470, 0.2),
     seq(1, 1.05*rel_max_525, 0.2),
     seq(1, 1.05*rel_max_625, 0.2)) %>%
  format(digits = 2)
xaxis_tick_labels_foot<-
  c( seq(1, 1.05*rel_max_430, 0.2),
     seq(1, 1.05*rel_max_525, 0.2),
     seq(1, 1.05*rel_max_545, 0.2)) %>%
  format(digits = 2)
xaxis_tick_labels_alfalfa<-
  c( seq(1, 1.05*rel_max_430, 0.2),
     seq(1, 1.05*rel_max_525, 0.2),
     seq(1, 1.05*rel_max_625, 0.2)) %>%
  format(digits = 2)

#Y limits and tick labels
ylim_large<-c(0, 900)
ylim_small<-c(0, 700)
yrange<-ylim_small[2]-ylim_small[1]
yaxis_ticks_large<-seq(0, ylim_large[2], 200)
yaxis_ticks_small<-seq(0, ylim_small[2], 200)

#Set figure width, margins and calculate plot sizes and figure height.
#All inputs in inches
f_width <- 7.5
m_bottom <- 1.0
m_left <- 0.5
m_top <- 0.29
m_right <- 0.05
subplot_height <- 1.81
subplot_width <- f_width-m_left-m_right
subplot_a_height <- subplot_height + m_top + m_bottom
subplot_b_height <- subplot_height + m_top + m_bottom
subplot_c_height <- ylim_large[2]/ylim_small[2] * subplot_height + m_top + m_bottom
subplot_d_height <- ylim_large[2]/ylim_small[2] * subplot_height + m_top + m_bottom
subplot_e_height <- ylim_large[2]/ylim_small[2] * subplot_height + m_top + m_bottom
f_height <- sum(subplot_a_height, subplot_b_height, subplot_c_height,
                subplot_d_height, subplot_e_height)
#calculates subplot y position
subplot_ypos <- c(1, 
                  1-(subplot_a_height) / f_height, 
                  1-(subplot_a_height + subplot_b_height) / f_height, 
                  1-(subplot_a_height + subplot_b_height + 
                       subplot_c_height) / f_height,
                  1-(subplot_a_height + subplot_b_height + 
                       subplot_c_height + subplot_d_height) / f_height, 
                  0)

axis_mgp <- c(3,0.35,0) #2nd value controls position of tick labels
tck1 <- -0.010 #Value to control tick marks length in a & b panels
tck2 <- -0.015 #Value to control tick marks length in c-e panels
n_trts<-length(levels(point.pred.CO2_VGrY$treatment_factor)) #number of treatments
CI_width <- 0.16*sp_btw_bars #half width of confidence interval rectangles in plot units
jitter <- 0.089*sp_btw_bars #half width of the jitter cloud
mean_lwd <- 2 #width of the mean line
#xaxis line positions
xaxis_pos <- c(xlim[1], 2.25*sp_btw_bars,
               2.75*sp_btw_bars + 0, 3.25*sp_btw_bars + 1,
               3.75*sp_btw_bars + 1, 4.25*sp_btw_bars + 2,
               4.75*sp_btw_bars + 2, 5.25*sp_btw_bars + 3,
               5.75*sp_btw_bars + 3, xlim[2])
#alternating ypos of sample size text
sample_size_ypos_large <- c(ylim_large[2]+0.02*yrange, ylim_large[2]+0.08*yrange)
sample_size_ypos_small <- c(ylim_small[2]+0.02*yrange, ylim_small[2]+0.08*yrange)
sample_size_cex <- 5/8 #size of sample size text
subplot_label_cex <- 12/8 #size of subplot levels
subplot_label_xpos <- 0.6 #x position of subplot labels in lines (0.2"/line)
subplot_label_ypos_large <- ylim_large[2] #y position of subplot labels in plot units
subplot_label_ypos_small <- ylim_small[2] #y position of subplot labels in plot units
cex.axis <- 6/8 #tick label font size
sample_size_col <- "gray70"
#space between the bars and the stars
sig_star_spacer_a <- 0.09*yrange
sig_star_spacer_b <- 0.11*yrange 
sig_star_spacer_c <- 0.10*yrange 
sig_star_spacer_d <- 0.08*yrange 
sig_star_spacer_e <- 0.10*yrange 

#odor diagram variables
od_xcenter <- -0.07*sp_btw_bars
od_xcenter_npc <- 
  (od_xcenter - (xlim[1]-m_left*xrange/subplot_width) ) / 
  (f_width*xrange/subplot_width)
od.x1 <- unit(od_xcenter_npc, "npc")
od.x2 <- unit(od_xcenter_npc-0.018, "npc")
od.x3 <- unit(od_xcenter_npc+0.029, "npc")
od.x4 <- unit(od_xcenter_npc+0.027, "npc")
od.x5 <- unit(od_xcenter_npc+0.022, "npc")
od.y.spacer <- -0.061
od.y1 <- unit(subplot_ypos[1] + od.y.spacer + 0.010, "npc")
od.y2 <- unit(subplot_ypos[2] + od.y.spacer + 0.010, "npc")
od.y3 <- unit(subplot_ypos[3] + od.y.spacer + 0.000, "npc")
od.y4 <- unit(subplot_ypos[3] + od.y.spacer + 0.005, "npc")
od.y5 <- unit(subplot_ypos[4] + od.y.spacer + 0.000, "npc")
od.y6 <- unit(subplot_ypos[4] + od.y.spacer + 0.000, "npc")
od.y7 <- unit(subplot_ypos[5] + od.y.spacer + 0.003, "npc")
od.y8 <- unit(subplot_ypos[5] + od.y.spacer + 0.000, "npc")
od.width <- unit(0.675, "inch")
od_label.x <- od_xcenter
od_label.y_small <- 467
od_label.y_large <- 600
od_label.cex <-1

#Variables for circle size and positioning
#incorporating yrange and ylim_small[1] keeps the circle position constant 
#when ylims change
trt_circle_y <- ylim_small[1]-0.300*yrange 
ctr_circle_y <- trt_circle_y-0.090*yrange #expressed relative to "trt_circle_y"
circle_text_x_spacer <- 0.12 #space between the circles and their labels plot units
circle_size<-3 #character expansion factor default=1
trt_circle_outline <- NA
ctr_circle_outline <- NA

#Variables for the positioning of the boxes around the circles
box_lwd <- 1.25
int_box_1_upper <- trt_circle_y+0.05*yrange
int_box_1_lower <- ctr_circle_y-0.05*yrange
int_box_1_left <- display.info.CO2$x.pos[5]-0.4*sp_btw_bars
int_box_1_right <- display.info.CO2$x.pos[12]+0.4*sp_btw_bars
int_box_2_upper <- trt_circle_y+0.05*yrange
int_box_2_lower <- ctr_circle_y-0.05*yrange
int_box_2_left <- display.info.CO2$x.pos[13]-0.4*sp_btw_bars
int_box_2_right <- display.info.CO2$x.pos[20]+0.4*sp_btw_bars
int_box_3_upper <- trt_circle_y+0.05*yrange
int_box_3_lower <- ctr_circle_y-0.05*yrange
int_box_3_left <- display.info.CO2$x.pos[21]-0.4*sp_btw_bars
int_box_3_right <- display.info.CO2$x.pos[28]+0.4*sp_btw_bars
odor_box_upper <- int_box_1_upper+0.016*yrange
odor_box_lower <- int_box_1_lower-0.016*yrange
odor_box_left <- display.info.alfalfa$x.pos[3]-0.4*sp_btw_bars
odor_box_right <- display.info.alfalfa$x.pos[28]+0.5*sp_btw_bars
clean_air_L_box_upper <- int_box_1_upper+0.016*yrange
clean_air_L_box_lower <- int_box_1_lower-0.016*yrange
clean_air_L_box_left <- display.info.alfalfa$x.pos[1]-0.35*sp_btw_bars
clean_air_L_box_right <- display.info.alfalfa$x.pos[2]+0.35*sp_btw_bars
clean_air_R_box_upper <- int_box_1_upper+0.016*yrange
clean_air_R_box_lower <- int_box_1_lower-0.016*yrange
clean_air_R_box_left <- display.info.alfalfa$x.pos[29]-0.35*sp_btw_bars
clean_air_R_box_right <- display.info.alfalfa$x.pos[29]+0.35*sp_btw_bars
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
box_col_gray <-  subset(display.info, trt_color=="gray" & trt_intensity=="1.00", 
                        select=trt.disp.col) %>% unique() %>%  as.character()
odor_box_col <- "gray70"
clean_air_box_col <- hsv(0.567, 0.50, 0.85)
box_label_y <- odor_box_lower-0.02*yrange



#-------- Execute code from here to end to generate graph

#Export figure to a pdf with the width and height from above
#also sets font type and size
pdf(file=paste0(output_dir,"figure S4 (activ).pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)

# #Export figure to a png with the width and height from above
# #also sets font type and size
# png(file=paste0(output_dir,"figure S4 (activ).png"),
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
  mai=c(m_bottom,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[2],subplot_ypos[1]),
  cex=1) #prevents font scaling in multi-panel figures

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim_small, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("A")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos_small, las=1, cex=subplot_label_cex)

#axis
axis(side=1, las=1, at=xaxis_ticks_CO2, labels=xaxis_tick_labels_CO2,
     mgp=axis_mgp, tck=tck1, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
lines(x=xaxis_pos[1:2], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[7:8], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[9:10], y=rep(ylim_small[1],2))
axis(side=2, las=1, at=yaxis_ticks_small, mgp=axis_mgp, tck=tck1, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim_small)
#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim_small, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim_small, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim_small, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim_small, col="grey90", lwd=1.5)

#axis labels
mtext("relative intensity", side=1, line=1.6)
mtext("mean number of trajectories", side=2, line=2.2)

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
stripchart(point_pred~treatment_factor, data=point.pred.CO2, 
           at=display.info.CO2$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.CO2$trt.disp.col.alpha)

#significance stars
text(x=display.info.CO2$x.pos,
     y=display.info.CO2$conf.high+sig_star_spacer_a,
     labels=display.info.CO2$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <-
  events %>% filter(odor == "CO2") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.CO2$x.pos, y=rep(sample_size_ypos_small,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#odor diagrams
grid.picture(CO2_SVG, x=od.x1, y=od.y1, just=c(0.5,0), width=od.width)
text(x=od_label.x, y=od_label.y_small, labels=bquote(bold(CO[2]~"alone")),
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
        y=c(int_box_1_upper, int_box_1_upper, int_box_1_lower, int_box_1_lower),
        col=NA, border=box_col_470, lwd=box_lwd)
polygon(x=c(int_box_2_left, int_box_2_right, int_box_2_right, int_box_2_left),
        y=c(int_box_2_upper, int_box_2_upper, int_box_2_lower, int_box_2_lower),
        col=NA, border=box_col_525, lwd=box_lwd)
polygon(x=c(int_box_3_left, int_box_3_right, int_box_3_right, int_box_3_left),
        y=c(int_box_3_upper, int_box_3_upper, int_box_3_lower, int_box_3_lower),
        col=NA, border=box_col_625, lwd=box_lwd)
polygon(x=c(odor_box_left, odor_box_right, odor_box_right, odor_box_left),
        y=c(odor_box_upper, odor_box_upper, odor_box_lower, odor_box_lower),
        col=NA, border=odor_box_col, lwd=box_lwd)
polygon(x=c(clean_air_L_box_left, clean_air_L_box_right, clean_air_L_box_right, clean_air_L_box_left),
        y=c(clean_air_L_box_upper, clean_air_L_box_upper, clean_air_L_box_lower, clean_air_L_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)
polygon(x=c(clean_air_R_box_left, clean_air_R_box_right, clean_air_R_box_right, clean_air_R_box_left),
        y=c(clean_air_R_box_upper, clean_air_R_box_upper, clean_air_R_box_lower, clean_air_R_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)

#box labels
text(x=c(int_box_1_right, int_box_2_right, int_box_3_right), y=box_label_y, adj=c(1,1),
     labels=c("478 nm intensity ramp","527 nm intensity ramp","621 nm intensity ramp"),
     col=c(box_col_470, box_col_525, box_col_625))
text(x=odor_box_left, y=box_label_y, adj=c(0,1), 
     labels=as.expression(bquote(CO[2]~" / odor period")), 
     col=odor_box_col)
text(x=clean_air_L_box_right, y=box_label_y, adj=c(1,1), 
     labels="clean air period", col=clean_air_box_col)


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
plot.window(xlim=xlim, ylim=ylim_small, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("B")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos_small, las=1, cex=subplot_label_cex)

#axis
axis(side=1, las=1, at=xaxis_ticks_CO2_VGrY, labels=xaxis_tick_labels_CO2_VGrY,
     mgp=axis_mgp, tck=tck1, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
lines(x=xaxis_pos[1:2], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[7:8], y=rep(ylim_small[1],2))
lines(x=xaxis_pos[9:10], y=rep(ylim_small[1],2))
axis(side=2, las=1, at=yaxis_ticks_small, mgp=axis_mgp, tck=tck1, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim_small)
#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim_small, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim_small, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim_small, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim_small, col="grey90", lwd=1.5)

#axis labels
mtext("relative intensity", side=1, line=1.6)
mtext("mean number of trajectories", side=2, line=2.2)

#for loop cycles through each set paired set of stimuli
for (i in seq(1,length(display.info.CO2_VGrY$x.pos),1)) {
  #confidence intervals
  polygon(x=c(-CI_width,CI_width,CI_width,-CI_width)
          + display.info.CO2_VGrY$x.pos[i],
          y=c(rep(display.info.CO2_VGrY$conf.low[i],2),
              rep(display.info.CO2_VGrY$conf.high[i],2)),
          col=display.info.CO2_VGrY$trt.disp.col.alpha[i],
          border=NA)
  #mean line
  lines(x=c(-CI_width,CI_width) + display.info.CO2_VGrY$x.pos[i],
        y=rep(display.info.CO2_VGrY$predicted[i],2),
        col=display.info.CO2_VGrY$trt.disp.col[i],
        lwd=mean_lwd, lend=1) # "lend=1" makes line ends flat
}

#set seed for jitter
set.seed(123456)

#stripchart
stripchart(point_pred~treatment_factor, data=point.pred.CO2_VGrY, 
           at=display.info.CO2_VGrY$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.CO2_VGrY$trt.disp.col.alpha)

#significance stars
text(x=display.info.CO2_VGrY$x.pos,
     y=display.info.CO2_VGrY$conf.high+sig_star_spacer_b,
     labels=display.info.CO2_VGrY$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <-
  events %>% filter(odor == "CO2_VGrY") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.CO2_VGrY$x.pos, y=rep(sample_size_ypos_small,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#odor diagrams
grid.picture(CO2_SVG, x=od.x1, y=od.y2, just=c(0.5,0), width=od.width)
text(x=od_label.x, y=od_label.y_small, labels=bquote(bold(CO[2]~"alone")),
     adj=c(0.5,0.5), cex=od_label.cex)

#Circle text labels
text(x=min(display.info.CO2_VGrY$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.CO2_VGrY$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.CO2_VGrY$x.pos, y=rep(trt_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.CO2_VGrY$trt.disp.col,
       col=NA)
points(x=display.info.CO2_VGrY$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.CO2_VGrY$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_1_left, int_box_1_right, int_box_1_right, int_box_1_left),
        y=c(int_box_1_upper, int_box_1_upper, int_box_1_lower, int_box_1_lower),
        col=NA, border=box_col_430, lwd=box_lwd)
polygon(x=c(int_box_2_left, int_box_2_right, int_box_2_right, int_box_2_left),
        y=c(int_box_2_upper, int_box_2_upper, int_box_2_lower, int_box_2_lower),
        col=NA, border=box_col_gray, lwd=box_lwd)
polygon(x=c(int_box_3_left, int_box_3_right, int_box_3_right, int_box_3_left),
        y=c(int_box_3_upper, int_box_3_upper, int_box_3_lower, int_box_3_lower),
        col=NA, border=box_col_545, lwd=box_lwd)
polygon(x=c(odor_box_left, odor_box_right, odor_box_right, odor_box_left),
        y=c(odor_box_upper, odor_box_upper, odor_box_lower, odor_box_lower),
        col=NA, border=odor_box_col, lwd=box_lwd)
polygon(x=c(clean_air_L_box_left, clean_air_L_box_right, clean_air_L_box_right, clean_air_L_box_left),
        y=c(clean_air_L_box_upper, clean_air_L_box_upper, clean_air_L_box_lower, clean_air_L_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)
polygon(x=c(clean_air_R_box_left, clean_air_R_box_right, clean_air_R_box_right, clean_air_R_box_left),
        y=c(clean_air_R_box_upper, clean_air_R_box_upper, clean_air_R_box_lower, clean_air_R_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)

#box labels
text(x=c(int_box_1_right, int_box_2_right, int_box_3_right), y=box_label_y, adj=c(1,1),
     labels=c("435 nm intensity ramp","gray intensity ramp","552 nm intensity ramp"),
     col=c(box_col_430, box_col_gray, box_col_545))
text(x=odor_box_left, y=box_label_y, adj=c(0,1), 
     labels=as.expression(bquote(CO[2]~" / odor period")), 
     col=odor_box_col)
text(x=clean_air_L_box_right, y=box_label_y, adj=c(1,1), 
     labels="clean air period", col=clean_air_box_col)



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
plot.window(xlim=xlim, ylim=ylim_large, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("C")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos_large, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks_tansy, labels=xaxis_tick_labels_tansy,
     mgp=axis_mgp, tck=tck2, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
lines(x=xaxis_pos[1:2], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[7:8], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[9:10], y=rep(ylim_large[1],2))
axis(side=2, las=1, at=yaxis_ticks_large, mgp=axis_mgp, tck=tck2, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim_large)
#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim_large, col="grey90", lwd=1.5)

#axis labels
mtext("relative intensity", side=1, line=1.6)
mtext("mean number of trajectories", side=2, line=2.2)

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
stripchart(point_pred~treatment_factor, data=point.pred.tansy, 
           at=display.info.tansy$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.tansy$trt.disp.col.alpha)

#significance stars
text(x=display.info.tansy$x.pos,
     y=display.info.tansy$conf.high+sig_star_spacer_c,
     labels=display.info.tansy$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <-
  events %>% filter(odor == "CO2 + tansy") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.tansy$x.pos, y=rep(sample_size_ypos_large,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y3, just=c(0.5,0), width=od.width)
grid.picture(tansy_SVG, x=od.x3, y=od.y4, just=c(0.5,0), height=od.width*0.8)
text(x=od_label.x, y=od_label.y_large, labels=bquote(bold("floral odor")),
     adj=c(0.5,0.5), cex=od_label.cex)

#Circle text labels
text(x=min(display.info.tansy$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.tansy$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.tansy$x.pos, y=rep(trt_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.tansy$trt.disp.col,
       col=NA)
points(x=display.info.tansy$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.tansy$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_1_left, int_box_1_right, int_box_1_right, int_box_1_left),
        y=c(int_box_1_upper, int_box_1_upper, int_box_1_lower, int_box_1_lower),
        col=NA, border=box_col_470, lwd=box_lwd)
polygon(x=c(int_box_2_left, int_box_2_right, int_box_2_right, int_box_2_left),
        y=c(int_box_2_upper, int_box_2_upper, int_box_2_lower, int_box_2_lower),
        col=NA, border=box_col_525, lwd=box_lwd)
polygon(x=c(int_box_3_left, int_box_3_right, int_box_3_right, int_box_3_left),
        y=c(int_box_3_upper, int_box_3_upper, int_box_3_lower, int_box_3_lower),
        col=NA, border=box_col_625, lwd=box_lwd)
polygon(x=c(odor_box_left, odor_box_right, odor_box_right, odor_box_left),
        y=c(odor_box_upper, odor_box_upper, odor_box_lower, odor_box_lower),
        col=NA, border=odor_box_col, lwd=box_lwd)
polygon(x=c(clean_air_L_box_left, clean_air_L_box_right, clean_air_L_box_right, clean_air_L_box_left),
        y=c(clean_air_L_box_upper, clean_air_L_box_upper, clean_air_L_box_lower, clean_air_L_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)
polygon(x=c(clean_air_R_box_left, clean_air_R_box_right, clean_air_R_box_right, clean_air_R_box_left),
        y=c(clean_air_R_box_upper, clean_air_R_box_upper, clean_air_R_box_lower, clean_air_R_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)

#box labels
text(x=c(int_box_1_right, int_box_2_right, int_box_3_right), y=box_label_y, adj=c(1,1),
     labels=c("478 nm intensity ramp","527 nm intensity ramp","621 nm intensity ramp"),
     col=c(box_col_470, box_col_525, box_col_625))
text(x=odor_box_left, y=box_label_y, adj=c(0,1), 
     labels=as.expression(bquote(CO[2]~" / odor period")), 
     col=odor_box_col)
text(x=clean_air_L_box_right, y=box_label_y, adj=c(1,1), 
     labels="clean air period", col=clean_air_box_col)



#---- Subplot d

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[5],subplot_ypos[4]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim_large, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("D")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos_large, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks_foot, labels=xaxis_tick_labels_foot,
     mgp=axis_mgp, tck=tck2, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
lines(x=xaxis_pos[1:2], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[7:8], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[9:10], y=rep(ylim_large[1],2))
axis(side=2, las=1, at=yaxis_ticks_large, mgp=axis_mgp, tck=tck2, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim_large)
#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim_large, col="grey90", lwd=1.5)

#axis labels
mtext("relative intensity", side=1, line=1.6)
mtext("mean number of trajectories", side=2, line=2.2)

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
stripchart(point_pred~treatment_factor, data=point.pred.foot, 
           at=display.info.foot$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.foot$trt.disp.col.alpha)

#significance stars
text(x=display.info.foot$x.pos,
     y=display.info.foot$conf.high+sig_star_spacer_d,
     labels=display.info.foot$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <-
  events %>% filter(odor == "CO2 + foot") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.foot$x.pos, y=rep(sample_size_ypos_large,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y5, just=c(0.5,0), width=od.width)
grid.picture(foot_SVG, x=od.x4, y=od.y6, just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y_large, labels=bquote(bold("host odor")),
     adj=c(0.5,0.5), cex=od_label.cex)

#Circle text labels
text(x=min(display.info.foot$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.foot$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.foot$x.pos, y=rep(trt_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.foot$trt.disp.col,
       col=NA)
points(x=display.info.foot$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.foot$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_1_left, int_box_1_right, int_box_1_right, int_box_1_left),
        y=c(int_box_1_upper, int_box_1_upper, int_box_1_lower, int_box_1_lower),
        col=NA, border=box_col_430, lwd=box_lwd)
polygon(x=c(int_box_2_left, int_box_2_right, int_box_2_right, int_box_2_left),
        y=c(int_box_2_upper, int_box_2_upper, int_box_2_lower, int_box_2_lower),
        col=NA, border=box_col_525, lwd=box_lwd)
polygon(x=c(int_box_3_left, int_box_3_right, int_box_3_right, int_box_3_left),
        y=c(int_box_3_upper, int_box_3_upper, int_box_3_lower, int_box_3_lower),
        col=NA, border=box_col_545, lwd=box_lwd)
polygon(x=c(odor_box_left, odor_box_right, odor_box_right, odor_box_left),
        y=c(odor_box_upper, odor_box_upper, odor_box_lower, odor_box_lower),
        col=NA, border=odor_box_col, lwd=box_lwd)
polygon(x=c(clean_air_L_box_left, clean_air_L_box_right, clean_air_L_box_right, clean_air_L_box_left),
        y=c(clean_air_L_box_upper, clean_air_L_box_upper, clean_air_L_box_lower, clean_air_L_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)
polygon(x=c(clean_air_R_box_left, clean_air_R_box_right, clean_air_R_box_right, clean_air_R_box_left),
        y=c(clean_air_R_box_upper, clean_air_R_box_upper, clean_air_R_box_lower, clean_air_R_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)

#box labels
text(x=c(int_box_1_right, int_box_2_right, int_box_3_right), y=box_label_y, adj=c(1,1),
     labels=c("435 nm intensity ramp","527 nm intensity ramp","552 nm intensity ramp"),
     col=c(box_col_430, box_col_525, box_col_545))
text(x=odor_box_left, y=box_label_y, adj=c(0,1), 
     labels=as.expression(bquote(CO[2]~" / odor period")), 
     col=odor_box_col)
text(x=clean_air_L_box_right, y=box_label_y, adj=c(1,1), 
     labels="clean air period", col=clean_air_box_col)



#---- Subplot e

par(
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[6],subplot_ypos[5]),
  #prevents font scaling in multi-panel figures
  cex=1, 
  new=T)

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim_large, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("E")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos_large, las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks_alfalfa, labels=xaxis_tick_labels_alfalfa,
     mgp=axis_mgp, tck=tck2, cex.axis=cex.axis, gap.axis = 0,
     col = NA, col.ticks = 1)
lines(x=xaxis_pos[1:2], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[7:8], y=rep(ylim_large[1],2))
lines(x=xaxis_pos[9:10], y=rep(ylim_large[1],2))
axis(side=2, las=1, at=yaxis_ticks_large, mgp=axis_mgp, tck=tck2, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim_large)
#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[6:7]),2), y=ylim_large, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[8:9]),2), y=ylim_large, col="grey90", lwd=1.5)

#axis labels
mtext("relative intensity", side=1, line=1.6)
mtext("mean number of trajectories", side=2, line=2.2)

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
stripchart(point_pred~treatment_factor, data=point.pred.alfalfa,
           at=display.info.alfalfa$x.pos, method = "jitter", jitter=jitter,
           vertical=TRUE, add=TRUE, pch=16,
           col=display.info.alfalfa$trt.disp.col.alpha)

#significance stars
text(x=display.info.alfalfa$x.pos,
     y=display.info.alfalfa$conf.high+sig_star_spacer_e,
     labels=display.info.alfalfa$sig_stars,
     adj=c(0.5,1))

#Sample size labels
labels <-
  events %>% filter(odor == "CO2 + alfalfa") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.alfalfa$x.pos, y=rep(sample_size_ypos_large,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y7, just=c(0.5,0), width=od.width)
grid.picture(alfalfa_SVG, x=od.x5, y=od.y8,just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y_large, labels=bquote(bold("plant infusion")),
     adj=c(0.5,0.5), cex=od_label.cex)

#Circle text labels
text(x=min(display.info.alfalfa$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.alfalfa$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.alfalfa$x.pos, y=rep(trt_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.alfalfa$trt.disp.col,
       col=NA)
points(x=display.info.alfalfa$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.alfalfa$ctr.disp.col,
       col=NA)

#Boxes around the circles
polygon(x=c(int_box_1_left, int_box_1_right, int_box_1_right, int_box_1_left),
        y=c(int_box_1_upper, int_box_1_upper, int_box_1_lower, int_box_1_lower),
        col=NA, border=box_col_430, lwd=box_lwd)
polygon(x=c(int_box_2_left, int_box_2_right, int_box_2_right, int_box_2_left),
        y=c(int_box_2_upper, int_box_2_upper, int_box_2_lower, int_box_2_lower),
        col=NA, border=box_col_525, lwd=box_lwd)
polygon(x=c(int_box_3_left, int_box_3_right, int_box_3_right, int_box_3_left),
        y=c(int_box_3_upper, int_box_3_upper, int_box_3_lower, int_box_3_lower),
        col=NA, border=box_col_625, lwd=box_lwd)
polygon(x=c(odor_box_left, odor_box_right, odor_box_right, odor_box_left),
        y=c(odor_box_upper, odor_box_upper, odor_box_lower, odor_box_lower),
        col=NA, border=odor_box_col, lwd=box_lwd)
polygon(x=c(clean_air_L_box_left, clean_air_L_box_right, clean_air_L_box_right, clean_air_L_box_left),
        y=c(clean_air_L_box_upper, clean_air_L_box_upper, clean_air_L_box_lower, clean_air_L_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)
polygon(x=c(clean_air_R_box_left, clean_air_R_box_right, clean_air_R_box_right, clean_air_R_box_left),
        y=c(clean_air_R_box_upper, clean_air_R_box_upper, clean_air_R_box_lower, clean_air_R_box_lower),
        col=NA, border=clean_air_box_col, lwd=box_lwd)

#box labels
text(x=c(int_box_1_right, int_box_2_right, int_box_3_right), y=box_label_y, adj=c(1,1),
     labels=c("435 nm intensity ramp","527 nm intensity ramp","621 nm intensity ramp"),
     col=c(box_col_430, box_col_525, box_col_625))
text(x=odor_box_left, y=box_label_y, adj=c(0,1), 
     labels=as.expression(bquote(CO[2]~" / odor period")), 
     col=odor_box_col)
text(x=clean_air_L_box_right, y=box_label_y, adj=c(1,1), 
     labels="clean air period", col=clean_air_box_col)

#Close Graph "Device"
dev.off()


