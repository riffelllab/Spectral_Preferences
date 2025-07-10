#Imports necessary packages
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(dplyr)
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
#output_dir<-"output files/color vs gray 0.50, air, CO2, water, CO2 + water/"
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

#Sets the sub directories to search for "resp_tracks.csv" track files
sub_dir_CO2<-"output files/color vs gray 0.50, CO2/"
sub_dir_CO2_tansy<-"output files/color vs gray 0.50, CO2 + tansy/"
sub_dir_CO2_foot<-"output files/color vs gray 0.50, CO2 + foot/"
sub_dir_CO2_alfalfa<-"output files/color vs gray 0.50, CO2 + alfalfa/"

#Extracts track level data from both groups of runs
events<-event_extract(sub_dir_CO2)
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
  mutate(odor = factor(odor, levels=c("CO2","CO2 + tansy",
                                      "CO2 + foot","CO2 + alfalfa"))) %>%
  #creates "treatment_factor"
  mutate(treatment_factor = case_when(
    #controls specified by individual stim_series values
    stim_series=="PreCO2Time" ~ "PreCO2Time",
    stim_series=="PreStim" ~ "PreStim",
    stim_series=="resp check, initial no odor" ~ "resp check, initial no odor",
    stim_series=="resp check, end of experiment" ~ "resp check, end of experiment",
    stim_series=="PostCO2Time" ~ "PostCO2Time",
    #end of experiment response check identified by intensity/color
    #of treatment and control as stim_series values has is different in 
    #current and older runs
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="1.00"  ~ "resp check, end of experiment",
    #black vs gray treatment assigned "black"
    trt_color=="all" & trt_intensity=="0.00" &
      ctr_color=="gray" & ctr_intensity=="0.50"  ~ "black",
    #gray treatment are assigned "gray-x.xx"
    trt_color=="gray" ~ paste0(trt_color,"-",trt_intensity),
    #all other assigned the trt color
    TRUE ~ trt_color), 
    .after = stim_series) %>%
  #Turns "treatment_factor" into a factor 
  mutate(treatment_factor = as.factor(treatment_factor)) %>%
  #Reorders the non color factor levels
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, 
                       "PreCO2Time", "resp check, initial no odor", "PreStim", 
                       "black", "gray-0.50", "gray-1.00", "385")) %>%
  mutate(treatment_factor = 
           fct_relevel(treatment_factor, 
                       "resp check, end of experiment", "PostCO2Time", after = Inf)) %>%
  #turns these variables numeric
  mutate_at(c("activity_index", "mean_duration", "mean_track_speed",
              "n_tracks", "n_tracks_resp", "p_tracks_resp", 
              "mean_pref_index",  "conf.low", "conf.high", "p_value"), 
            as.numeric) %>%
  #turns rep into a factor
  mutate(rep = as.factor(rep))
  

#Displays the total number of reps for each treatment_factor & odors
table(list(events$treatment_factor, events$odor))




#------------------------------------------------------
#------- mixed model analysis of recruitment ----------
#------------------------------------------------------

#Create two column vector of trajectories recruited to the stimuli and total trajectories
resp<-cbind(events$n_tracks_resp, 
            events$n_tracks)

#all stimuli model used for predictions
recruit.model <- glmmTMB(resp ~ 0 + odor + treatment_factor:odor + (1|rep), 
                        data = events, family = betabinomial(link = "logit"))
summary(recruit.model)

#Plot residuals to see if model fits
recruit.model_simres<-simulateResiduals(recruit.model)
plot(recruit.model_simres)
testOutliers(recruit.model_simres, type="bootstrap")

#Null Model
recruit.null <- glmmTMB(resp ~ (1|rep), data = events, 
                     family = betabinomial(link = "logit"))
summary(recruit.null)

#Run Anova for overall model fit
anova(recruit.model, recruit.null)



#---- Analysis with only spectral sweep
events.sub<-
  events %>%
  #Filters to only spectral sweeps
  filter(str_detect(stim_series, "spectral")) %>%
  mutate(treatment_factor = droplevels(treatment_factor))

#Create two column vector of trajectories recruited to the stimuli and total trajectories
resp<-cbind(events.sub$n_tracks_resp, 
            events.sub$n_tracks)

#Full Model
recruit.model.1 <- glmmTMB(resp ~ treatment_factor + #odor + 
                          treatment_factor:odor + (1|rep), 
                        data = events.sub, 
                        family = betabinomial(link = "logit"))
summary(recruit.model.1)

#Plot residuals to see if model fits
recruit.model.1_simres<-simulateResiduals(recruit.model.1)
plot(recruit.model.1_simres)
testOutliers(recruit.model.1_simres, type="bootstrap")


#no interaction model
recruit.model.2 <- glmmTMB(resp ~ treatment_factor + odor + (1|rep), 
                        data = events.sub, 
                        family = betabinomial(link = "logit"))
summary(recruit.model.2)

#anova testing the interaction
anova(recruit.model.1, recruit.model.2)


#no odor model
recruit.model.3 <- glmmTMB(resp ~ treatment_factor + (1|rep), 
                        data = events.sub, 
                        family = betabinomial(link = "logit"))
summary(recruit.model.3)

#anova testing the effect of odor
anova(recruit.model.1, recruit.model.3)


#no color model
recruit.model.4 <- glmmTMB(resp ~ odor + (1|rep), 
                        data = events.sub, 
                        family = betabinomial(link = "logit"))
summary(recruit.model.4)

#anova testing the effect of color
anova(recruit.model.1, recruit.model.4)


#null model
recruit.null <- glmmTMB(resp ~ (1|rep), 
                        data = events.sub, 
                        family = betabinomial(link = "logit"))
summary(recruit.null)

#anova testing overall model fit
anova(recruit.model.1, recruit.null)




#---------------------------------------------------------------
#------------ Extract recruitment model predictions ------------
#---------------------------------------------------------------

#creates list of unique combinations of "treatment_factor", "odor", "rep"
newdata <- events %>% select(treatment_factor, odor, rep) %>% distinct()

#point predictions
point.pred <-
  predict(recruit.model, re.form = NULL, type="response", newdata=newdata) %>%
  #combine with "newdata" frame
  cbind(newdata, .) %>%
  #rename prediction column
  rename(point_pred = ".")


#Extracts overall model predictions, SE and confidence intervals
model.pred<-
  predict_response(recruit.model, terms=c("odor","treatment_factor"), 
                   type="fixed", ci_level = 0.95) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(treatment_factor=group, odor=x) %>%
  #drops SE column
  select(!std.error) %>%
  #relocate rep column
  relocate(c(odor), .after=treatment_factor) %>%
  #sorts to match order of p-values in "recruit.model"
  arrange(treatment_factor, odor) %>%
  #extracts and adds the p-value from model summary
  mutate(p_value=summary(recruit.model)$coefficients$cond[,4]) %>%
  #puts things back to the original order
  arrange(odor, treatment_factor)




#------------------------------------------------------
#-------- Collate display info for the graph ----------
#------------------------------------------------------

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
  #----This section determines fine x positions from measured peak wavelengths
  #creates a column of wavelength (x) positions
  #creates NAs as intended, warning supressed
  mutate(x.pos=suppressWarnings(as.numeric(as.character(trt_color)))) %>%
  #adds a columns of peak wavelengths (from "functions.R")
  left_join(wavelengths, by=c("x.pos"="wl")) %>%
  #replace x.pos with peak_wl
  select(!x.pos) %>%
  rename(x.pos = peak_wl) %>%
  #Edits "x.pos" to place the non-color stimuli in the UV 
  #before the color stimuli
  mutate(x.pos=case_when(
    treatment_factor == "PreCO2Time" ~ 255,
    treatment_factor == "resp check, initial no odor" ~ 270,
    treatment_factor == "PreStim" ~ 300,
    treatment_factor == "black" ~ 330,
    treatment_factor == "gray-0.50" ~ 345,
    treatment_factor == "gray-1.00" ~ 360,
    treatment_factor == "resp check, end of experiment" ~ 780,
    treatment_factor == "PostCO2Time" ~ 810,
    TRUE ~ x.pos )) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=trt_color) %>%
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
display.info.tansy<-subset(display.info, odor=="CO2 + tansy")
display.info.foot<-subset(display.info, odor=="CO2 + foot")
display.info.alfalfa<-subset(display.info, odor=="CO2 + alfalfa")




#------------------------------------------------------
#--------- mixed model analysis of activation ---------
#------------------------------------------------------

#all stimuli model used for predictions
activ.model <- glmmTMB(n_tracks ~ odor + treatment_factor + treatment_factor:odor + (1|rep), 
                      data = events, family = nbinom1(link = "log"))
summary(activ.model)

#Plot residuals to see if model fits
activ.model_simres<-simulateResiduals(activ.model)
plot(activ.model_simres)
testOutliers(activ.model_simres, type="bootstrap")

#Null Model
activ.null <- glmmTMB(n_tracks ~ (1|rep), data = events, 
                     family = nbinom1(link = "log"))
summary(activ.null)

#Run Anova for overall model fit
anova(activ.model, activ.null)

#Calculate emmeans for all of the "treatment_factor:odor" treatment combinations
activ.model.emm<-regrid(emmeans(activ.model, specs = ~ treatment_factor:odor))

#specify different combination of treatment   
#combinations to compare with vectors of 1s and 0s
trt_lvls <- summary(activ.model.emm) %>% as.data.frame()
PreCO2Time <- as.numeric(trt_lvls$treatment_factor == "PreCO2Time" )
resp_check_init <- as.numeric(trt_lvls$treatment_factor == "resp check, initial no odor" )
PreStim <- as.numeric(trt_lvls$treatment_factor == "PreStim" )
stim_black <- as.numeric(trt_lvls$treatment_factor == "black" )
stim_gray_0.50 <- as.numeric(trt_lvls$treatment_factor == "gray-0.50" )
stim_gray_1.00 <- as.numeric(trt_lvls$treatment_factor == "gray-1.00" )
stim_430 <- as.numeric(trt_lvls$treatment_factor == "430" )
stim_525 <- as.numeric(trt_lvls$treatment_factor == "525" )
stim_625 <- as.numeric(trt_lvls$treatment_factor == "625" )
resp_check_end <- as.numeric(trt_lvls$treatment_factor == "resp check, end of experiment" )
PostCO2Time <- as.numeric(trt_lvls$treatment_factor == "PostCO2Time" )


#calculate the significance of the contrast between the different levels
contrasts<-
  contrast(activ.model.emm, adjust = "sidak", method = 
             list("resp_check_init - PreCO2Time"     = resp_check_init - PreCO2Time,
                  "PreStim         - PreCO2Time"     = PreStim         - PreCO2Time,
                  "stim_black      - PreStim"        = stim_black      - PreStim,
                  "stim_gray_0.50  - stim_black"     = stim_gray_0.50  - stim_black,
                  "stim_gray_1.00  - stim_gray_0.50" = stim_gray_1.00  - stim_gray_0.50,
                  "resp_check_end  - stim_black"     = resp_check_end  - stim_black,
                  "PostCO2Time     - PreCO2Time"     = PostCO2Time     - PreCO2Time,
                  "PostCO2Time     - resp_check_end" = resp_check_end  - stim_black
             )) %T>% print() %>%
  as.data.frame() %>%
  mutate(sig_stars=case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ "n.s." ))



#---- Analysis with only spectral sweep
events.sub<-
  events %>%
  #Filters to only spectral sweeps
  filter(str_detect(stim_series, "spectral")) %>%
  mutate(treatment_factor = droplevels(treatment_factor))

#Full Model
activ.model.1 <- glmmTMB(n_tracks ~ treatment_factor + odor + 
                          treatment_factor:odor + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.1)

#Plot residuals to see if model fits
activ.model.1_simres<-simulateResiduals(activ.model.1)
plot(activ.model.1_simres)
testOutliers(activ.model.1_simres, type="bootstrap")


#no interaction model
activ.model.2 <- glmmTMB(n_tracks ~ treatment_factor + odor + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.2)

#anova testing the interaction
anova(activ.model.1, activ.model.2)


#no odor model
activ.model.3 <- glmmTMB(n_tracks ~ treatment_factor + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.3)

#anova testing the effect of odor
anova(activ.model.1, activ.model.3)


#no color model
activ.model.4 <- glmmTMB(n_tracks ~ odor + (1|rep), 
                        data = events.sub, 
                        family = nbinom1(link = "log"))
summary(activ.model.4)

#anova testing the effect of color
anova(activ.model.1, activ.model.4)


#null model
activ.null <- glmmTMB(n_tracks ~ (1|rep), 
                     data = events.sub, 
                     family = nbinom1(link = "log"))
summary(activ.null)

#anova testing overall model fit
anova(activ.model.1, activ.null)




#--------------------------------------------------------------
#----------------- code to generate graph ---------------------
#--------------------------------------------------------------

#------ this section has variables used in the graph call below

#Set figure width, margins and calculate plot sizes and figure height.
#All inputs in inches
f_width <- 7.5
f_height <- 9.45
m_bottom <- 0.01
m_bottom_outer <- 1.0
m_left <- 0.5
m_top <- 0.29
m_right <- 0.05
subplot_width <- f_width-m_left-m_right
subplot_height <- (f_height - 4*m_top - 3*m_bottom - m_bottom_outer) / 4
subplot_1_height <- subplot_height + m_top + m_bottom
subplot_2_height <- subplot_height + m_top + m_bottom
subplot_3_height <- subplot_height + m_top + m_bottom
subplot_4_height <- subplot_height + m_top + m_bottom_outer
#calculates subplot y position
subplot_ypos <- c(1, 
                  1-(subplot_1_height)/f_height, 
                  1-(subplot_1_height+subplot_2_height)/f_height, 
                  1-(subplot_1_height+subplot_2_height+subplot_3_height)/f_height, 
                  0)

#X and Y limits
xlim<-c(min(display.info$x.pos)-30, 
        max(display.info$x.pos)+15)
xrange<-xlim[2]-xlim[1]
ylim<-c(0.00, 0.135)
yrange<-ylim[2]-ylim[1]
xaxis_ticks<-seq(400, 750, 50)
yaxis_ticks<-seq(0, ylim[2], 0.05)

axis_mgp <- c(3,0.35,0) #2nd value controls position of tick labels
tck <- -0.015 #Value to control tick mark length default is -0.01
n_trts<-length(levels(display.info$treatment_factor)) #number of treatment colors
CI_width <- 4.5 #half width of confidence interval rectangles in nm
jitter <- 2.25 #half width of the jitter cloud
mean_lwd <- 2 #width of the mean line
#xaxis line positions
xaxis_pos <- c(xlim[1], display.info.alfalfa$x.pos[6]+10,
               display.info.alfalfa$x.pos[7]-6, 750,
               display.info.alfalfa$x.pos[24]-7, display.info.alfalfa$x.pos[25]+8)
#alternating ypos of sample size text
sample_size_ypos <- c(ylim[2]+0.02*yrange, ylim[2]+0.08*yrange) 
sample_size_cex <- 5/8 #size of sample size text
subplot_label_cex <- 12/8 #size of subplot levels
subplot_label_xpos <- 0.6 #x position of subplot labels in lines (0.2"/line)
subplot_label_ypos <- ylim[2] #y position of subplot labels in plot units
cex.axis <- 6/8 #tick label font size
sample_size_col <- "gray70"
sig_star_spacer <- 0.10*yrange #space between the bars and the stars
sig_star_cex <- 1 #9.5/8 #size of significance stars

#odor diagram variables
od_xcenter <- 272
od_xcenter_npc <- 
  (od_xcenter - (xlim[1]-m_left*xrange/subplot_width) ) / 
  (f_width*xrange/subplot_width)
od.x1 <- unit(od_xcenter_npc, "npc")
od.x2 <- unit(od_xcenter_npc-0.020, "npc")
od.x3 <- unit(od_xcenter_npc+0.032, "npc")
od.x4 <- unit(od_xcenter_npc+0.030, "npc")
od.x5 <- unit(od_xcenter_npc+0.024, "npc")
od.y.spacer <- -0.11
od.y1 <- unit(subplot_ypos[1] + od.y.spacer + 0.000, "npc")
od.y2 <- unit(subplot_ypos[2] + od.y.spacer + 0.000, "npc")
od.y3 <- unit(subplot_ypos[2] + od.y.spacer + 0.005, "npc")
od.y4 <- unit(subplot_ypos[3] + od.y.spacer + 0.000, "npc")
od.y5 <- unit(subplot_ypos[3] + od.y.spacer + 0.000, "npc")
od.y6 <- unit(subplot_ypos[4] + od.y.spacer + 0.003, "npc")
od.y7 <- unit(subplot_ypos[4] + od.y.spacer - 0.000, "npc")
od.width <- unit(0.75, "inch")
od_label.x <- od_xcenter
od_label.y <- 0.082
od_label.cex <-1

#Variables for circle size and positioning
#incorporating yrange and ylim[1] keeps the circle position constant 
#when ylims change
trt_circle_y <- ylim[1]-0.30*yrange 
ctr_circle_y <- trt_circle_y-0.09*yrange #expressed relative to "trt_circle_y"
circle_text_x_spacer <- 15 #space between the circles and their labels in nm
circle_size<-3 #character expansion factor default=1
trt_circle_outline <- NA
ctr_circle_outline <- NA

#Variables for the positioning of the boxes around the circles
box_lwd <- 1.25
ss_box_upper <- trt_circle_y+0.05*yrange
ss_box_lower <- ctr_circle_y-0.05*yrange
ss_box_left <- display.info.alfalfa$x.pos[7]-10
ss_box_right <- display.info.alfalfa$x.pos[23]+10
odor_box_upper <- ss_box_upper+0.016*yrange
odor_box_lower <- ss_box_lower-0.016*yrange
odor_box_left <- display.info.alfalfa$x.pos[3]-10
odor_box_right <- display.info.alfalfa$x.pos[24]+10
clean_air_L_box_upper <- odor_box_upper
clean_air_L_box_lower <- odor_box_lower
clean_air_L_box_left <- display.info.alfalfa$x.pos[1]-10
clean_air_L_box_right <- display.info.alfalfa$x.pos[2]+10
clean_air_R_box_upper <- odor_box_upper
clean_air_R_box_lower <- odor_box_lower
clean_air_R_box_left <- display.info.alfalfa$x.pos[25]-10
clean_air_R_box_right <- display.info.alfalfa$x.pos[25]+10
ss_box_col <- "black"
odor_box_col <- "gray70"
clean_air_box_col <- hsv(0.567, 0.50, 0.85)
box_label_y <- odor_box_lower-0.02*yrange



#-------- Execute code from here to end to generate graph

# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure S3.pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)

#Export figure to a png with the width and height from above
#also sets font type and size
# png(file=paste0(output_dir,"figure S3.png"),
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
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("A")), side=2, line=subplot_label_xpos, 
      at=subplot_label_ypos, las=1, cex=subplot_label_cex)

#axis
lines(x=xaxis_pos[1:2], y=rep(ylim[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim[1],2))
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)
#gray lines dividing colors
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)

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
stripchart(point_pred~treatment_factor, data=subset(point.pred, odor=="CO2"), 
           at=display.info.CO2$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.CO2$trt.disp.col.alpha)

#significance stars
text(x=display.info.CO2$x.pos,
     y=display.info.CO2$conf.high+sig_star_spacer,
     labels=display.info.CO2$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  events %>% filter(odor == "CO2") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.CO2$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

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

#axis
lines(x=xaxis_pos[1:2], y=rep(ylim[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim[1],2))
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)

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
stripchart(point_pred~treatment_factor, data=subset(point.pred, odor=="CO2 + tansy"), 
           at=display.info.tansy$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.tansy$trt.disp.col.alpha)

#significance stars
text(x=display.info.tansy$x.pos,
     y=display.info.tansy$conf.high+sig_star_spacer,
     labels=display.info.tansy$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  events %>% filter(odor == "CO2 + tansy") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
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

#axis
lines(x=xaxis_pos[1:2], y=rep(ylim[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim[1],2))
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)

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
stripchart(point_pred~treatment_factor, data=subset(point.pred, odor=="CO2 + foot"), 
           at=display.info.foot$x.pos, method = "jitter", jitter=jitter, 
           vertical=TRUE, add=TRUE, pch=16, 
           col=display.info.foot$trt.disp.col.alpha)

#significance stars
text(x=display.info.foot$x.pos,
     y=display.info.foot$conf.high+sig_star_spacer,
     labels=display.info.foot$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#Sample size labels
labels <-
  events %>% filter(odor == "CO2 + foot") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.foot$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

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
axis(side=1, las=1, at=xaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=xaxis_pos[1:2], y=rep(ylim[1],2))
lines(x=xaxis_pos[3:4], y=rep(ylim[1],2))
lines(x=xaxis_pos[5:6], y=rep(ylim[1],2))
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=ylim)
lines(x=rep(mean(xaxis_pos[2:3]),2), y=ylim, col="grey90", lwd=1.5)
lines(x=rep(mean(xaxis_pos[4:5]),2), y=ylim, col="grey90", lwd=1.5)

#axis labels
mtext("wavelength (nm)", side=1, line=1.6)
mtext("proportion recruited", 
      side=2, line=2.2, at=0.305)

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
stripchart(point_pred~treatment_factor, data=subset(point.pred, odor=="CO2 + alfalfa"),
           at=display.info.alfalfa$x.pos, method = "jitter", jitter=jitter,
           vertical=TRUE, add=TRUE, pch=16,
           col=display.info.alfalfa$trt.disp.col.alpha)

#significance stars
text(x=display.info.alfalfa$x.pos,
     y=display.info.alfalfa$conf.high+sig_star_spacer,
     labels=display.info.alfalfa$sig_stars,
     adj=c(0.5,1), cex=sig_star_cex)

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y6, just=c(0.5,0), width=od.width)
grid.picture(alfalfa_SVG, x=od.x5, y=od.y7,just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("plant infusion")), 
     adj=c(0.5,0.5), cex=od_label.cex)

#Sample size labels
labels <-
  events %>% filter(odor == "CO2 + alfalfa") %>% 
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.alfalfa$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

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
polygon(x=c(ss_box_left, ss_box_right, ss_box_right, ss_box_left),
        y=c(ss_box_upper, ss_box_upper, ss_box_lower, ss_box_lower),
        col=NA, border=ss_box_col, lwd=box_lwd)
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
text(x=c(odor_box_left, ss_box_left), y=rep(box_label_y,2), adj=c(0,1),
     labels=as.expression(c(bquote(CO[2]~" / odor period"), "spectral sweep")),
     col=c(odor_box_col, ss_box_col))
text(x=clean_air_L_box_right, y=box_label_y, adj=c(1,1), 
     labels="clean air period", col=clean_air_box_col)

#Close Graph "Device"
dev.off()



