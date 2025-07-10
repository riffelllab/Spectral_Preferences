#Imports necessary packages
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(stringr)
library(forcats)
library(readr)
library(grImport2)
library(emmeans)
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
#-------- mixed model analysis of recruitment ---------
#------------------------------------------------------

#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(events$n_tracks_resp, 
            events$n_tracks)

#all stimuli model used for predictions
recruit.model <- glmmTMB(resp ~ odor + treatment_factor + treatment_factor:odor + (1|rep), 
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



#Calculate emmeans for all of the "treatment_factor:odor" treatment combinations
recruit.model.emm<-regrid(emmeans(recruit.model, specs = ~ treatment_factor:odor))

#specify different combination of treatment   
#combinations to compare with vectors of 1s and 0s
trt_lvls <- summary(recruit.model.emm) %>% as.data.frame()
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
  contrast(recruit.model.emm, adjust = "sidak", method = 
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




#---------------------------------------------------------------
#----------- Extract recruitment model predictions -------------
#---------------------------------------------------------------

#Extracts overall model predictions, SE and confidence intervals
#Prediction throws an expected error about interactions
recruit.model.pred<-
  #odor is omitted and predictions are averaged across all odors
  predict_response(recruit.model, terms=c("treatment_factor"), 
                   margin = "marginalmeans", type="fixed", ci_level = 0.95) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(treatment_factor=x, odor=group) %>%
  #drops SE column
  select(!c(odor, std.error)) %>%
  #sorts to match order of p-values in "recruit.model"
  arrange(treatment_factor)




#------------------------------------------------------
#-------- Collate display info for the graph ----------
#------------------------------------------------------

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info.recruit<-
  #uses tracks as input
  events %>%
  #limits to CO2 only runs as all odors are combined in this figure
  filter(odor=="CO2") %>%
  #----this section summarizes by trt_color and odor
  group_by(treatment_factor) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_color=first(trt_color),
            trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups="drop") %>%
  #sort by "odor" and then "trt_color"
  arrange(treatment_factor) %>%
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
  #creates a column of x positions
  mutate(x.pos=case_when(
    treatment_factor == "PreCO2Time" ~ 255,
    treatment_factor == "resp check, initial no odor" ~ 270,
    treatment_factor == "PreStim" ~ 300,
    treatment_factor == "black" ~ 330,
    treatment_factor == "gray-0.50" ~ 345,
    treatment_factor == "gray-1.00" ~ 360,
    treatment_factor == "430" ~ 390,
    treatment_factor == "525" ~ 405,
    treatment_factor == "625" ~ 420,
    treatment_factor == "resp check, end of experiment" ~ 450,
    treatment_factor == "PostCO2Time" ~ 480,
    TRUE ~ NA )) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=trt_color) %>%
  #----this section joins model predictions
  left_join(recruit.model.pred, by="treatment_factor") %>%
  filter(!is.na(x.pos))




#--------------------------------------------------------------
#----------------- code to generate graph ---------------------
#--------------------------------------------------------------

#------ this section has variables used in the graph call below

#Set figure width, margins and calculate plot sizes and figure height.
#All inputs in inches
f_width <- 3.5
f_height <- 3.5
m_bottom <- 0.88
m_left <- 0.5
m_top <- 0.15
m_right <- 0.05
subplot_width <- f_width-m_left-m_right
subplot_height <- f_height - m_top - m_bottom

#X and Y limits
xlim<-c(min(display.info.recruit$x.pos)-15, 
        max(display.info.recruit$x.pos)+15)
xrange<-xlim[2]-xlim[1]
ylim<-c(0.00, 0.135)
yrange<-ylim[2]-ylim[1]
xaxis_ticks<-seq(400, 750, 50)
yaxis_ticks<-seq(0, ylim[2], 0.05)

axis_mgp <- c(3,0.35,0) #2nd value controls position of tick labels
tck <- -0.015 #Value to control tick mark length default is -0.01
n_trts<-length(unique(display.info.recruit$treatment_factor)) #number of treatment colors
CI_width <- 6 #half width of confidence interval rectangles in nm
jitter <- 2.25 #half width of the jitter cloud
mean_lwd <- 2 #width of the mean line
sb.lwd <- 1.25 #width of the significance bars
#xaxis line positions
xaxis_pos <- c(xlim[1], display.info.recruit$x.pos[6]+10,
               display.info.recruit$x.pos[7]-6, 750,
               display.info.recruit$x.pos[24]-7, display.info.recruit$x.pos[25]+8)
#alternating ypos of sample size text
sample_size_ypos <- c(ylim[2]-0.04*yrange, ylim[2]+0.02*yrange) 
sample_size_cex <- 5/8 #size of sample size text
cex.axis <- 6/8 #tick label font size
sample_size_col <- "gray70"

#Variables for circle size and positioning
#incorporating yrange and ylim[1] keeps the circle position constant 
#when ylims change
trt_circle_y <- ylim[1]-0.15*yrange 
ctr_circle_y <- trt_circle_y-0.065*yrange #expressed relative to "trt_circle_y"
circle_text_x_spacer <- 15 #space between the circles and their labels in nm
circle_size<-3 #character expansion factor default=1
trt_circle_outline <- NA
ctr_circle_outline <- NA

#significance bars positions
sig_star_spacer <- 0.03*yrange #space between the bars and the stars
bar.1.y <- display.info.recruit$conf.high[1] + sig_star_spacer
bar.7.y <- bar.1.y + 2.2*sig_star_spacer
bar.3.y <- display.info.recruit$conf.high[4] + sig_star_spacer
bar.2.y <- bar.3.y - 2.2*sig_star_spacer
bar.4.y <- bar.3.y + 2.2*sig_star_spacer
bar.5.y <- bar.3.y - 2.2*sig_star_spacer
bar.6.y <- bar.4.y + 2.2*sig_star_spacer
bar.8.y <- bar.3.y - 2.2*sig_star_spacer
bar_vl.col <- "gray85"
bar_vl.lty <- "11"

#Variables for the positioning of the boxes around the circles
box_lwd <- 1.25
ss_box_upper <- trt_circle_y+0.035*yrange
ss_box_lower <- ctr_circle_y-0.035*yrange
ss_box_left <- display.info.recruit$x.pos[7]-10
ss_box_right <- display.info.recruit$x.pos[9]+10
odor_box_upper <- ss_box_upper+0.011*yrange
odor_box_lower <- ss_box_lower-0.011*yrange
odor_box_left <- display.info.recruit$x.pos[3]-10
odor_box_right <- display.info.recruit$x.pos[10]+10
clean_air_L_box_upper <- odor_box_upper
clean_air_L_box_lower <- odor_box_lower
clean_air_L_box_left <- display.info.recruit$x.pos[1]-10
clean_air_L_box_right <- display.info.recruit$x.pos[2]+10
clean_air_R_box_upper <- odor_box_upper
clean_air_R_box_lower <- odor_box_lower
clean_air_R_box_left <- display.info.recruit$x.pos[11]-10
clean_air_R_box_right <- display.info.recruit$x.pos[11]+10
ss_box_col <- "black"
odor_box_col <- "gray70"
clean_air_box_col <- hsv(0.567, 0.50, 0.85)
box_label_y <- odor_box_lower-0.013*yrange



#-------- Execute code from here to end to generate graph

# #Export figure to a pdf with the width and height from above
# #also sets font type and size
pdf(file=paste0(output_dir,"figure 3.pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)

#Export figure to a pdf with the width and height from above
#also sets font type and size
# png(file=paste0(output_dir,"figure 3.png"),
#     width=f_width, height=f_height, units="in", res=150,
#     family="sans", pointsize=8)

par(
  #specifies the outer (whole plot) margin in inches
  omi=c(0,0,0,0),
  #specifies the inner (subplot) margin in inches
  mai=c(m_bottom,m_left,m_top,m_right),
  #allows drawing outside of the plotting area
  xpd = TRUE,
  #prevents font scaling in multi-panel figures
  cex=1)

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#axes
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
box(bty="l")

#axis labels
mtext("wavelength (nm)", side=1, line=0.45)
mtext("proportion recruited", side=2, line=2.3)

#bar vertical lines
lines( x = rep(display.info.recruit$x.pos[1], 2), 
       y = c(display.info.recruit$conf.high[1], bar.2.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[2], 2), 
       y = c(display.info.recruit$conf.high[2], bar.1.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[3], 2), 
       y = c(display.info.recruit$conf.high[3], bar.3.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[4], 2), 
       y = c(display.info.recruit$conf.high[4], bar.6.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[5], 2), 
       y = c(display.info.recruit$conf.high[5], bar.4.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[6], 2), 
       y = c(display.info.recruit$conf.high[6], bar.5.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[10], 2), 
       y = c(display.info.recruit$conf.high[10], bar.6.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)
lines( x = rep(display.info.recruit$x.pos[11], 2), 
       y = c(display.info.recruit$conf.high[11], bar.8.y), 
       lwd = sb.lwd, lty = bar_vl.lty, col = bar_vl.col)

#significance bars
#bar 1
lines( x = c( rep(display.info.recruit$x.pos[1], 2), rep(display.info.recruit$x.pos[2], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.1.y, 
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(1,2)] ) + 2.5,
      y = bar.1.y  + 1*sig_star_spacer,
      labels = contrasts$sig_stars[1], adj = c(0.5,0.5) )
#bar 2
lines( x = c( rep(display.info.recruit$x.pos[1], 2), rep(display.info.recruit$x.pos[3], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.2.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(1,3)] ),
      y = bar.2.y  + 0.5*sig_star_spacer,
      labels = contrasts$sig_stars[2], adj = c(0.5,0.5) )
#bar 3
lines( x = c( rep(display.info.recruit$x.pos[3], 2), rep(display.info.recruit$x.pos[4], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.3.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(3,4)] ),
      y = bar.3.y  + 0.5*sig_star_spacer,
      labels = contrasts$sig_stars[3], adj = c(0.5,0.5) )
#bar 4
lines( x = c( rep(display.info.recruit$x.pos[4], 2), rep(display.info.recruit$x.pos[5], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.4.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(4,5)] ),
      y = bar.4.y  + 0.5*sig_star_spacer,
      labels = contrasts$sig_stars[4], adj = c(0.5,0.5) )
#bar 5
lines( x = c( rep(display.info.recruit$x.pos[5], 2), rep(display.info.recruit$x.pos[6], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.5.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(5,6)] ) + 2.5,
      y = bar.5.y  + 1*sig_star_spacer,
      labels = contrasts$sig_stars[5], adj = c(0.5,0.5) )
#bar 6
lines( x = c( rep(display.info.recruit$x.pos[4], 2), rep(display.info.recruit$x.pos[10], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.6.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(4,10)] ),
      y = bar.6.y  + 0.5*sig_star_spacer,
      labels = contrasts$sig_stars[6], adj = c(0.5,0.5) )
#bar 7
lines( x = c( rep(display.info.recruit$x.pos[1], 2), rep(display.info.recruit$x.pos[11], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.7.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(1,11)] ),
      y = bar.7.y  + 0.5*sig_star_spacer,
      labels = contrasts$sig_stars[7], adj = c(0.5,0.5) )
#bar 8
lines( x = c( rep(display.info.recruit$x.pos[10], 2), rep(display.info.recruit$x.pos[11], 2) ),
       y = c( -0.5*sig_star_spacer, 0, 0, -0.5*sig_star_spacer ) + bar.8.y,
       lwd = sb.lwd )
text( x = mean( display.info.recruit$x.pos[c(10,11)] ),
      y = bar.8.y  + 0.5*sig_star_spacer,
      labels = contrasts$sig_stars[6], adj = c(0.5,0.5) )

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

#Sample size labels
labels <-
  events %>%
  select(treatment_factor, n_tracks_resp) %>% 
  group_by(treatment_factor) %>% 
  #calculates the total number of recruited trajectories
  summarise(n_tracks_resp = sum(n_tracks_resp), .groups="drop") %>%
  #filters out treatment levels not used in "display.info.recruit"
  filter(treatment_factor %in% display.info.recruit$treatment_factor) %>%
  pull(n_tracks_resp) %>%
  paste0("(", ., ")" )
text(x=display.info.recruit$x.pos, y=rep(sample_size_ypos,length.out=n_trts),
     adj=c(0.5,1), labels=labels, cex=sample_size_cex, col=sample_size_col)

#Circle text labels
text(x=min(display.info.recruit$x.pos-circle_text_x_spacer), 
     y=trt_circle_y, adj=c(1,0.5), labels="test")
text(x=min(display.info.recruit$x.pos-circle_text_x_spacer), 
     y=ctr_circle_y, adj=c(1,0.5), labels="control")

#Treatment and control circles
points(x=display.info.recruit$x.pos, y=rep(trt_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.recruit$trt.disp.col,
       col=NA)
points(x=display.info.recruit$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.recruit$ctr.disp.col,
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



