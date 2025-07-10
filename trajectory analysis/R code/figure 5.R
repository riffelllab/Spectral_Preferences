#Imports necessary packages
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(emmeans)
library(dplyr)
library(magrittr)
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
tansy_SVG <- file.path(odor_source_dir, "tansy-cairo v2.svg") %>% readPicture()
foot_SVG <- file.path(odor_source_dir, "foot-cairo.svg") %>% readPicture()
alfalfa_SVG <- file.path(odor_source_dir, "plant infusion-cairo.svg") %>% readPicture()

#Sets the sub directories to search for "resp_tracks.csv" track files
sub_dir_CO2<-"output files/color vs gray 0.50, CO2/"
sub_dir_CO2_tansy<-"output files/color vs gray 0.50, CO2 + tansy/"
sub_dir_CO2_foot<-"output files/color vs gray 0.50, CO2 + foot/"
sub_dir_CO2_alfalfa<-"output files/color vs gray 0.50, CO2 + alfalfa/"

#Extracts track level data from each groups of runs
tracks<-track_extract(sub_dir_CO2)
tracks<-track_extract(sub_dir_CO2_tansy, tracks)
tracks<-track_extract(sub_dir_CO2_foot, tracks)
tracks<-track_extract(sub_dir_CO2_alfalfa, tracks)


#Creates a subset of data and processes it for analysis and display
tracks<-
  #starts with "tracks"
  tracks %>%
  #Filters to only spectral sweeps
  filter(str_detect(stim_series, "spectral")) %>%
  #Turns "trt_color" to a factor
  mutate(trt_color = as.factor(trt_color),
         rep = as.factor(rep)) %>%
  #turns these variables numeric
  mutate_at(c("pref_index", "duration", "trt_time", 
              "ctr_time", "stim_time"), as.numeric) %>%
  mutate(odor = factor(odor, levels=c("CO2","CO2 + tansy",
                                      "CO2 + foot","CO2 + alfalfa")))


#Displays the total number of trajectories for each treatment color
#and odor across all the runs
table.summary<-table(list(tracks$trt_color, tracks$rep))
#Displays the row and column means
round(as.data.frame(rowMeans(table.summary)), 0)
round(as.data.frame(colMeans(table.summary)), 0)




#----------------------------------------------------------
#---------- mixed model analysis of preference ------------
#----------------------------------------------------------

#Convert times to milliseconds and round to get the required integer inputs
resp<-cbind(round(tracks$trt_time*1000,0), 
            round(tracks$ctr_time*1000,0))


pref.model <- 
  glmmTMB(resp ~ 0 + trt_color:odor + (1|trt_color:rep), 
          data = tracks, family = betabinomial(link = "logit"))
summary(pref.model)

#Plot residuals to see if model fits
pref.model_simres<-simulateResiduals(pref.model)
plot(pref.model_simres)


#Odor color interaction model
pref.model.1 <- 
  glmmTMB(resp ~ 0 + trt_color + trt_color:odor + (1|trt_color:rep), 
          data = tracks, family = betabinomial(link = "logit"))
summary(pref.model.1)

#Plot residuals to see if model fits
pref.model.1_simres<-simulateResiduals(pref.model.1)
plot(pref.model.1_simres)


#Odor effect model
pref.model.2 <- 
  glmmTMB(resp ~ 0 + trt_color + odor + (1|trt_color:rep), 
          data = tracks, family = betabinomial(link = "logit"))
summary(pref.model.2)

#Anova for treatment effect of interaction
anova(pref.model.1, pref.model.2)


#No odor effect model
pref.model.3 <-
  glmmTMB(resp ~ 0 + trt_color + (1|trt_color:rep), 
          data = tracks, family = betabinomial(link = "logit"))
summary(pref.model.3)

#Anova for treatment effect of odor
anova(pref.model.1, pref.model.3)


#No trt_color effect model
pref.model.4 <-
  glmmTMB(resp ~ 0 + odor + (1|trt_color:rep), 
          data = tracks, family = betabinomial(link = "logit"))
summary(pref.model.4)

#Anova for treatment effect of color
anova(pref.model.1, pref.model.4)


#Null Model
pref.null <- glmmTMB(resp ~ 0 + (1|rep:trt_color), data = tracks, 
                        family = betabinomial(link = "logit"))
summary(pref.null)

#Anova for the overall model fit
anova(pref.model.1, pref.null)


#Calculate emmeans for all of the "trt_color:odor" treatment combinations
pref.model.emm<-regrid(emmeans(pref.model, specs = ~ trt_color:odor))

#specify different combination of treatment   
#combinations to compare with vectors of 1s and 0s
leds <- levels(tracks$trt_color) %>% as.numeric() 
blues <- (leds <= 450) %>% as.numeric()
greens<- (leds >= 470 & leds <= 570 ) %>% as.numeric()
reds <- (leds >= 590) %>% as.numeric()


blue_CO2      <- c(blues,     rep(0,17), rep(0,17), rep(0,17))
blue_tansy    <- c(rep(0,17), blues,     rep(0,17), rep(0,17))
blue_foot     <- c(rep(0,17), rep(0,17), blues,     rep(0,17))
blue_alfalfa  <- c(rep(0,17), rep(0,17), rep(0,17), blues    )
green_CO2     <- c(greens,    rep(0,17), rep(0,17), rep(0,17))
green_tansy   <- c(rep(0,17), greens,    rep(0,17), rep(0,17))
green_foot    <- c(rep(0,17), rep(0,17), greens,    rep(0,17))
green_alfalfa <- c(rep(0,17), rep(0,17), rep(0,17), greens   )
red_CO2       <- c(reds,      rep(0,17), rep(0,17), rep(0,17))
red_tansy     <- c(rep(0,17), reds,      rep(0,17), rep(0,17))
red_foot      <- c(rep(0,17), rep(0,17), reds,      rep(0,17))
red_alfalfa   <- c(rep(0,17), rep(0,17), rep(0,17), reds     )


#calculate the significance of the contrast between the different levels
color_contrast<-
  contrast(pref.model.emm, adjust = "none", method = 
             list("blue_tansy    - blue_CO2"  = blue_tansy    - blue_CO2,
                  "green_tansy   - green_CO2" = green_tansy   - green_CO2,
                  "red_tansy     - red_CO2"   = red_tansy     - red_CO2,
                  "blue_foot     - blue_CO2"  = blue_foot     - blue_CO2,
                  "green_foot    - green_CO2" = green_foot    - green_CO2,
                  "red_foot      - red_CO2"   = red_foot      - red_CO2,
                  "blue_alfalfa  - blue_CO2"  = blue_alfalfa  - blue_CO2,
                  "green_alfalfa - green_CO2" = green_alfalfa - green_CO2,
                  "red_alfalfa   - red_CO2"   = red_alfalfa   - red_CO2
             )) %T>% print() %>%
  as.data.frame() %>%
  mutate(sig_stars=case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    #p.value < 0.10  ~ ".",
    TRUE ~ "n.s." ))




#---------------------------------------------------------------
#------------ Extract preference model predictions -------------
#---------------------------------------------------------------

#Extracts overall model predictions, SE and confidence intervals
model.pred<-
  predict_response(pref.model.1, terms=c("odor","trt_color"), 
                   type="fixed", ci_level = 0.95) %>%
  #convert to dataframe
  as.data.frame() %>%
  #renames calculated columns to original terms
  rename(trt_color=group, odor=x) %>%
  #drops SE column
  select(!std.error) %>%
  #relocate rep column
  relocate(c(odor), .after=trt_color) %>%
  #converts from proportion to preference index
  mutate(across(!c(trt_color, odor), ~ 2 * .x - 1))




#------------------------------------------------------
#-------- Collate display info for the graph ----------
#------------------------------------------------------

#Create a list of treatment and control colors and 
#wavelength positions for graphic display
display.info<-
  #uses tracks as input
  tracks %>%
  #----this section summarizes by trt_color and odor
  group_by(trt_color, odor) %>%
  #takes the first value of the variables below when
  #summarizing all rows with the same "event" value
  summarise(trt_color=first(trt_color),
            trt_intensity=first(trt_intensity),
            ctr_color=first(ctr_color),
            ctr_intensity=first(ctr_intensity),
            .groups="drop") %>%
  #sort by "odor" and then "trt_color"
  arrange(odor, trt_color) %>%
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
  #----This section determines fine x positions from measured peak wavelengths
  #creates a column of wavelength (x) positions
  #creates NAs as intended, warning suppressed
  mutate(x.pos=suppressWarnings(as.numeric(as.character(trt_color)))) %>%
  #adds a columns of peak wavelengths (from "functions.R")
  left_join(wavelengths, by=c("x.pos"="wl")) %>%
  #replace x.pos with peak_wl
  select(!x.pos) %>%
  rename(x.pos = peak_wl) %>%
  #moves x.pos after treatment factor
  relocate(x.pos, .after=trt_color) %>%
  #joins model predictions
  left_join(model.pred, by=c("trt_color", "odor"))

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
f_width <- 3.5
f_height <- 7
m_bottom <- 0.9
m_left <- 0.4
m_top <- 0.2
m_right <- 0.025
subplot_width <- f_width-m_left-m_right
subplot_height <- (f_height - 3*m_top - m_bottom) / 3
subplot_1_height <- subplot_height + m_top
subplot_2_height <- subplot_height + m_top
subplot_3_height <- subplot_height + m_top + m_bottom
#calculates subplot y position
subplot_ypos <- c(1, 
                  1-(subplot_1_height)/f_height, 
                  1-(subplot_1_height+subplot_2_height)/f_height,
                  0)

#X and Y limits
xlim<-c(min(display.info$x.pos)-15, 760)
xrange<-xlim[2]-xlim[1]
ylim<-c(-0.015, 0.55)
yrange<-ylim[2]-ylim[1]
xaxis_ticks<-seq(400, 750, 50)
yaxis_ticks<-seq(0, ylim[2], 0.25)

axis_mgp <- c(3,0.35,0) #2nd value controls position of tick labels
tck <- -0.015 #Value to control tick mark length default is -0.01
n_trts<-length(levels(display.info$trt_color)) #number of treatment colors
sample_size_ypos <- c(ylim[2]+0.08*yrange, ylim[2]+0.02*yrange) #alternating ypos of sample size text
sample_size_cex <- 5/8 #size of sample size text
subplot_label_cex <- 12/8 #size of subplot levels
cex.axis <- 6/8 #tick label font size
sample_size_col <- "gray70"
sig_star_spacer <- 0.10*yrange #space between the bars and the stars
sb.pos <- yrange * 0.03

#colors
np_line_col <- "gray70"
CO2_col     <- "gray40"
tansy_col   <- rgb(255, 240,  10, maxColorValue = 255)
foot_col    <- rgb(200, 145, 105, maxColorValue = 255)
alfalfa_col <- rgb(100, 155,  60, maxColorValue = 255)

#line styles
CO2.lty     <- "32" #CO2 line type
np_line.lty <- 2 #no preference line type
np_line.lwd <- 1.5 #no preference line width
mean.lwd <- 2 #mean line width
CO2.lwd  <- 1.5 #CO2 line width
sb.lwd   <- 1.25 #significance bar line width
np_line.lend <- 1 #line end type "butt"
sb.lend  <- 1 #significance bar line end type "butt"

#odor diagram variables
od_xcenter <- 695
od_xcenter_npc <- 
  (od_xcenter - (xlim[1]-m_left*xrange/subplot_width) ) / 
  (f_width*xrange/subplot_width)
od.x2 <- unit(od_xcenter_npc-0.050, "npc")
od.x3 <- unit(od_xcenter_npc+0.055, "npc")
od.x4 <- unit(od_xcenter_npc+0.050, "npc")
od.x5 <- unit(od_xcenter_npc+0.044, "npc")
od.y.spacer <- -0.275
od.y2 <- unit(subplot_ypos[1] + od.y.spacer + 0.000, "npc")
od.y3 <- unit(subplot_ypos[1] + od.y.spacer + 0.005, "npc")
od.y4 <- unit(subplot_ypos[2] + od.y.spacer + 0.000, "npc")
od.y5 <- unit(subplot_ypos[2] + od.y.spacer + 0.000, "npc")
od.y6 <- unit(subplot_ypos[3] + od.y.spacer + 0.003, "npc")
od.y7 <- unit(subplot_ypos[3] + od.y.spacer - 0.000, "npc")
od.width <- unit(0.75, "inch")
od_label.x <- od_xcenter
od_label.y <- 0.03
od_label.cex <-1

#Variables for circle size and positioning
#incorporating yrange and ylim[1] keeps the circle position constant 
#when ylims change
trt_circle_y <- ylim[1]-0.27*yrange 
ctr_circle_y <- trt_circle_y-0.09*yrange #expressed relative to "trt_circle_y"
circle_text_x_spacer <- 20 #space between the circles and their labels in nm
circle_size<-2.5 #character expansion factor default=1
trt_circle_outline <- NA
ctr_circle_outline <- NA

#Variables for the positioning of the boxes around the circles
ss_box_upper <- trt_circle_y+0.045*yrange
ss_box_lower <- ctr_circle_y-0.045*yrange
ss_box_left <- display.info.alfalfa$x.pos[1]-10
ss_box_right <- display.info.alfalfa$x.pos[17]+10
ss_box_lwd <- 1.25
ss_box_col <- "black"
box_label_y <- ss_box_lower-0.02*yrange



#-------- Execute code from here to end to generate graph

#Export figure to a pdf with the width and height from above
#also sets font type and size
pdf(file=paste0(output_dir,"figure 5 (org).pdf"),
    width=f_width, height=f_height, family="sans", pointsize=8)

par(
  #specifies the outer (whole plot) margin in inches
  omi=c(0,0,0,0),
  #allows drawing outside of the plotting area
  xpd = TRUE)



#---- Subplot a

par(
  #specifies the inner (subplot) margin in inches
  mai=c(0,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[2],subplot_ypos[1]),
  cex=1) #prevents font scaling in multi-panel figures

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("A")), side=2, line=0.6, at=ylim[2], las=1, cex=subplot_label_cex)

#axis
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=c(0,ylim[2]))
#no preference line
lines(x=xlim, y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#mean polygon
polygon(x = c(display.info.CO2$x.pos, rev(display.info.CO2$x.pos)),
        y = c(display.info.CO2$predicted, rev(display.info.tansy$predicted)),
        border=NA, col="gray90")

#mean lines
lines(x = display.info.CO2$x.pos, display.info.CO2$predicted, 
      col=CO2_col, lwd=CO2.lwd, lty=CO2.lty)
lines(x = display.info.tansy$x.pos, display.info.tansy$predicted, 
      col=tansy_col, lwd=mean.lwd)

#significance bars
lines(x=display.info.CO2$x.pos[as.logical(blues)], 
      y=rep(ylim[2], sum(blues)), lwd=sb.lwd, lend=1)
lines(x=display.info.CO2$x.pos[as.logical(greens)], 
      y=rep(ylim[2], sum(greens)), lwd=sb.lwd, lend=1)
lines(x=display.info.CO2$x.pos[as.logical(reds)], 
      y=rep(ylim[2], sum(reds)), lwd=sb.lwd, lend=1)
#stars
text(x=c(mean(display.info.CO2$x.pos[as.logical(blues)]),
         mean(display.info.CO2$x.pos[as.logical(greens)]),
         mean(display.info.CO2$x.pos[as.logical(reds)])),
     y=ylim[2] + sb.pos*c(1,0.25,1), adj=c(0.5,0),
     labels=color_contrast$sig_stars[1:3])

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y2, just=c(0.5,0), width=od.width)
grid.picture(tansy_SVG, x=od.x3, y=od.y3, just=c(0.5,0), height=od.width*0.8)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("floral odor")), 
     adj=c(0.5,0.5), cex=od_label.cex)



#---- Subplot b

par(
  #specifies the inner (subplot) margin in inches
  mai=c(0,m_left,m_top,m_right),
  #specifies subplot position
  fig=c(0,1,subplot_ypos[3],subplot_ypos[2]),
  cex=1, #prevents font scaling in multi-panel figures
  new=T) #adds another subplot in combination with "fig"

#plots blank graph
plot.new()
#xaxs="i" & yaxs="i" remove extra 4% of space of either side of y and xlims
plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

#subtext label
mtext(bquote(bold("B")), side=2, line=0.6, at=ylim[2], las=1, cex=subplot_label_cex)

#axis
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=c(0,ylim[2]))
#no preference line
lines(x=xlim, y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#mean polygon
polygon(x = c(display.info.CO2$x.pos, rev(display.info.CO2$x.pos)),
        y = c(display.info.CO2$predicted, rev(display.info.foot$predicted)), 
        border=NA, col="gray90")

#mean lines
lines(x = display.info.CO2$x.pos, display.info.CO2$predicted, 
      col=CO2_col, lwd=CO2.lwd, lty=CO2.lty)
lines(x = display.info.foot$x.pos, display.info.foot$predicted, 
      col=foot_col, lwd=mean.lwd)

#significance bars
lines(x=display.info.CO2$x.pos[as.logical(blues)], 
      y=rep(ylim[2], sum(blues)), lwd=sb.lwd, lend=1)
lines(x=display.info.CO2$x.pos[as.logical(greens)], 
      y=rep(ylim[2], sum(greens)), lwd=sb.lwd, lend=1)
lines(x=display.info.CO2$x.pos[as.logical(reds)], 
      y=rep(ylim[2], sum(reds)), lwd=sb.lwd, lend=1)
#stars
text(x=c(mean(display.info.CO2$x.pos[as.logical(blues)]),
         mean(display.info.CO2$x.pos[as.logical(greens)]),
         mean(display.info.CO2$x.pos[as.logical(reds)])),
     y=ylim[2] + sb.pos*c(0.25,0.25,0.25), adj=c(0.5,0),
     labels=color_contrast$sig_stars[4:6])

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y4, just=c(0.5,0), width=od.width)
grid.picture(foot_SVG, x=od.x4, y=od.y5, just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("host odor")), 
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
mtext(bquote(bold("C")), side=2, line=0.6, at=ylim[2], las=1, cex=subplot_label_cex)

#axes
axis(side=1, las=1, at=xaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis,
     col = NA, col.ticks = 1)
axis(side=2, las=1, at=yaxis_ticks, mgp=axis_mgp, tck=tck, cex.axis=cex.axis)
lines(x=rep(xlim[1],2), y=c(0,ylim[2]))
#no preference line
lines(x=xlim, y=rep(0,2), col=np_line_col, 
      lty=np_line.lty, lwd=np_line.lwd, lend=np_line.lend)

#axis labels
mtext("wavelength (nm)", side=1, line=1.6)
mtext("preference index", side=2, line=2.0, at=0.93)

#mean polygon
polygon(x = c(display.info.CO2$x.pos, rev(display.info.CO2$x.pos)),
        y = c(display.info.CO2$predicted, rev(display.info.alfalfa$predicted)), 
        border=NA, col="gray90")

#mean lines
lines(x = display.info.CO2$x.pos, display.info.CO2$predicted, 
      col=CO2_col, lwd=CO2.lwd, lty=CO2.lty)
lines(x = display.info.alfalfa$x.pos, display.info.alfalfa$predicted, 
      col=alfalfa_col, lwd=mean.lwd)

#significance bars
lines(x=display.info.CO2$x.pos[as.logical(blues)], 
      y=rep(ylim[2], sum(blues)), lwd=sb.lwd, lend=1)
lines(x=display.info.CO2$x.pos[as.logical(greens)], 
      y=rep(ylim[2], sum(greens)), lwd=sb.lwd, lend=1)
lines(x=display.info.CO2$x.pos[as.logical(reds)], 
      y=rep(ylim[2], sum(reds)), lwd=sb.lwd, lend=1)
#stars
text(x=c(mean(display.info.CO2$x.pos[as.logical(blues)]),
         mean(display.info.CO2$x.pos[as.logical(greens)]),
         mean(display.info.CO2$x.pos[as.logical(reds)])),
     y=ylim[2] + sb.pos*c(1,0.25,1), adj=c(0.5,0),
     labels=color_contrast$sig_stars[7:9])

#odor diagrams
grid.picture(CO2_SVG, x=od.x2, y=od.y6, just=c(0.5,0), width=od.width)
grid.picture(alfalfa_SVG, x=od.x5, y=od.y7,just=c(0.5,0), height=od.width)
text(x=od_label.x, y=od_label.y, labels=bquote(bold("plant infusion")), 
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
       #circle outline color, NA = no outline
       col=trt_circle_outline)
points(x=display.info.alfalfa$x.pos, y=rep(ctr_circle_y, n_trts), 
       pch=21, cex=circle_size,
       #circle fill color
       bg=display.info.alfalfa$ctr.disp.col,
       #circle outline color, NA = no outline
       col=ctr_circle_outline)

#Boxes around the circles
polygon(x=c(ss_box_left, ss_box_right, ss_box_right, ss_box_left),
        y=c(ss_box_upper, ss_box_upper, ss_box_lower, ss_box_lower),
        col=NA, border=ss_box_col, lwd=ss_box_lwd)

#box labels
text(x=ss_box_left, y=box_label_y, adj=c(0,1),
     labels="spectral sweep", col=ss_box_col)

#Close Graph "Device"
dev.off()

