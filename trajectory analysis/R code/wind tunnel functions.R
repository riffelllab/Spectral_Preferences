#Loads required packages
library(dplyr)

#Loads custom color functions
source("R code/color functions.R")


#Function to read in trajectory positions (spots) from H5 files
h5read<-function(h5_loc) {
  #loads "hdf5r" package
  require(hdf5r)
  #Opens connection to H5 files
  file.h5<-H5File$new(h5_loc, mode = "r")
  #Creates R object for the "kalman_estimates" in the H5 file
  dataset<-file.h5[["kalman_estimates/"]]
  #Determines the row number / dimensions of the dataset
  dataset_dim<-dataset$dims
  #Uses this index to coerce "dataset" to a dataframe of that length
  dataset<-dataset[1:dataset_dim]
  #Closes connection to H5 files
  file.h5$close_all()
  #Keeps only "obj_id", frame","timestamp","x","y", and "z" columns
  dataset<-with(dataset,data.frame(obj_id,frame,timestamp,x,y,z))
  #Returns "dataset"
  dataset
}


#loads spectra data associating LED label wavelengths with 
#their measured peak wavelengths
wavelengths<-
  read.csv("../spectroscopy/irradiance/LED synth channels/led integrals full.csv", 
  colClasses = "numeric") %>%
  #drops other columns
  select(wl, peak_wl)


#Read in relative intensities of each color ramp.
#While the amount of light from the LEDs scales linearly
#with the setting (IE 0.50 or 1.00) this is displayed on top
#of the background illumination which is not spectrally flat
#for this reason the same LED intensity settings modifying
#isoquantal stimuli will give different relative intensities
#when background illumination is taken into account.
rel_intensity<-read.csv("R code/relative intensities.csv", colClasses="character")

  
#Function to produce a barplot of the "yvar" variable from "event_data_export"
event_barplot<-function(yvar, event_data, export_folder, ylab=NULL) {
  
  #loads ggplot for graph
  library(ggplot2)
  
  #subsets and processes the event data for the graph
  event_data_graph<-
    #Inputs "event_data"
    event_data %>%
    #adds a hex.color column for graph display
    mutate(hex.color = color_lookup(trt_color, trt_intensity)) %>%
    #adds an alpha (50%) to the color
    mutate(hex.color.alpha = color.add.alpha(hex.color, 70))
  
  #If "ylab" is not set us the column name from event_data_sub
  if(is.null(ylab)) { ylab<-enquo(yvar) }
  
  #Graph call with the various geometries
  bar.plot<-
    #Sets data and the variables defining the geometries below
    ggplot(event_data_graph, aes(x = event, y = {{yvar}}, fill = event)) +
    #Add a y axis label
    ylab(ylab) + 
    #sets the amount of extra space to leave around the data for the x and y axes
    scale_y_continuous(expand=expansion(mult = c(0, .10))) +
    scale_x_discrete(expand=expansion(add = c(0.8, 0.8))) +
    #creates facets for each stim_series (IE spectral downsweep)
    facet_grid(cols = vars(stim_series), scales = "free_x", space = "free_x") +
    #adds dotted line at 1 for the activity index graph
    {if (as_label(enquo(yvar))=="activity_index") 
      geom_hline(yintercept = 1, linetype = 1, color="grey80")} +
    #Lines added for "activity_index" barplot
    {if (as_label(enquo(yvar))=="activity_index" & 
         #Lines excluded if "PreCO2Time" is not an event
         "PreCO2Time" %in% event_data_graph$event)
      geom_hline(data=subset(event_data_graph, stim_series == "PreCO2Time"), 
                 aes(yintercept = mean(subset(event_data_graph, 
                            !(stim_series %in% c("PreCO2Time",
                                                 "resp check, initial no odor",
                                                 "PostCO2Time")),
                            activity_index, drop=TRUE))), linetype = 1, color="red")} +
    #adds dotted line at 5 for the number of tracks respondng graph
    {if (as_label(enquo(yvar))=="n_tracks_resp") 
      geom_hline(yintercept = 5, linetype = 1, color="grey80")} +
    #creates the bars
    geom_bar(stat = 'identity') +
    #sets the colors for the barplot
    scale_fill_manual(breaks = levels(event_data_graph$event), 
                      values = event_data_graph$hex.color.alpha)
    
  #Adds styling to graph
  bar.plot<-
    bar.plot +
    #sets base theme
    theme_minimal() +
    theme(
      #omits grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #omits x axis label
      axis.title.x = element_blank(),
      #rotates x-axis text
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      #adds margin for y axis title
      axis.title.y = element_text(margin = margin(t = 0, r = 0.4, b = 0, l = 0,
                                                  unit = "cm")),
      #adds additional white space above the plot
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
      #adds axis line "L"
      axis.line = element_line(colour = "black", linewidth = 0.25,
                               linetype = "solid", lineend="square"),
      #adds tick marks
      axis.ticks = element_line(colour = "black", linewidth = 0.25,
                                linetype = "solid", lineend="square"),
      axis.ticks.length = unit(0.15, unit = "cm"),
      #omits legend
      legend.position="none",
      #omits facet labels on top of the graph
      strip.text.x = element_blank() 
    )
  
  #exports to png, 
  #"as_label(enquo(" parses yvar as a character in the file name
  png(file.path(export_folder, paste0("barplot ",as_label(enquo(yvar)),".png")),
      width=12, height=6, units="in", res=150, pointsize=10, family="sans")
  #graph call, the call does need to be wrapped in a print command
  print(bar.plot)
  #Close Graph dev
  dev.off()
  
}


#Function to produce a violin plot of preference index
preference_index_violin_plot<-function(track_data, event_data,
                                       export_folder, ylab="preference index") {
  #loads ggplot for graph
  library(ggplot2)
  
  #Adds color and significance star columns
  event_data_graph<-
    #Inputs "event_data"
    event_data %>%
    #adds a hex.color column for graph display, variable name must be in "quotations"
    mutate(hex.color = color_lookup(trt_color, trt_intensity)) %>%
    #adds an alpha (50%) to the color
    mutate(hex.color.alpha = color.add.alpha(hex.color, 50)) %>%
    mutate(sig.stars = case_when(
      p_value >= 0.05 ~  "",
      p_value < 0.05 & p_value >= 0.01 ~ "*",
      p_value < 0.01 & p_value >= 0.001 ~ "**",
      p_value < 0.001 ~  "***"
    ))
  
  #Adds color and significance star columns
  track_data_graph<-
    #Inputs "event_data_export"
    track_data %>%
    #adds a hex.color column for graph display, variable name must be in "quotations"
    mutate(hex.color = color_lookup(trt_color, trt_intensity)) %>%
    #adds an alpha (0-30%) to the color
    mutate(hex.color.alpha = color.add.alpha(hex.color, 50))
  
  #scales the dotsize used depending on the max number of dots in 
  #a single event. This prevents overlap between events.
  dotsize<-40/max(event_data_graph$n_tracks_resp, na.rm = TRUE)
  if(dotsize > 10) {dotsize<-10}
  
  #Graph call with the various geometries
  violin.plot<-
    #Sets data and the variables defining the geometries below
    ggplot(event_data_graph, aes(x = event, y = pref_index, fill = event)) +
    #y limits & label
    ylim(-1.000,1.10) + ylab(ylab) +
    #creates facets for each stim_series (IE spectral downsweep)
    facet_grid(cols = vars(stim_series), scales = "free_x", space = "free_x") +
    #adds dotted line at 0
    geom_hline(yintercept = 0, linetype = 2, color="grey80") +
    #Add mean and confidence interval bar
    geom_crossbar(data=event_data_graph,
                  aes(x = event, y = mean_pref_index, 
                      ymin = conf.low, ymax = conf.high),
                  colour = event_data_graph$hex.color.alpha, fill = NA,
                  width = 0.25, linewidth = 0.5) +
    geom_text(data=event_data_graph, 
              aes(label = sig.stars, group = event, y = 1.00),
              vjust=0, hjust=0.5, size = 3) +
    geom_text(data=event_data_graph,
              aes(label = paste0("(",n_tracks_resp,")"), group = event, y = 1.10),
              vjust=0, hjust=0.5, size = 2.5, color = "grey70")+
    #violin plot
    geom_violin(data=track_data_graph, trim = TRUE, scale = "width", color = NA, 
                show.legend = FALSE) +
    #sets the colors for the violin plot
    scale_fill_manual(breaks = levels(track_data_graph$event), 
                      values = event_data_graph$hex.color.alpha) +
    #stacked dot plot
    geom_dotplot(
      data=track_data_graph,
      #multiplying the preference index by 0.99 allows the violin
      #plot to completely cover the points at 1 and -1.
      aes(x = event, y = pref_index*0.99),
      binaxis = "y", stackdir = "center", stackratio = 0.6,
      fill = track_data_graph$hex.color.alpha, stroke=NA, binwidth = 0.01, 
      dotsize=dotsize) 

  
  #Adds styling to graph
  violin.plot<-
    violin.plot +
    #sets base theme
    theme_minimal() +
    theme(
      #omits grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #omits x axis label
      axis.title.x = element_blank(),
      #rotates x-axis text
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      #adds margin for y axis title
      axis.title.y = element_text(margin = margin(t = 0, r = 0.4, b = 0, l = 0,
                                                  unit = "cm")),
      #adds axis line "L"
      axis.line = element_line(colour = "black", linewidth = 0.25,
                               linetype = "solid", lineend="square"),
      #adds tick marks
      axis.ticks = element_line(colour = "black", linewidth = 0.25,
                                linetype = "solid", lineend="square"),
      axis.ticks.length = unit(0.15, unit = "cm"),
      #omits legend
      legend.position="none",
      #omits facet labels on top of the graph
      strip.text.x = element_blank() 
    )
  
  #exports to png
  png(file.path(export_folder,"vioplot mean_pref_index.png"),
      width=12, height=6, units="in", res=150, pointsize=10, family="sans")
  
  #wrapping the graph call in the suppress warnings silences warning when
  #events have 1 or 0 trajectories
  suppressWarnings(print(
    violin.plot
  ))
  
  #Close Graph dev
  dev.off()
}


#Function to produce a heatmap of mosquito position from the "v_event"
event_heatmap<-function(v_event, spot_data, export_folder, 
                        cuesSetup, color.table, lim, cue_ann=TRUE) {
  
  #Creates a subset of "spot_data" from the specified "event"
  event_sub<-subset(spot_data, order==v_event)
  
  #--- Create a matrix of frequencies at all two dimensional xy coordinates
  #Set number of x bins
  x.nbins <- 601
  #Set number of y bins in proportion to the relative
  #lengths of "lim$x" and "lim$y"
  y.nbins <- ceiling(x.nbins*lim$y/lim$x)
  #Set number of y bins in proportion to the relative
  #lengths of "lim$x" and "lim$Z"
  z.nbins <- ceiling(x.nbins*lim$z/(2*lim$x))
  #Calculate x and y bins
  x.bin <- seq(-lim$x, lim$x, length=x.nbins)
  y.bin <- seq(-lim$y, lim$y, length=y.nbins)
  z.bin <- seq(     0, lim$z, length=z.nbins)
  #Creates a list of xy & xz frequencies
  #findInterval lists the bin index for each x/y/z value
  #table then cross tabulates the x and y or x and z frequencies
  xy_freq <-  as.data.frame(table( findInterval(event_sub$x, x.bin), 
                                   findInterval(event_sub$y, y.bin) ))
  xz_freq <-  as.data.frame(table( findInterval(event_sub$x, x.bin), 
                                   findInterval(event_sub$z, z.bin) ))
  #Converts from frequency to probability density / occupancy
  xy_freq$Freq<-xy_freq$Freq/sum(xy_freq$Freq)
  xz_freq$Freq<-xz_freq$Freq/sum(xz_freq$Freq)
  #Converts the factors back to numeric data
  xy_freq[,1] <- as.numeric(xy_freq[,1])
  xy_freq[,2] <- as.numeric(xy_freq[,2])
  xz_freq[,1] <- as.numeric(xz_freq[,1])
  xz_freq[,2] <- as.numeric(xz_freq[,2])
  #Creates an empty matrix to write to
  xy_prob2D <- matrix(data=0, nrow=x.nbins, ncol=y.nbins)
  xz_prob2D <- matrix(data=0, nrow=x.nbins, ncol=z.nbins)
  #Folds the table up into the matrix
  xy_prob2D[cbind(xy_freq[,1], xy_freq[,2])] <- xy_freq[,3]
  xz_prob2D[cbind(xz_freq[,1], xz_freq[,2])] <- xz_freq[,3]
  #Sets the min and max for for the color ramp breaks
  cr_lim<-c(0,0.0001)
  #Sets a series of breaks for the color ramp
  #All values greater 0.01% probability density will get max color value
  xy_breaks<-c(seq(cr_lim[1],cr_lim[2],cr_lim[2]/10),max(xy_prob2D))
  xz_breaks<-c(seq(cr_lim[1],cr_lim[2],cr_lim[2]/10),max(xz_prob2D))
  
  
  #Set figure width and calculate plot sizes and figure height
  #All inputs in inches
  f_width<-7.5
  m_bottom<-0.5
  m_left<-0.6
  m_top<-0.1
  m_right<-0.55
  m_inner<-0.3
  m_legend<-0.15
  #Relative widths of main plot and the legend plot (sum to 1)
  main_pwidth<-0.94
  legend_pwidth<-1-main_pwidth
  #width of mainplot
  plot_width<-(f_width-m_left-2*m_legend-m_right)*main_pwidth
  #scales plot size to maintain xy aspect ratio of mainplot
  plot_height<-plot_width*(y.nbins/x.nbins)
  f_height<-2*plot_height+2*m_inner+m_bottom+m_top
  
  #tick length
  tcl<-(-0.4)
  #tick label position
  mgp.x<-c(3, 0.4, 0)
  mgp.y<-c(3, 0.6, 0)
  #axis label margin position in "lines" (1 line = 0.2")
  m.line.x<-2
  m.line.y<-2.7
  #legend label x position
  legend.x<-5.5
  #line width for border and axes
  lwd<-0.75
  #legend text size
  leg.cex=0.9
  #label text size
  label.m.cex=1.1
  label.t.cex=1.3
  label.t2.cex=1.1
  
  #Creates a color ramp for heatmaps using ColorBrewer
  rf <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,'GnBu')))
  r <- rf(length(xy_breaks)-1)
  #Creates a 2x121 matrix of pixel values for the matrix
  #The matrix is >100 in the long dimension to give a smooth gradient
  #The matrix goes beyond the limit of the color ramp breaks in order to
  #give space for the max value of the probability matrix which is 
  #plotted at occupancy value at 0.012%
  leg_breaks<-
    matrix( rep( seq(cr_lim[1], cr_lim[2]*1.2,  length.out=121), 2), 
            nrow=2, byrow=TRUE)
  #Color ramp for legend. It stops at 0.01% so that all values beyond that show
  #the color at the top of the ramp. The -1 is because there needs to be one less
  #color than the number of breaks
  leg_r <- rf(length(leg_breaks[1,1:101])-1)
  
  #---Sets stimulus colors
  #Sets odor color to 50% grey
  odor.color<-hsv(0,0,0.5,0.7)
  
  #Looks up "a" and "b" colors from "color.table" defined with functions
  #The "[1]" index uses the first row for the look up.
  #a 70% alpha is then added
  a.color <- 
    color_lookup(event_sub$a_color[1], event_sub$a_intensity[1]) %>%
    color.add.alpha(70)
  b.color <- 
    color_lookup(event_sub$b_color[1], event_sub$b_intensity[1]) %>%
    color.add.alpha(70)
  
  #Sets file path and file name for pdf 
  filepath<-file.path(export_folder,
                      paste0("heatmap-",event_sub$order[1],
                             " ",event_sub$event[1],".png" ))
  #Export figure to a pdf with the width and height from above
  png(file=filepath, width=f_width, height=f_height, 
      units="in", res=300, type="cairo", pointsize=10, family="sans")
  
  #Margins
  par(mai=c(0,0,m_inner,m_legend), #inner margins in inches
      omi=c(m_bottom,m_left,m_top,m_right), #outer margins in inches
      cex=1) #Ensures the font size stays the same with multiple plots
  
  #Creates a layout of the four plots, with the two heat maps
  #taking up 94% of the width, and the two legends taking up 6%
  layout(matrix(c(1:4),2,2, byrow=T), 
         width=c(main_pwidth, legend_pwidth, main_pwidth, legend_pwidth))
  
  
  #Creates xy image from matrix, each cell of matrix is one pixel
  image(x.bin, y.bin, xy_prob2D, zlim=cr_lim, breaks=xy_breaks, col=r,
        xaxt= "n", yaxt= "n", xlab="", ylab="") #omits axes
  #------------------------------------------------------
  #Only added if cue annotations are turned on
  if(cue_ann==TRUE){
  #Labels C02 source
  points(x=cuesSetup["odor","x"], y=cuesSetup["odor","y"], 
         pch=21, col="black", bg=odor.color, lwd=1, cex=3)
  #Labels LED stimuli, creates a square polygon with a black border, and
  #a semi-transparent color from the "color.table". The position and size
  #are identified above.
  plotrix::draw.circle(x=cuesSetup["a","x"],
                       y=cuesSetup["a","y"],
                       radius=cuesSetup["a","size"],
                       border="black", col=a.color, lwd=1)
  plotrix::draw.circle(x=cuesSetup["b","x"],
                       y=cuesSetup["b","y"],
                       radius=cuesSetup["b","size"],
                       border="black", col=b.color, lwd=1)
  }#-------------------------------------------------------
  #Y axis
  axis(2, at=round(seq(-0.3, 0.3, 0.1),2), mgp=mgp.y, tcl=tcl, las=1,
       lwd=lwd, cex.axis=leg.cex)
  #top view label
  text(x=-0.91, y=0.32, "top view", adj=c(0,0), cex=label.t2.cex, xpd=NA)
  #Y axis label
  mtext("Y (m)", side=2, line=m.line.y, cex=label.m.cex)
  #draws border
  box(lwd=lwd)
  
  
  #Draws the color gradient for the legend
  image(c(0,1), leg_breaks[1,], leg_breaks, zlim=cr_lim, 
        ylim=c(cr_lim[1],cr_lim[2]*1.2), 
        breaks=c(leg_breaks[1,1:100],max(xy_prob2D)), col=leg_r,
        xaxt= "n", yaxt= "n", xlab="", ylab="", #omits axes
        bty="n") #omits border
  #Adds axis for the legend
  #Legend shows the max occupancy value at 0.012% with eveything above 0.01%
  #being colored with the color from the top of the color ramp
  axis(4, at=seq(cr_lim[1], cr_lim[2]*1.2, cr_lim[2]*0.2), 
       labels=sprintf("%.3f", 100*c(
         seq(cr_lim[1], cr_lim[2], cr_lim[2]*0.2),
         max(xy_prob2D) ) ),
       mgp=mgp.y, tcl=tcl, las=1, lwd=lwd, cex.axis=leg.cex)
  #Axis break lines
  lines(x=c(1.3,1.7),y=c(cr_lim[2]*1.10,cr_lim[2]*1.13), lwd=lwd, xpd=NA)
  lines(x=c(1.3,1.7),y=c(cr_lim[2]*1.07,cr_lim[2]*1.10), lwd=lwd, xpd=NA)
  #Axis label for the graph, text is used here to allow for 270° rotation
  #xpd=NA allows this text to be drawn anywhere on the figure
  text(legend.x, y=cr_lim[2]/2, "occupancy (%)", adj=c(0.5,0), 
       cex=label.t.cex, srt=270, xpd=NA)
  #draws border
  box(lwd=lwd)
  
  
  #Creates xz image from matrix, each cell of matrix is one pixel
  image(x.bin, z.bin, xz_prob2D, zlim=cr_lim, breaks=xz_breaks, col=r,
        xaxt= "n", yaxt= "n", xlab="", ylab="")  #omits axes
  #------------------------------------------------------
  #Only added if cue annotations are turned on
  if(cue_ann==TRUE){
  #Labels C02 source
  points(x=cuesSetup["odor","x"], y=cuesSetup["odor","z"], 
         pch=21, col="black", bg=odor.color, lwd=1, cex=3)
  #Labels LED stimuli, creates a square polygon with a black border, and
  #a semi-transparent color from the "color.table". The position and size
  #are identified above. Only the a stimulus is shown in the side view
  polygon(x=c(cuesSetup["a","x"]-cuesSetup["a","size"], 
              cuesSetup["a","x"]+cuesSetup["a","size"],
              cuesSetup["a","x"]+cuesSetup["a","size"], 
              cuesSetup["a","x"]-cuesSetup["a","size"]),
          #stimuli have no size in Z dimension but are depicted here 2 cm tall
          #centered on the x=axis, "xpd=NA" allows the polygon to be drawn
          #outside the plot
          y=c(cuesSetup["a","z"], cuesSetup["a","z"], 
              -cuesSetup["a","z"], -cuesSetup["a","z"]), 
          border="black", col=a.color, lwd=1, xpd=NA)
  }#-------------------------------------------------------
  #Y(Z) axis
  axis(2, at=seq(0, 0.6, 0.1), mgp=mgp.y, tcl=tcl, las=1,
       lwd=lwd, cex.axis=leg.cex)
  #X axis 
  axis(1, at=seq(-0.9, 0.9, 0.1), mgp=mgp.x, tcl=tcl,
       lwd=lwd, cex.axis=leg.cex)
  #side view label
  text(x=-0.91, y=0.625, "side view", adj=c(0,0), cex=label.t2.cex, xpd=NA)
  #Y(Z) axis label
  mtext("Z (m)", side=2, line=m.line.y, cex=label.m.cex)
  #X axis label
  mtext("X (m)", side=1, line=m.line.x, cex=label.m.cex)
  #Airflow label
  text(x=-0.92, y=-0.11, "airflow", adj=c(0,0.5), cex=label.t2.cex, xpd=NA)
  #Airflow arrow
  lines(x=c(-0.79,-0.71), y=c(-0.114, -0.114), lwd=1.5, xpd=NA)
  lines(x=c(-0.72, -0.71, -0.72), y=c(-0.104, -0.114, -0.124), lwd=1.5, xpd=NA)
  #draws border
  box(lwd=lwd)
  
  
  #Draws the color gradient for the legend
  image(c(0,1), leg_breaks[1,], leg_breaks, zlim=cr_lim, 
        ylim=c(cr_lim[1],cr_lim[2]*1.2),
        breaks=c(leg_breaks[1,1:100],max(xz_prob2D)), col=leg_r, 
        xaxt= "n", yaxt= "n", xlab="", ylab="", #omits axes
        bty="n") #omits border
  #Adds axis for the legend
  #Legend shows the max occupancy value at 0.012% with eveything above 0.01%
  #being colored with the color from the top of the color ramp
  axis(4, at=seq(cr_lim[1], cr_lim[2]*1.2, cr_lim[2]*0.2), 
       labels=sprintf("%.3f", 100*c(
         seq(cr_lim[1], cr_lim[2], cr_lim[2]*0.2),
         max(xz_prob2D) ) ),
       mgp=mgp.y, tcl=tcl, las=1, lwd=lwd, cex.axis=leg.cex)
  #Axis break lines
  lines(x=c(1.3,1.7),y=c(cr_lim[2]*1.10,cr_lim[2]*1.13), lwd=lwd, xpd=NA)
  lines(x=c(1.3,1.7),y=c(cr_lim[2]*1.07,cr_lim[2]*1.10), lwd=lwd, xpd=NA)
  #Axis label for the graph, text is used here to allow for 270° rotation
  #xpd=NA allows this text to be drawn anywhere on the figure
  text(legend.x, y=cr_lim[2]/2, "occupancy (%)", adj=c(0.5,0), 
       cex=label.t.cex, srt=270, xpd=NA)
  #draws border
  box(lwd=lwd) 
  
  #Concatenates title from values in the first row of "event_sub"
  title<-paste0(event_sub$stim_series[1],", odor: ",event_sub$odor[1],
                ", a: ",event_sub$a_color[1]," nm - ",event_sub$a_intensity[1],
                " vs b: ",event_sub$b_color[1]," nm - ",event_sub$b_intensity[1])
  #Removes "nm" from titles with "all"
  title<-gsub(pattern="all nm", replacement="all", x=title)
  
  #main title for the figure, add more descriptive label
  mtext(title, side=3, outer=TRUE, line=-1.4, cex=label.m.cex)
  
  #Close Graph dev
  dev.off()
}

#Function to extract track level data from multiple runs and append it to
#a dataframe. For one group (IE CO2) of runs a null dataframe is used.
#For multiple groups (IE CO2 and CO2 + Tansy) the function can be run twice with
#The output from the first run used as the input for the second run.
track_extract<-function(sub_dir, tracks_out=NULL) {
  
  #---- Track level files
  #Creates a list of "tracks data (resp).csv" track files
  track.files<-
    list.files(sub_dir, pattern="tracks data \\(resp\\).csv", recursive = TRUE)
  #Filters out track files moved into the "~exclude" folder
  track.files<-
    grep("~exclude", track.files, value=T, invert=T)
  
  #Reads in and rbinds each csv file
  for(i in seq(1, length(track.files), 1)){
    tracks_i<-
      readr::read_csv(paste0(sub_dir,track.files[i]), 
               #keeps all columns as characters
               col_types = readr::cols(.default = "c"))
    tracks_out<-bind_rows(tracks_out, tracks_i)
    remove(tracks_i)
  }
  #output
  tracks_out
}

#Function to extract vent level data from multiple runs and append it to
#a dataframe. For one group (IE CO2) of runs a null dataframe is used.
#For multiple groups (IE CO2 and CO2 + Tansy) the function can be run twice with
#The output from the first run used as the input for the second run.
event_extract<-function(sub_dir, events_out=NULL) {
  
  #---- event level files
  #Creates a list of "event data.csv" event files
  event.files<-
    list.files(sub_dir, pattern="event data.csv", recursive = TRUE)
  #Filters out event files moved into the "~exclude" folder
  event.files<-
    grep("~exclude", event.files, value=T, invert=T)
  
  #Reads in and rbinds each csv file
  for(i in seq(1, length(event.files), 1)){
    events_i<-
      readr::read_csv(paste0(sub_dir,event.files[i]), 
               #keeps all columns as characters
               col_types = readr::cols(.default = "c"))
    events_out<-bind_rows(events_out, events_i)
    remove(events_i)
  }
  #output
  events_out
}

#Function to extract treatment level data from multiple experiments and append 
#it to a dataframe. For one group (IE CO2) of runs a null dataframe is used.
#For multiple groups (IE CO2 and CO2 + Tansy) the function can be run twice with
#The output from the first run used as the input for the second run.
trt.point_extract<-function(sub_dir, trt.points_out=NULL) {
  
  #---- Track level files
  #Creates a list of "treatment level point predictions.csv" files
  trt.point.files<-
    list.files(sub_dir, pattern="treatment level point predictions.csv", recursive = TRUE)
  
  #Reads in and rbinds each csv file
  for(i in seq(1, length(trt.point.files), 1)){
    trt.points_i<-
      readr::read_csv(paste0(sub_dir,trt.point.files[i]), 
               #keeps all columns as characters
               col_types = readr::cols(.default = "c"))
    trt.points_out<-bind_rows(trt.points_out, trt.points_i)
    remove(trt.points_i)
  }
  #output
  trt.points_out
}
